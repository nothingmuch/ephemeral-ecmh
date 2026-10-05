//! PARI/GP point counting for `prove` (feature `pari`).
//!
//! Which algorithm PARI 2.17.3 runs:
//!
//! - GF(2^n), n = 127, 109 and 122: `ellcard` goes to `F2xq_ellcard`,
//!   which is Harley's variant of Mestre's AGM (`F2xq_elltrace_Harley`): a
//!   Newton lift of the AGM equation to 2-adic precision (n + 1)/2 + 2,
//!   then a norm. There is no SEA in characteristic 2, and `ellsea` falls
//!   back to the same code, ignoring `tors`, so binary candidates get no
//!   early abort. GF(2^122) is the Rust tower GF(2^61)\[u\]/(u^2 + u + 1)
//!   embedded in a field of PARI's choosing (`State::tower`); GLS curves,
//!   whose B lies in GF(2^61), are counted there too, over GF(2^122), not
//!   through the subfield.
//! - p = 2^127 - 1: `ellcard` tests for CM, then runs SEA
//!   (`Fp_ellcard_SEA`, from 56 bits up). SEA uses the modular polynomials
//!   of the `seadata` package when it can find them, and computes them
//!   (`polmodular`) otherwise, which is slower (see `seadata`).
//! - Quadratic odd fields: `ellcard` counts over the explicit extension
//!   in the Rust basis; `ellsea` uses SEA there with tors-based early abort.
//!
//! `ellsea(E, tors)` returns 0 once some Elkies prime l not dividing `tors`
//! divides #E, which is `Count::order_early_abort`. tors = 2 (Edwards) or
//! 1 (prime order) never aborts an order that `prove` would accept, but it
//! reveals neither l nor #E, so only `find` uses it. With tors = 2 it
//! cannot abort on 8 | #E either; `Family::quick_reject` settles that
//! before any count.
//!
//! PARI keeps its stack in thread-locals, so every call runs on one
//! dedicated thread, started on first use. A PARI error exits the process.

use crate::curvegen::criteria::Count;
use std::ffi::{CStr, CString, c_char, c_long, c_ulong};
use std::path::PathBuf;
use std::sync::OnceLock;
use std::sync::mpsc;

use crate::curve::{binary127, edwards127, weier127};
use crate::curvegen::prove::Factor;
use crate::field::fp127::P;
use crate::field::gf2_127;

type Gen = *mut c_long;

const INIT_JMPM: c_ulong = 1;
const INIT_DFTM: c_ulong = 4;
const INIT_NOIMTM: c_ulong = 16;
const D_SILENT: c_long = 0;
/// Real precision for `ellinit`, in bits; exact curves don't use it.
const DEFAULTPREC: c_long = 64;

// Linked by build.rs.
unsafe extern "C" {
    static gen_0: Gen;
    static gen_1: Gen;
    /// reassigned by `setdefault("datadir")`
    static mut pari_datadir: *const c_char;

    fn pari_init_opts(parisize: usize, maxprime: c_ulong, opts: c_ulong);
    fn paristack_setsize(rsize: usize, vsize: usize);
    fn setdefault(s: *const c_char, v: *const c_char, flag: c_long) -> Gen;
    fn get_avma() -> c_ulong;
    fn set_avma(av: c_ulong);
    fn gclone(x: Gen) -> Gen;

    fn uutoi(hi: c_ulong, lo: c_ulong) -> Gen;
    fn stoi(x: c_long) -> Gen;
    fn itou(x: Gen) -> c_ulong;
    fn shifti(x: Gen, n: c_long) -> Gen;
    fn remi2n(x: Gen, n: c_long) -> Gen;
    fn mkintmod(x: Gen, y: Gen) -> Gen;
    fn mkvec5(a: Gen, b: Gen, c: Gen, d: Gen, e: Gen) -> Gen;
    fn gmul(x: Gen, y: Gen) -> Gen;
    fn binaire(x: Gen) -> Gen;
    fn gtopoly(x: Gen, v: c_long) -> Gen;
    fn ffgen(t: Gen, v: c_long) -> Gen;
    fn Fq_to_FF(x: Gen, ff: Gen) -> Gen;

    fn ellinit(x: Gen, d: Gen, prec: c_long) -> Gen;
    fn ellcard(e: Gen, p: Gen) -> Gen;
    fn ellsea(e: Gen, tors: c_long) -> Gen;
    fn Z_factor(n: Gen) -> Gen;
}

/// The fields' parameters, cloned off the stack once.
struct State {
    p: Gen,
    /// generator of F_2[t]/(t^127 + t^63 + 1)
    g: Gen,
}

type Job = Box<dyn FnOnce(&State) + Send>;

fn int(v: u128) -> Gen {
    unsafe { uutoi((v >> 64) as c_ulong, v as c_ulong) }
}

/// The polynomial over F_2 whose coefficient of x^i is bit i of v, with
/// integer coefficients, in variable 0.
fn bits_poly(v: u128) -> Gen {
    unsafe { gtopoly(binaire(int(v)), 0) }
}

/// A cloned generator of F_2[t]/(f), f given by its bits.
fn binary_generator(f: u128) -> Gen {
    unsafe { gclone(ffgen(gmul(bits_poly(f), mkintmod(gen_1, stoi(2))), -1)) }
}

/// Non-negative t_INT below 2^128.
fn to_int(n: Gen) -> u128 {
    unsafe {
        let hi = itou(shifti(n, -64)) as u128;
        let lo = itou(remi2n(n, 64)) as u128;
        hi << 64 | lo
    }
}

/// gel(x, i): the i-th component of a vector, column or matrix.
unsafe fn gel(x: Gen, i: usize) -> Gen {
    unsafe { *x.add(i) as Gen }
}

fn start() -> mpsc::Sender<Job> {
    let (tx, rx) = mpsc::channel::<Job>();
    std::thread::Builder::new()
        .name("pari".into())
        .stack_size(64 << 20)
        .spawn(move || {
            let state = unsafe {
                // no INIT_SIGm: leave the signal handlers to Rust
                pari_init_opts(8 << 20, 1 << 20, INIT_JMPM | INIT_DFTM | INIT_NOIMTM);
                // grows on demand up to 4 GiB of address space
                paristack_setsize(8 << 20, 4 << 30);
                if std::env::var_os("GP_DATA_DIR").is_none()
                    && let Some(dir) = option_env!("GP_DATA_DIR")
                {
                    let (k, v) = (c"datadir", CString::new(dir).unwrap());
                    setdefault(k.as_ptr(), v.as_ptr(), D_SILENT);
                }
                State {
                    p: gclone(int(P)),
                    g: binary_generator(1 << 127 | 1 << 63 | 1),
                }
            };
            for job in rx {
                let av = unsafe { get_avma() };
                job(&state);
                unsafe { set_avma(av) };
            }
        })
        .expect("spawn the PARI thread");
    tx
}

/// Runs `f` on the PARI thread; its stack is reset afterwards, so nothing
/// on it may escape.
fn run<T: Send + 'static>(f: impl FnOnce(&State) -> T + Send + 'static) -> T {
    static JOBS: OnceLock<mpsc::Sender<Job>> = OnceLock::new();
    let (tx, rx) = mpsc::sync_channel(1);
    let job: Job = Box::new(move |s| tx.send(f(s)).unwrap());
    JOBS.get_or_init(start).send(job).unwrap();
    rx.recv().expect("the PARI thread died")
}

/// The seadata directory PARI will read, if it is there: without it SEA
/// still works, but computes its modular polynomials itself, which is
/// slower.
pub fn seadata() -> Option<PathBuf> {
    let dir = run(|_| unsafe { CStr::from_ptr(pari_datadir).to_string_lossy().into_owned() });
    let sea = PathBuf::from(dir).join("seadata");
    sea.join("sea0").exists().then_some(sea)
}

/// PARI's `ellcard`, `ellsea` and `factor`.
#[derive(Clone, Copy, Debug, Default)]
pub struct Pari;

/// y^2 + xy = x^3 + x^2 + B over F_2[t]/(f), g a generator, B as bits.
fn binary_ell(g: Gen, big_b: u128) -> Gen {
    unsafe {
        let b = Fq_to_FF(bits_poly(big_b), g);
        ellinit(
            mkvec5(gen_1, gen_1, gen_0, gen_0, b),
            std::ptr::null_mut(),
            DEFAULTPREC,
        )
    }
}

/// y^2 = x^3 + a2 x^2 + a4 x + a6 over F_p.
fn fp127_ell(s: &State, a2: u128, a4: u128, a6: u128) -> Gen {
    unsafe {
        let a = mkvec5(gen_0, int(a2), gen_0, int(a4), int(a6));
        ellinit(a, s.p, DEFAULTPREC)
    }
}

fn edwards127_ell(s: &State, c: &edwards127::Curve) -> Gen {
    fp127_ell(s, c.a2.value(), c.a4.value(), 0)
}

fn weier127_ell(s: &State, c: &weier127::Curve) -> Gen {
    fp127_ell(s, 0, P - 3, c.b.value())
}

fn sea(e: Gen, tors: u32) -> Option<u128> {
    let n = to_int(unsafe { ellsea(e, tors.into()) });
    (n != 0).then_some(n)
}

impl Count<binary127::Curve> for Pari {
    fn order(&mut self, c: &binary127::Curve) -> u128 {
        let big_b = gf2_127::to_u128(c.big_b);
        run(move |s| to_int(unsafe { ellcard(binary_ell(s.g, big_b), std::ptr::null_mut()) }))
    }
}

impl Count<edwards127::Curve> for Pari {
    fn order(&mut self, c: &edwards127::Curve) -> u128 {
        let c = *c;
        run(move |s| to_int(unsafe { ellcard(edwards127_ell(s, &c), std::ptr::null_mut()) }))
    }

    fn order_early_abort(&mut self, c: &edwards127::Curve, tors: u32) -> Option<u128> {
        let c = *c;
        run(move |s| sea(edwards127_ell(s, &c), tors))
    }
}

impl Count<weier127::Curve> for Pari {
    fn order(&mut self, c: &weier127::Curve) -> u128 {
        let c = *c;
        run(move |s| to_int(unsafe { ellcard(weier127_ell(s, &c), std::ptr::null_mut()) }))
    }

    fn order_early_abort(&mut self, c: &weier127::Curve, tors: u32) -> Option<u128> {
        let c = *c;
        run(move |s| sea(weier127_ell(s, &c), tors))
    }
}

impl Factor for Pari {
    fn smallest_prime_factor(&mut self, m: u128) -> Option<u128> {
        // the primes of the factorization matrix come sorted
        Some(run(move |_| {
            to_int(unsafe { gel(gel(Z_factor(int(m)), 1), 1) })
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curve::weier;

    unsafe extern "C" {}

    /// Pairs of field elements as bits under `mask`, from a fixed xorshift
    /// stream, with the all-ones element first: products that reduce.
    fn samples(mask: u128) -> Vec<(u128, u128)> {
        let mut x = 0x9e37_79b9_7f4a_7c15_f39c_c060_5ced_c834_u128;
        let mut next = move || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x & mask
        };
        let mut v = vec![(mask, mask)];
        v.extend((0..16).map(|_| (next(), next())));
        v
    }

    /// x^3 - 3x + b has a root, a point of order 2, exactly when #E is
    /// even: `OddCurve::new` against PARI's count, over F_p and GF(p^2).
    #[test]
    fn odd_curve_matches_the_parity_of_the_count() {
        let mut parities = [0; 2];
        let mut check = |odd: bool, n: u128| {
            assert_eq!(odd, n % 2 == 1, "#E = {n}");
            parities[odd as usize] += 1;
        };
        for (b, _) in samples(crate::field::fp127::P) {
            let Some(c) = weier127::Curve::new(crate::field::fp127::Fp::new(b)) else {
                continue; // the all-ones sample is p, so b = 0
            };
            check(weier::OddCurve::new(c).is_some(), Pari.order(&c));
        }
        assert!(parities.iter().all(|&k| k > 0), "{parities:?}");
    }
}
