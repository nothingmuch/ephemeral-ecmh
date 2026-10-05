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

use crate::curve::{
    binary109, binary122, binary127, edwards, edwards127, twisted, weier, weier127,
};
use crate::curvegen::prove::Factor;
use crate::field::fp61x2::{self, Fq};
use crate::field::fp127::P;
use crate::field::gf2_122::gf2_61::{self, MASK61};
use crate::field::{OddField, fp64x2, gf2_109, gf2_122, gf2_127, goldilocks2};

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
    fn mkvec2(a: Gen, b: Gen) -> Gen;
    fn mkvec3(a: Gen, b: Gen, c: Gen) -> Gen;
    fn mkvec5(a: Gen, b: Gen, c: Gen, d: Gen, e: Gen) -> Gen;
    fn gmul(x: Gen, y: Gen) -> Gen;
    fn gadd(x: Gen, y: Gen) -> Gen;
    fn binaire(x: Gen) -> Gen;
    fn gtopoly(x: Gen, v: c_long) -> Gen;
    fn poleval(x: Gen, y: Gen) -> Gen;
    fn fetch_user_var(s: *const c_char) -> c_long;
    fn ffinit(p: Gen, n: c_long, v: c_long) -> Gen;
    fn ffgen(t: Gen, v: c_long) -> Gen;
    fn Fq_to_FF(x: Gen, ff: Gen) -> Gen;
    fn FFX_roots(f: Gen, ff: Gen) -> Gen;

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
    /// generator of F_2[t]/(t^109 + t^5 + t^4 + t^2 + 1)
    g109: Gen,
    /// z and u of the GF(2^122) tower, as roots of z^61 + z^23 + z^15 +
    /// z^5 + 1 and u^2 + u + 1 in a degree-122 field from `ffinit`. Any
    /// pair of roots embeds the tower; which one PARI returns first does
    /// not change a count.
    tower: [Gen; 2],
    /// generator i of `F_p[i]/(i² + 1)`, p = 2^61 - 1
    fp61x2: Gen,
    /// generator i of `F_p[i]/(i² - 2)`, p = 2^64 - 59
    fp64x2: Gen,
    /// generator i of `F_p[i]/(i² - 7)`, p = 2^64 - 2^32 + 1
    goldilocks2: Gen,
}

/// Exact quadratic bases supported by the PARI bridge. All three fields'
/// Hasse upper endpoints fit u128, as required by `Count`.
trait Quadratic: OddField + Send + 'static {
    const P: u64;
    /// i², as a canonical base-field coefficient.
    const NU: u64;
    fn generator(s: &State) -> Gen;
    fn coefficients(self) -> (u64, u64);
}

macro_rules! quadratic {
    ($m:ident, $nu:expr) => {
        impl Quadratic for $m::Fq {
            const P: u64 = $m::P;
            const NU: u64 = $nu;
            fn generator(s: &State) -> Gen {
                s.$m
            }
            fn coefficients(self) -> (u64, u64) {
                (self.a.value(), self.b.value())
            }
        }
    };
}

quadratic!(fp61x2, fp61x2::P - 1);
quadratic!(fp64x2, 2);

impl Quadratic for goldilocks2::Fq {
    const P: u64 = goldilocks2::P;
    const NU: u64 = 7;
    fn generator(s: &State) -> Gen {
        s.goldilocks2
    }
    fn coefficients(self) -> (u64, u64) {
        use p3_field::PrimeField64;
        let (a, b) = goldilocks2::parts(self);
        (a.as_canonical_u64(), b.as_canonical_u64())
    }
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

/// Cloned roots z of the base modulus and u of u^2 + u + 1 in GF(2^122).
fn tower_generators() -> [Gen; 2] {
    unsafe {
        // a variable below x in priority, as FFX_roots needs for its field
        let t = fetch_user_var(c"t".as_ptr());
        let k = ffgen(ffinit(stoi(2), 122, t), -1);
        let f61 = 1 << 61 | 1 << 23 | 1 << 15 | 1 << 5 | 1;
        [f61, 0b111].map(|f| gclone(gel(FFX_roots(bits_poly(f), k), 1)))
    }
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

/// A cloned generator of `F_p[i]/(i² - NU)`, initialized on the worker.
fn quadratic_generator<F: Quadratic>() -> Gen {
    unsafe {
        let polynomial = gtopoly(mkvec3(gen_1, gen_0, int((F::P - F::NU).into())), 0);
        let domain = gmul(polynomial, mkintmod(gen_1, int(F::P.into())));
        gclone(ffgen(domain, -1))
    }
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
                    g109: binary_generator(1 << 109 | 1 << 5 | 1 << 4 | 1 << 2 | 1),
                    tower: tower_generators(),
                    fp61x2: quadratic_generator::<fp61x2::Fq>(),
                    fp64x2: quadratic_generator::<fp64x2::Fq>(),
                    goldilocks2: quadratic_generator::<goldilocks2::Fq>(),
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

/// a0 + a1 u, from `gf2_122::to_u128`'s a0 | a1 << 64, in PARI's field.
fn tower_value(s: &State, v: u128) -> Gen {
    let [z, u] = s.tower;
    let (a0, a1) = (v & MASK61 as u128, v >> 64);
    unsafe {
        gadd(
            poleval(bits_poly(a0), z),
            gmul(poleval(bits_poly(a1), z), u),
        )
    }
}

/// y^2 + xy = x^3 + u x^2 + B over GF(2^122), B as `gf2_122::to_u128`.
fn binary122_ell(s: &State, big_b: u128) -> Gen {
    let [_, u] = s.tower;
    unsafe {
        let b = tower_value(s, big_b);
        ellinit(
            mkvec5(gen_1, u, gen_0, gen_0, b),
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

/// The canonical coefficients a + b i, in the same basis as Rust.
fn quadratic_value<F: Quadratic>(s: &State, x: F) -> Gen {
    unsafe {
        let (a, b) = x.coefficients();
        let coefficients = mkvec2(int(b.into()), int(a.into()));
        Fq_to_FF(gtopoly(coefficients, 0), F::generator(s))
    }
}

/// y² = x³ + a2 x² + a4 x + a6, with an explicit extension domain even
/// when every coefficient lies in the base field.
fn quadratic_ell<F: Quadratic>(s: &State, a2: F, a4: F, a6: F) -> Gen {
    unsafe {
        let a = mkvec5(
            gen_0,
            quadratic_value(s, a2),
            gen_0,
            quadratic_value(s, a4),
            quadratic_value(s, a6),
        );
        ellinit(a, F::generator(s), DEFAULTPREC)
    }
}

fn edwards61x2_ell(s: &State, c: &edwards::Curve<Fq>) -> Gen {
    quadratic_ell(s, c.a2, c.a4, Fq::ZERO)
}

fn weier61x2_ell(s: &State, c: &weier::Curve<Fq>) -> Gen {
    quadratic_ell(s, Fq::ZERO, -(Fq::ONE + Fq::ONE + Fq::ONE), c.b)
}

/// The a=-1 Edwards model is birational to y²=x³+a2*x²+a4*x,
/// a2=(d-1)/2 and a4=(d+1)²/16. Count E, not its encoded quotient.
fn twisted_ell<F: Quadratic>(s: &State, c: &twisted::Curve<F>) -> Gen {
    let half = (F::ONE + F::ONE).inv();
    let a2 = (c.d - F::ONE) * half;
    let a4 = ((c.d + F::ONE) * half.square()).square();
    quadratic_ell(s, a2, a4, F::ZERO)
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

impl Count<binary109::Curve> for Pari {
    fn order(&mut self, c: &binary109::Curve) -> u128 {
        let big_b = gf2_109::to_u128(c.big_b);
        run(move |s| to_int(unsafe { ellcard(binary_ell(s.g109, big_b), std::ptr::null_mut()) }))
    }
}

impl Count<binary122::Dense> for Pari {
    fn order(&mut self, c: &binary122::Dense) -> u128 {
        let big_b = gf2_122::to_u128(c.big_b);
        run(move |s| to_int(unsafe { ellcard(binary122_ell(s, big_b), std::ptr::null_mut()) }))
    }
}

impl Count<binary122::Gls> for Pari {
    fn order(&mut self, c: &binary122::Gls) -> u128 {
        let big_b = gf2_61::to_u64(c.big_b).into();
        run(move |s| to_int(unsafe { ellcard(binary122_ell(s, big_b), std::ptr::null_mut()) }))
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

impl Count<edwards::Curve<Fq>> for Pari {
    fn order(&mut self, c: &edwards::Curve<Fq>) -> u128 {
        let c = *c;
        run(move |s| to_int(unsafe { ellcard(edwards61x2_ell(s, &c), std::ptr::null_mut()) }))
    }

    fn order_early_abort(&mut self, c: &edwards::Curve<Fq>, tors: u32) -> Option<u128> {
        let c = *c;
        run(move |s| sea(edwards61x2_ell(s, &c), tors))
    }
}

impl Count<weier::Curve<Fq>> for Pari {
    fn order(&mut self, c: &weier::Curve<Fq>) -> u128 {
        let c = *c;
        run(move |s| to_int(unsafe { ellcard(weier61x2_ell(s, &c), std::ptr::null_mut()) }))
    }

    fn order_early_abort(&mut self, c: &weier::Curve<Fq>, tors: u32) -> Option<u128> {
        let c = *c;
        run(move |s| sea(weier61x2_ell(s, &c), tors))
    }
}

impl<F: Quadratic> Count<twisted::Curve<F>> for Pari {
    fn order(&mut self, c: &twisted::Curve<F>) -> u128 {
        let c = *c;
        run(move |s| to_int(unsafe { ellcard(twisted_ell(s, &c), std::ptr::null_mut()) }))
    }

    fn order_early_abort(&mut self, c: &twisted::Curve<F>, tors: u32) -> Option<u128> {
        let c = *c;
        run(move |s| sea(twisted_ell(s, &c), tors))
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
    use crate::field::Packed;
    use fp61x2::{Fp, P};

    unsafe extern "C" {
        fn FF_to_FpXQ(x: Gen) -> Gen;
        fn polcoef(x: Gen, degree: c_long, variable: c_long) -> Gen;
        fn gequal(x: Gen, y: Gen) -> std::ffi::c_int;
    }

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

    /// `binary_ell`'s map of GF(2^109) bits into PARI's field multiplies
    /// as `gf2_109` does: the same modulus, the same bit order.
    #[test]
    fn gf2_109_embeds() {
        let f = |v| gf2_109::from_u128(v);
        for (a, b) in samples(gf2_109::MASK109) {
            let ab = gf2_109::to_u128(f(a) * f(b));
            let same = run(move |s| unsafe {
                let x = |v| Fq_to_FF(bits_poly(v), s.g109);
                gequal(gmul(x(a), x(b)), x(ab)) != 0
            });
            assert!(same, "{a:#x} * {b:#x}");
        }
    }

    /// `tower_value` is a ring homomorphism from `gf2_122` into PARI's
    /// field, so `binary122_ell` counts the curve the tower defines.
    #[test]
    fn gf2_122_tower_embeds() {
        let f = |v| gf2_122::from_u128(v);
        let mut pairs = samples(gf2_122::MASK122);
        // u * u = u + 1, and z * z^60 = z^61 reduces in the base field
        pairs.extend([(1 << 64, 1 << 64), (2, 1 << 60)]);
        for (a, b) in pairs {
            let ab = gf2_122::to_u128(f(a) * f(b));
            let same = run(move |s| unsafe {
                let x = |v| tower_value(s, v);
                gequal(gmul(x(a), x(b)), x(ab)) != 0
            });
            assert!(same, "{a:#x} * {b:#x}");
        }
    }

    fn fq(a: u64, b: u64) -> Fq {
        Fq::new(Fp::new(a), Fp::new(b))
    }

    fn parts(x: Gen) -> (u64, u64) {
        unsafe {
            let polynomial = FF_to_FpXQ(x);
            (
                itou(polcoef(polynomial, 0, 0)),
                itou(polcoef(polynomial, 1, 0)),
            )
        }
    }

    #[test]
    fn fp61x2_basis_and_coefficients() {
        let values = [fq(0, 0), fq(1, 0), fq(0, 1), fq(P - 1, P - 2), fq(2, 3)];
        let got = run(move |s| values.map(|x| parts(quadratic_value(s, x))));
        assert_eq!(got, values.map(|x| (x.a.value(), x.b.value())));

        let (x, y) = (fq(2, 3), fq(5, 7));
        let product =
            run(move |s| unsafe { parts(gmul(quadratic_value(s, x), quadratic_value(s, y))) });
        assert_eq!(product, (P - 11, 29)); // (2+3i)(5+7i), i² = −1
    }

    // PARI 2.17.3 fixtures over i²+1, p=2^61−1. These test the binding;
    // agreement with GP is not an independent point-counting algorithm.
    // In GP, i=ffgen(Mod(1,2^61-1)*(x^2+1)), d=4+i. Edwards counts use
    // ellcard(ellinit([0,(1+d)/2,0,(1-d)^2/16,0],i)), with d or -d;
    // Weierstrass counts use ellcard(ellinit([-3,b],i)), with b=1+i or 1.
    const EDWARDS_PLUS: u128 = 5316911983139663486184127572655626632;
    const EDWARDS_MINUS: u128 = 5316911983139663488214545196698422472;
    const WEIER: u128 = 5316911983139663485373691078963717768;

    #[test]
    fn fp61x2_edwards_counts_match_gp() {
        let d = fq(4, 1);
        assert_eq!(Pari.order(&edwards::Curve::new(d).unwrap()), EDWARDS_PLUS);
        // a=1, parameter −d is isomorphic to a=−1, parameter d via u'=iu.
        assert_eq!(Pari.order(&edwards::Curve::new(-d).unwrap()), EDWARDS_MINUS);
    }

    #[test]
    fn fp61x2_weier_count_matches_gp() {
        assert_eq!(Pari.order(&weier::Curve::new(fq(1, 1)).unwrap()), WEIER);
    }

    #[test]
    fn fp61x2_base_coefficients_keep_the_extension_domain() {
        let base_order = 2305843011750340688_u128;
        let p1 = P as u128 + 1;
        let trace = p1.abs_diff(base_order);
        let expected = p1 * p1 - trace * trace;
        assert_eq!(expected, 5316911983139663485180651577861924608);
        assert_eq!(Pari.order(&weier::Curve::new(fq(1, 0)).unwrap()), expected);
    }

    #[test]
    fn fp61x2_sea_full_and_early_abort() {
        let edwards = edwards::Curve::new(fq(4, 1)).unwrap();
        let weier = weier::Curve::new(fq(1, 1)).unwrap();
        assert_eq!(Pari.order_early_abort(&edwards, 0), Some(EDWARDS_PLUS));
        assert_eq!(Pari.order_early_abort(&weier, 0), Some(WEIER));
        assert_eq!(Pari.order_early_abort(&edwards, 2), None);
        assert_eq!(Pari.order_early_abort(&weier, 1), None);
    }

    /// Check the exact basis, including canonical exports at the modulus
    /// boundary. Equal point counts alone cannot detect conjugating i.
    fn quadratic_basis<F: Quadratic + Packed>(values: [F; 5], p: u64, nu: u64) {
        let got = run(move |s| values.map(|x| parts(quadratic_value(s, x))));
        assert_eq!(got, [(0, 0), (1, 0), (0, 1), (p - 1, p - 2), (2, 3)]);

        let width = F::BITS / 2;
        let (x, y) = (values[4], F::unpack(5 | 7 << width).unwrap());
        let product =
            run(move |s| unsafe { parts(gmul(quadratic_value(s, x), quadratic_value(s, y))) });
        assert_eq!(product, (((10 + 21 * nu as u128) % p as u128) as u64, 29));
    }

    #[test]
    fn quadratic_fp61_basis() {
        quadratic_basis(
            [fq(0, 0), fq(1, 0), fq(0, 1), fq(P - 1, P - 2), fq(2, 3)],
            P,
            P - 1,
        );
    }

    #[test]
    fn quadratic_fp64_basis() {
        use fp64x2::{Fp, Fq, P};
        let fq = |a, b| Fq::new(Fp::new(a), Fp::new(b));
        quadratic_basis(
            [fq(0, 0), fq(1, 0), fq(0, 1), fq(P - 1, P - 2), fq(2, 3)],
            P,
            2,
        );
    }

    #[test]
    fn quadratic_goldilocks_basis() {
        use goldilocks2::{Fp, P, new};
        let fq = |a, b| new(Fp::new(a), Fp::new(b));
        quadratic_basis(
            [fq(0, 0), fq(1, 0), fq(0, 1), fq(P - 1, P - 2), fq(2, 3)],
            P,
            7,
        );
        // The backend also permits redundant base-field representatives.
        let x = fq(P, P + 1);
        assert_eq!(run(move |s| parts(quadratic_value(s, x))), (0, 1));
    }

    /// The returned order is that of E, while the public point operations
    /// and identity test are in G = E/⟨T⟩, of order #E/2.
    fn quotient_order<F: Quadratic + Packed>(d: F, p: u64) -> u128 {
        let c = twisted::Curve::new(d).unwrap();
        let n = Pari.order(&c);
        let p = p as u128;
        assert!(((p - 1).pow(2)..=(p + 1).pow(2)).contains(&n));
        assert_eq!(n % 4, 0);
        let salt = crate::hash::Salted::new(b"quadratic-count-binding", &[0; 32]);
        for index in 0u32..3 {
            let a = c.hash_to_curve(&salt, &index.to_le_bytes());
            assert!(c.mul(&c.from_affine(&a), n / 2).is_identity());
        }
        assert_eq!(Pari.order_early_abort(&c, 0), Some(n));
        n
    }

    #[test]
    fn quadratic_twisted61_order() {
        let d = fq(4, 1);
        assert_eq!(quotient_order(d, P), EDWARDS_MINUS);
        assert_eq!(
            Pari.order_early_abort(&twisted::Curve::new(d).unwrap(), 2),
            None
        );
    }

    #[test]
    fn quadratic_twisted64_order() {
        use fp64x2::{Fp, Fq, P};
        let d = Fq::new(Fp::ZERO, Fp::ONE);
        // GP 2.17.3, i=ffgen(Mod(1,2^64-59)*(x^2-2)), d=i:
        // ellcard(ellinit([0,(d-1)/2,0,(d+1)^2/16,0],i)).
        let n = 340282366920938461281225221704229487896;
        assert_eq!(quotient_order(d, P), n);
        // tors=2 tolerates the entire 2-part, including 8 | #E.
        assert_eq!(n % 16, 8);
        assert_eq!(
            Pari.order_early_abort(&twisted::Curve::new(d).unwrap(), 2),
            Some(n)
        );
    }

    #[test]
    fn quadratic_twisted_goldilocks_order() {
        use goldilocks2::{Fp, P, new};
        let d = new(Fp::new(0), Fp::new(1));
        // The same GP expression with i²=7 and p=2^64-2^32+1.
        assert_eq!(
            quotient_order(d, P),
            340282366762482138486535933895052498568
        );
        assert_eq!(
            Pari.order_early_abort(&twisted::Curve::new(d).unwrap(), 2),
            None
        );
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
        for (b0, b1) in samples(u128::from(P)) {
            let Some(c) = weier::Curve::new(fq(b0 as u64, b1 as u64)) else {
                continue;
            };
            check(weier::OddCurve::new(c).is_some(), Pari.order(&c));
        }
        assert!(parities.iter().all(|&k| k > 0), "{parities:?}");
    }
}
