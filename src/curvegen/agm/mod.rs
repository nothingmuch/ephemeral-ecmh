//! Point counting on E_B : y^2 + xy = x^3 + x^2 + B over GF(2^m),
//! without PARI. The counter is generic over the field's polynomial
//! basis (`Modulus`); `binary127` has a = 1, the model counted here.

//!
//! Satoh–Skjernaa–Taguchi lifting on the AGM modular equation (Gaudry's
//! MSST). For the AGM step (a, b) -> ((a + b)/2, sqrt(ab)), the Landen
//! variable x = (a - b)/(4(a + b)) satisfies Phi(x, x') = 0 with
//!
//!   Phi(X, Y) = Y - X^2 (1 + 4Y)^2,  so Y = X^2 mod 2.
//!
//! On the canonical lift of y^2 + xy = x^3 + a6 the step is Frobenius, so
//! x solves Phi(x, sigma x) = 0 with x = sqrt(a6) mod 2 (a6 itself gives the
//! conjugate curve, same trace), and Mestre's a/a' = 1 + 4x makes the unit
//! root of Frobenius the norm N(1 + 4x).
//!
//! - Lift: Newton doubles the precision of x, 1 -> 64 bits in 6 steps. Each
//!   solves sigma(d) = a d + b with 2 | a, by d <- sigma^-1(a d + b), one
//!   bit per round.
//! - Norm: N(1 + 4x) = exp(Tr log(1 + 4x)), with (1 + 4x)^(2^6) taken first
//!   to shorten the series. Both series are rescaled so that no step
//!   divides by 2: x mod 2^64 gives the unit root mod 2^66.
//! - Trace: t = N + q/N mod 2^66, and |t| <= 2 sqrt(q) < 2^65 makes that
//!   exact for every m <= 127 (`CAP`). For m >= 66, as at m = 109, 122
//!   and 127, q/N = 2^m / N vanishes mod 2^66; below, as at m = 61, N
//!   being a 2-adic unit, `trace` computes it as 2^m N^-1.
//!
//! Cost per count, m = 127 (the `operation_counts` test): 107
//! multiplications, 12 squarings and 69 Frobenius matrix products (6 sigma,
//! 63 sigma^-1), each m^2 = 16129 64-bit multiply-adds (half that for a
//! squaring): about 2.9 million 64-bit multiply-adds in all, 59% in
//! multiplications and 38% in matrix products. No elapsed time is inferred
//! from these counts.
//!
//! The plain AGM iteration (Mestre, Harley) takes a step per bit, each a
//! 2-adic square root and inverse (~13 multiplications), ~800 in all.

mod zq;

pub use zq::{CAP, Frobenius, Gf127, Modulus, Zq};

use crate::curve::{binary, binary127};
use crate::curvegen::criteria::Count;
use crate::field::gf2_127;
use crate::hash::Salted;
use std::sync::OnceLock;

/// Squarings of 1 + 4x before the log series.
const K: u32 = 6;
const LOG: [u64; 8] = log_coefs();
const EXP: [u64; 64] = exp_coefs();
const MASK66: u128 = (1 << 66) - 1;

/// floor(2 sqrt(2^127)), as in select.rs.
const HASSE: u128 = 26087635650665564424;
const TAG_CHECK: &[u8] = b"ephemeral-ecmh/agm/check";

const fn inv_odd(a: u64) -> u64 {
    // a^2 = 1 mod 8; each Newton step doubles the bits
    let mut x = a;
    let mut i = 0;
    while i < 5 {
        x = x.wrapping_mul(2u64.wrapping_sub(a.wrapping_mul(x)));
        i += 1;
    }
    x
}

/// log(1 + 2^(K+2) w) / 2^(K+2) = w sum_n c_n (-w)^(n-1), with
/// c_n = 2^((K+2)(n-1)) / n: 0 mod 2^64 from n = 9 on.
const fn log_coefs() -> [u64; 8] {
    let mut c = [0; 8];
    let mut n = 1;
    while n <= 8 {
        let v = (n as u64).trailing_zeros();
        c[n - 1] = inv_odd((n as u64) >> v) << ((K + 2) * (n as u32 - 1) - v);
        n += 1;
    }
    c
}

/// exp(4s) = 1 + 4s sum_n d_n s^(n-1), with d_n = 4^(n-1) / n!:
/// v(d_n) = n - 2 + popcount(n) >= 64 for all n > 64 (and some below).
const fn exp_coefs() -> [u64; 64] {
    let mut d = [0; 64];
    let mut odd = 1u64;
    let mut n = 1;
    while n <= 64 {
        let k = n as u64;
        odd = odd.wrapping_mul(k >> k.trailing_zeros());
        let v = 2 * (n as u32 - 1) - (n as u32 - k.count_ones());
        d[n - 1] = if v < 64 { inv_odd(odd) << v } else { 0 };
        n += 1;
    }
    d
}

/// Tables for one field; a count then allocates nothing.
pub struct Counter<F: Modulus> {
    frob: Frobenius<F>,
}

impl<F: Modulus> Default for Counter<F> {
    fn default() -> Self {
        Self::new()
    }
}

impl<F: Modulus> Counter<F> {
    pub fn new() -> Self {
        Self {
            frob: Frobenius::new(),
        }
    }

    /// The trace t of y^2 + xy = x^3 + a6 (a6 as field bits, nonzero),
    /// which has q + 1 - t points; t = 1 mod 4.
    pub fn trace(&self, a6: u128) -> i128 {
        assert!(a6 != 0, "y^2 + xy = x^3 is singular");
        let g = self.unit_root(&self.lift(a6));
        let q_g = if F::M >= 66 { 0 } else { inv_mod66(g) << F::M };
        let t = (g + q_g) & MASK66;
        t as i128 - if t >> 65 == 1 { 1 << 66 } else { 0 }
    }

    /// The x with Phi(x, sigma x) = 0 and x = a6 mod 2, mod 2^64.
    fn lift(&self, a6: u128) -> Zq<F> {
        let one = Zq::ONE;
        let mut x = Zq::from_bits(a6);
        // x is right mod 2^n
        let mut n = 1;
        while n < 64 {
            let y = self.frob.sigma(&x);
            let e = y.scale(4) + one;
            let g = x * e;
            let phi = y - g.square();
            // Phi_X = -2 g e, Phi_Y = 1 - 8 g x = 1 mod 8
            let gx8 = (g * x).scale(8);
            let phi_y = one - gx8;
            let (mut iota, mut prec) = (one + gx8, 6);
            while prec < n {
                iota = iota * (Zq::scalar(2) - phi_y * iota);
                prec *= 2;
            }
            // x + 2^n d: Phi + 2^n (Phi_X d + Phi_Y sigma d) = 0 mod 2^2n, so
            // sigma d = a d + b mod 2^n.
            let a = (g * e).scale(2) * iota;
            let b = -(phi.shr(n) * iota);
            let mut d = self.frob.sigma_inv(&b);
            for _ in 1..n {
                d = self.frob.sigma_inv(&(a * d + b));
            }
            x = x + d.shl(n);
            n *= 2;
        }
        x
    }

    /// N(1 + 4x) mod 2^66.
    fn unit_root(&self, x: &Zq<F>) -> u128 {
        // (1 + 4x)^(2^j) = 1 + 2^(j+2) w
        let mut w = *x;
        for j in 0..K {
            w = w + w.square().shl(j + 1);
        }
        let mw = -w;
        let mut acc = Zq::scalar(LOG[LOG.len() - 1]);
        for &c in LOG[..LOG.len() - 1].iter().rev() {
            acc = acc * mw + Zq::scalar(c);
        }
        // s = Tr(log(1 + 4x)) / 4
        let s = self.frob.trace(&(w * acc));
        let e = EXP
            .iter()
            .rev()
            .fold(0u64, |e, &d| e.wrapping_mul(s).wrapping_add(d));
        1 + 4 * s.wrapping_mul(e) as u128
    }
}

fn inv_mod66(g: u128) -> u128 {
    let mut x = g;
    for _ in 0..6 {
        x = x.wrapping_mul(2u128.wrapping_sub(g.wrapping_mul(x)));
    }
    x & MASK66
}

/// n = #E for a binary curve, after a sanity check on a hashed point P.
/// P has Tr(x) = Tr(a2) = 1, so it lies in 2E, the subgroup of order #E/2
/// (#E = 2 mod 4), and its order divides #E/2: (n/2) P = O if n is right.
/// The check is not a proof of the count: it passes any n for which n/2
/// is a multiple of P's order. Panics if (n/2) P != O.
fn checked<M: binary::Model>(c: &binary::Curve<M>, n: u128) -> u128 {
    let p = c.hash_to_curve(&Salted::new(TAG_CHECK, &[0; 32]), &[]);
    assert!(
        c.mul(&c.from_affine(&p), n / 2).is_identity(),
        "point count fails the order check"
    );
    n
}

/// q + 1 + t, q = 2^m, for the trace t of the a2 = 0 curve, once |t| is
/// within `hasse` = floor(2 sqrt q).
fn twist_order(m: usize, t: i128, hasse: u128) -> u128 {
    assert!(
        t.unsigned_abs() <= hasse,
        "trace outside the Hasse interval"
    );
    (1u128 << m | 1).checked_add_signed(t).unwrap()
}

impl Counter<Gf127> {
    /// #E_B. For odd m, Tr(a2 = 1) = 1 makes E_B the quadratic twist of the
    /// a2 = 0 curve: #E_B = q + 1 + t = 2 mod 4. Panics if B = 0, or if the
    /// count fails `checked`'s sanity check.
    pub fn order(&self, c: &binary127::Curve) -> u128 {
        let t = self.trace(gf2_127::to_u128(c.big_b));
        checked(c, twist_order(127, t, HASSE))
    }
}

/// The GF(2^127) counter, built on first use (~750 multiplications).
pub fn counter() -> &'static Counter<Gf127> {
    static C: OnceLock<Counter<Gf127>> = OnceLock::new();
    C.get_or_init(Counter::new)
}

/// Native point-counting backend for the binary candidates.
pub struct Agm;

impl Count<binary127::Curve> for Agm {
    fn order(&mut self, c: &binary127::Curve) -> u128 {
        counter().order(c)
    }
}

#[cfg(test)]
mod tests {
    use super::zq::ops;
    use super::*;
    use crate::curve::binary127::tests::curve_and_points;
    use crate::field::gf2_127::is_zero;
    use proptest::prelude::*;

    #[test]
    fn series_coefficients() {
        // n c_n = 2^(8(n-1)), n! d_n = 4^(n-1), mod 2^64
        for n in 1..=8u64 {
            let want = 1u64.checked_shl(8 * (n as u32 - 1)).unwrap_or(0);
            assert_eq!(LOG[n as usize - 1].wrapping_mul(n), want);
        }
        let mut f = 1u64;
        for n in 1..=64u64 {
            f = f.wrapping_mul(n);
            let want = 1u64.checked_shl(2 * (n as u32 - 1)).unwrap_or(0);
            assert_eq!(EXP[n as usize - 1].wrapping_mul(f), want);
        }
    }

    #[test]
    fn lift_is_canonical() {
        let c = counter();
        let x = c.lift(0x1234_5678_9abc_def0_1357_9bdf_2468_ace0);
        let y = c.frob.sigma(&x);
        assert_eq!(y, x.square() * (y.scale(4) + Zq::ONE).square());
    }

    #[test]
    fn operation_counts() {
        let c = counter();
        let before = ops::get();
        c.trace(0x1234_5678_9abc_def0_1357_9bdf_2468_ace0);
        let after = ops::get();
        let d: Vec<u64> = (0..3).map(|i| after[i] - before[i]).collect();
        assert_eq!(d[ops::MUL], 107);
        assert_eq!(d[ops::SQR], 12);
        assert_eq!(d[ops::FROB], 69);
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(16))]

        #[test]
        fn order_kills_points((c, ps) in curve_and_points(2)) {
            let n = counter().order(&c);
            prop_assert_eq!(n % 4, 2);
            prop_assert!(n.abs_diff(1 << 127 | 1) <= HASSE);
            for p in ps {
                prop_assert!(is_zero(c.mul(&c.from_affine(&p), n).x));
                prop_assert!(is_zero(c.mul(&c.from_affine(&p), n / 2).x));
            }
        }
    }
}
