//! Small-l torsion sieve: reject a candidate curve before counting points.
//!
//! l | #E(F_q) iff E(F_q) has a point of order l, iff the l-division
//! polynomial psi_l has a root x0 in F_q whose y is in F_q too (and not
//! only on the quadratic twist). The point is exactly a certificate
//! rejection (`select`), so the sieve's witness needs no further work.
//!
//! - g = gcd(psi_l, x^q - x) collects the roots in F_q (`poly::field_roots`):
//!   log2 q squarings mod psi_l, plus products by x for set bits of q.
//!   Usually g = 1. Otherwise all roots of g are tried: a root on the
//!   quadratic twist fails to decode. Binary decoding requires
//!   Tr(x) = 1 and Tr(b/x) = 0; odd-field decoding requires a square
//!   right-hand side. Frobenius eigenvalues 1 and -1 on E\[l\] require
//!   q = -1 mod l. This occurs at q = 2^127 and l = 3.

//!
//! The binary, prime-field Edwards and Weierstrass sieves return encoded
//! small-order witnesses. The model-specific first tests are 8 | #E for
//! Edwards and 2 | #E for Weierstrass.

//!
//! # Cost model
//!
//! Trying l costs one gcd(psi_l, x^q - x), for any curve: psi_l has degree
//! d = (l^2 - 1)/2, and 127 squarings mod psi_l take 127 d^2/2 M over
//! GF(2^127) (a table of x^(2i) mod psi_l makes squaring linear) and about
//! 127 (3/2) d^2 M over F_p. Building psi_l and the gcd itself add O(d^2).
//! Splitting the gcd, when it is not 1, works mod polynomials of degree
//! about l and is left out. Counted with `poly::tests::Counted`, and
//! converted at 0.98 / 0.51 / 330 ns per M / S / I over GF(2^127) (the
//! criterion field-op baselines); over F_p, timed on an M4 (the median of
//! 5 candidates, single-shot). Every timing in this module is a
//! design-time measurement made outside the benchmark runs:
//!
//! | l  |   d | GF(2^127) M |  S  | ~time   | F_p M   |  S  | time    |
//! |----|-----|-------------|-----|---------|---------|-----|---------|
//! |  3 |   4 |        1063 | 509 | 2.0 us  |    2285 | 508 | 14 us   |
//! |  5 |  12 |        9406 | 1.5k| 11 us   |   24.9k | 1.5k| 88 us   |
//! |  7 |  24 |       37.4k | 3.0k| 39 us   |    103k | 3.0k| 0.31 ms |
//! | 11 |  60 |        232k | 7.4k| 231 us  |    654k | 7.4k| 1.82 ms |
//! | 13 |  84 |        452k |  10k| 449 us  |   1.28M |  10k| 3.53 ms |
//! | 17 | 144 |       1.31M |  18k| 1.30 ms |   3.76M |  18k| 10.3 ms |
//! | 19 | 180 |       2.06M |  22k| 2.03 ms |   5.87M |  22k| 15.9 ms |
//! | 23 | 264 |       4.41M |  32k| 4.34 ms |   12.6M |  32k| 35 ms   |
//!
//! Trying l, last in the sieve, saves a point count of time T on the fraction
//! Pr[l | #E] of the candidates that reach it, so it is worthwhile iff
//! cost(l) < Pr[l | #E] T. With Frobenius a random element of GL_2(F_l) of
//! determinant q, Pr[l | #E] = l/(l^2 - 1) if q = 1 mod l, else 1/(l - 1),
//! which sage/sieve.sage checks on sampled candidates. For q = 2^127,
//! no sieved odd prime divides q - 1; for p = 2^127 - 1, the relevant
//! primes dividing p - 1 are 3, 7 and 19.

//!
//! - GF(2^127): `agm` counts points in 0.43 ms, so l = 7 is worthwhile
//!   (71 us saved against 39, 49 measured) and 11 is not (43 us against
//!   0.23 ms): `GF2_127_L_MAX` = 7, the optimal bound for any T from 0.24
//!   to 2.2 ms. It leaves
//!   0.31 of the candidates to count, against 0.51 for 3. PARI's 3.5 ms
//!   count would call for 11 (0.35 ms saved against 0.23), which leaves
//!   0.29.
//!
//! Over F_p, Karatsuba products with Barrett reduction (not implemented in
//! `poly`) would replace the 3/2 d^2 M of a squaring mod psi_l by about
//! three d-by-d products of O(d^1.58) each. Over GF(2^127) squaring is
//! already a linear map at d^2/2 M; modular composition is an alternative
//! for computing x^(2^127).

use crate::curve::encoding::Signed;
use crate::curve::{binary127, edwards, weier};
use crate::curvegen::poly::{Field, Poly, field_roots, find_root};
use crate::curvegen::select;
use crate::field::OddField;
use crate::field::gf2_127::{Gf, to_u128};
use std::collections::BTreeMap;

/// b-invariants (b2, b4, b6, b8) of a Weierstrass model: its division
/// polynomials depend on nothing else.
pub type BInvariants<F> = [F; 4];

/// y^2 + xy = x^3 + x^2 + B: a1 = a2 = 1, a6 = B.
pub fn gf2_127_invariants(c: &binary127::Curve) -> BInvariants<Gf> {
    [Gf::ONE, Gf::ZERO, Gf::ZERO, c.big_b]
}

/// A field whose points `edwards` and `weier` encode in 16 bytes.
pub trait Signed16: Field + Signed<Bytes = [u8; 16]> {}

impl<F: Field + Signed<Bytes = [u8; 16]>> Signed16 for F {}

/// The Montgomery model y^2 = x^3 + a2 x^2 + a4 x, in which `edwards`
/// encodes points.
pub fn edwards_invariants<F: Signed16>(c: &edwards::Curve<F>) -> BInvariants<F> {
    let two = F::small(2);
    [two * two * c.a2, two * c.a4, F::ZERO, -c.a4.square()]
}

/// y^2 = x^3 - 3x + b.
pub fn weier_invariants<F: Signed16>(c: &weier::Curve<F>) -> BInvariants<F> {
    [F::ZERO, -F::small(6), F::small(4) * c.b, -F::small(9)]
}

/// psi_2^2 = 4x^3 + b2 x^2 + 2 b4 x + b6, a polynomial in x on the curve.
fn psi2_squared<F: Field>(b: &BInvariants<F>) -> Poly<F> {
    let [b2, b4, b6, _] = *b;
    Poly::new(vec![b6, F::small(2) * b4, b2, F::small(4)])
}

/// f_n = psi_n for odd n and psi_n / psi_2 for even n: polynomials in x
/// for any Weierstrass model, with psi_2 = 2y + a1 x + a3 the only factor
/// that needs y. deg f_n = (n^2 - 1)/2 for odd n, leading coefficient n.
pub fn division_polynomial<F: Field>(b: &BInvariants<F>, n: usize) -> Poly<F> {
    let [b2, b4, b6, b8] = *b;
    let k = F::small;
    let mut memo = BTreeMap::from([
        (0, Poly::zero()),
        (1, Poly::constant(F::ONE)),
        (2, Poly::constant(F::ONE)),
        (3, Poly::new(vec![b8, k(3) * b6, k(3) * b4, b2, k(3)])),
        (
            4,
            Poly::new(vec![
                b4 * b8 - b6.square(),
                b2 * b8 - b4 * b6,
                k(10) * b8,
                k(10) * b6,
                k(5) * b4,
                b2,
                k(2),
            ]),
        ),
    ]);
    let ff = psi2_squared(b).square();
    division_rec(&mut memo, &ff, n)
}

/// The recurrences for psi_n, with psi_2^2 folded into ff = psi_2^4
/// wherever psi_2 appears to an even power; only the ~5 indices around n/2
/// at each of log n levels are computed.
fn division_rec<F: Field>(memo: &mut BTreeMap<usize, Poly<F>>, ff: &Poly<F>, n: usize) -> Poly<F> {
    if let Some(f) = memo.get(&n) {
        return f.clone();
    }
    let m = n / 2;
    let mut f = |i| division_rec(memo, ff, i);
    let out = if n % 2 == 1 {
        // psi_2m+1 = psi_m+2 psi_m^3 - psi_m-1 psi_m+1^3
        let (fm2, fm, fm1, fp1) = (f(m + 2), f(m), f(m - 1), f(m + 1));
        let mut a = &fm2 * &(&fm.square() * &fm);
        let mut b = &fm1 * &(&fp1.square() * &fp1);
        if m.is_multiple_of(2) {
            a = &a * ff;
        } else {
            b = &b * ff;
        }
        &a - &b
    } else {
        // psi_2 psi_2m = psi_m (psi_m+2 psi_m-1^2 - psi_m-2 psi_m+1^2)
        let (fm, fm2, fm1, fmm2, fp1) = (f(m), f(m + 2), f(m - 1), f(m - 2), f(m + 1));
        &fm * &(&(&fm2 * &fm1.square()) - &(&fmm2 * &fp1.square()))
    };
    memo.insert(n, out.clone());
    out
}

/// The first root of m in F that `decode` accepts.
fn torsion<F: Field>(m: &Poly<F>, decode: impl Fn(F) -> Option<[u8; 16]>) -> Option<[u8; 16]> {
    find_root(&field_roots(m), &mut |x| decode(x))
}

/// An encoded P != O with l P = O (so of order l for prime l), if E(F_q)
/// has one; l odd. Odd-order points are in 2E = {Tr(x) = 1}, which bit 0
/// of the encoding holds.
pub fn gf2_127_torsion(c: &binary127::Curve, l: u32) -> Option<[u8; 16]> {
    assert!(l >= 3 && l % 2 == 1);
    let psi = division_polynomial(&gf2_127_invariants(c), l as usize);
    torsion(&psi, |x| {
        c.decode(to_u128(x).to_le_bytes())
            .filter(|p| !p.is_identity())
            .map(|p| p.encode())
    })
}

/// `gf2_127_torsion` for the Edwards families, on their Montgomery model;
/// l odd, or 8 for a point of order exactly 8 (`edwards_order8`).
pub fn edwards_torsion<F: Signed16>(c: &edwards::Curve<F>, l: u32) -> Option<[u8; 16]> {
    if l == 8 {
        return edwards_order8(c);
    }
    assert!(l >= 3 && l % 2 == 1);
    let psi = division_polynomial(&edwards_invariants(c), l as usize);
    torsion(&psi, |x: F| {
        c.decode(F::to_bytes(x.pack())).map(|p| p.encode())
    })
}

/// The Montgomery x of points of order 8, if any, on a u^2 + v^2 =
/// 1 + d u^2 v^2 with a = ±1 and d and ad non-squares: two square roots,
/// no division polynomial. y may be on the twist.
///
/// The Montgomery model is y^2 = x (x^2 + a2 x + e^2), a2 = (a + d)/2 and
/// e = (a - d)/4. Its only rational point of order 2 is (0, 0), since
/// a2^2 - 4e^2 = ad is a non-square, so the 2-part of #E is cyclic; the
/// points of order 4 are those over x = ±e, and only over e is y^2 = a e^2
/// a square (over -e it is d e^2). So the points P of order 8 are those
/// with x(2P) = e, which solves to x + e^2/x = 2(e + r) with r^2 = ae,
/// so x = r (ar + 1 + t), t^2 = 2ar + 1, for either r: no r or no t, no
/// such P. The other t gives x(P + (0, 0)), rational with P.
fn order8_x<F: Field + OddField>(a: F, d: F) -> impl Iterator<Item = F> {
    let e = (a - d) * F::small(4).inv();
    let r = (a * e).sqrt();
    r.into_iter().flat_map(|r| [r, -r]).filter_map(move |r| {
        let t = (a * (r + r) + F::ONE).sqrt()?;
        Some(r * (a * r + F::ONE + t))
    })
}

/// A point of order 8 if 8 | #E, for the rule `select` checks on the
/// Edwards side (`order8_x`, a = 1). Over F_p, p = 3 mod 4, 8 | #E iff
/// e = (1 - d)/4 is a square, by 2-descent (P is a double iff x(P) is a
/// square) on T = (e, e); sieve::tests checks it against f_8.
pub fn edwards_order8<F: Signed16>(c: &edwards::Curve<F>) -> Option<[u8; 16]> {
    order8_x(F::ONE, c.d).find_map(|x| c.decode(F::to_bytes(x.pack())).map(|p| p.encode()))
}

/// `gf2_127_torsion` for the Weierstrass families; l = 2 or odd. The
/// rational 2-torsion points are (x, 0) for the roots x in F of
/// x^3 - 3x + b.
pub fn weier_torsion<F: Signed16>(c: &weier::Curve<F>, l: u32) -> Option<[u8; 16]> {
    assert!(l == 2 || (l >= 3 && l % 2 == 1));
    let m = if l == 2 {
        Poly::new(vec![c.b, -F::small(3), F::ZERO, F::ONE])
    } else {
        division_polynomial(&weier_invariants(c), l as usize)
    };
    torsion(&m, |x: F| {
        c.decode(F::to_bytes(x.pack())).map(|p| p.encode())
    })
}

/// A certificate rejection (l, encoded P), as `select` verifies them.
pub type Rejection = (u128, [u8; 16]);

/// The odd l to sieve with over GF(2^127), from the cost model.
pub const GF2_127_L_MAX: u32 = 7;

/// The largest odd l to sieve with before each caller's point count:
/// `criteria::find`'s early abort saves less than `prove::prove`'s full count
/// (see the cost model).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Bounds {
    pub find: u32,
    pub prove: u32,
}

pub const EDWARDS127: Bounds = Bounds { find: 5, prove: 13 };
pub const WEIER127: Bounds = Bounds { find: 5, prove: 13 };

fn odd_primes(l_max: u32) -> impl Iterator<Item = u32> {
    (3..=l_max).step_by(2).filter(|&n| {
        (3..)
            .step_by(2)
            .take_while(|d| d * d <= n)
            .all(|d| n % d != 0)
    })
}

/// The smallest odd prime l <= l_max dividing #E, with a point of order l;
/// `None` if the candidate needs a point count.
pub fn gf2_127(c: &binary127::Curve, l_max: u32) -> Option<Rejection> {
    odd_primes(l_max).find_map(|l| Some((l as u128, gf2_127_torsion(c, l)?)))
}

/// `gf2_127` for the Edwards families, after 8 (as the prover, sage/kat.sage,
/// tries it first).
pub fn edwards<F: Signed16>(c: &edwards::Curve<F>, l_max: u32) -> Option<Rejection> {
    edwards_order8(c)
        .map(|p| (8, p))
        .or_else(|| edwards_odd(c, l_max))
}

/// `edwards` without the rule for 8, for a prover whose
/// `Family::quick_reject` has applied it already.
pub fn edwards_odd<F: Signed16>(c: &edwards::Curve<F>, l_max: u32) -> Option<Rejection> {
    odd_primes(l_max).find_map(|l| Some((l as u128, edwards_torsion(c, l)?)))
}

/// `gf2_127` for the Weierstrass families, after 2.
pub fn weier<F: Signed16>(c: &weier::Curve<F>, l_max: u32) -> Option<Rejection> {
    core::iter::once(2)
        .chain(odd_primes(l_max))
        .find_map(|l| Some((l as u128, weier_torsion(c, l)?)))
}

/// Sieves candidates from `start` on, pushing their rejections, and
/// returns the first one the sieve passes: the next to count points on.
/// A prover alternates this with point counting, pushing a rejection of
/// its own for each counted candidate whose r is composite.
pub fn next_gf2_127(seed: &[u8; 32], start: u32, l_max: u32, rej: &mut Vec<Rejection>) -> u32 {
    for j in start.. {
        // B = 0 (probability 2^-127) is no curve, but `verify_gf2_127` still
        // takes, and skips, an entry for it.
        let Some(c) = select::gf2_127_candidate(seed, j) else {
            rej.push((0, [0; 16]));
            continue;
        };
        match gf2_127(&c, l_max) {
            Some(r) => rej.push(r),
            None => return j,
        }
    }
    unreachable!()
}

/// `next_gf2_127` for the Edwards family: incomplete candidates take no entry.
pub fn next_fp127(seed: &[u8; 32], start: u32, l_max: u32, rej: &mut Vec<Rejection>) -> u32 {
    for j in start.. {
        let Some(c) = select::fp127_candidate(seed, j) else {
            continue;
        };
        match edwards(&c, l_max) {
            Some(r) => rej.push(r),
            None => return j,
        }
    }
    unreachable!()
}

/// `next_gf2_127` for the Weierstrass family: singular candidates take no entry.
pub fn next_weier127(seed: &[u8; 32], start: u32, l_max: u32, rej: &mut Vec<Rejection>) -> u32 {
    for j in start.. {
        let Some(c) = select::weier127_candidate(seed, j) else {
            continue;
        };
        match weier(&c, l_max) {
            Some(r) => rej.push(r),
            None => return j,
        }
    }
    unreachable!()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curve::{edwards127, weier127};
    use crate::curvegen::poly::tests::{Counted, Ops, count};
    use crate::field::fp127::Fp;
    use proptest::prelude::*;

    /// x(nP) = x - psi_n-1 psi_n+1 / psi_n^2, from the f_n.
    fn x_multiple<F: Field>(b: &BInvariants<F>, x: F, n: usize) -> F {
        let f = |i| division_polynomial(b, i).eval(x);
        let p2 = psi2_squared(b).eval(x);
        let (num, den) = if n % 2 == 1 {
            (p2 * f(n - 1) * f(n + 1), f(n).square())
        } else {
            (f(n - 1) * f(n + 1), p2 * f(n).square())
        };
        x - num * den.inv()
    }

    fn degrees<F: Field + core::fmt::Debug>(b: &BInvariants<F>) {
        for n in [3usize, 5, 7, 9, 11, 13] {
            let f = division_polynomial(b, n);
            assert_eq!(f.degree(), Some((n * n - 1) / 2));
            let lead = *f.coeffs().last().unwrap();
            assert!((lead - F::small(n as u64)).is_zero());
        }
        // even: psi_n / psi_2 has degree (n^2 - 4)/2 in odd characteristic
        if !F::CHAR2 {
            for n in [4usize, 6, 8] {
                assert_eq!(division_polynomial(b, n).degree(), Some((n * n - 4) / 2));
            }
        }
    }

    proptest! {
        #[test]
        fn gf2_127_multiples((c, ps) in binary127::tests::curve_and_points(1)) {
            let p = ps[0];
            let b = gf2_127_invariants(&c);
            let mut q = p;
            for n in 2..12 {
                q = q.add(&p);
                prop_assert!(crate::field::gf2_127::eq(x_multiple(&b, p.x, n), q.x), "n = {}", n);
            }
            degrees(&b);
        }

        #[test]
        fn edwards127_multiples((c, ps) in edwards127::tests::curve_and_points(1)) {
            let p = ps[0];
            let b = edwards_invariants(&c);
            let mut q = p;
            for n in 2..12 {
                q = c.add(&q, &p);
                prop_assert_eq!(x_multiple(&b, p.x, n), q.x, "n = {}", n);
            }
            degrees(&b);
        }

        #[test]
        fn weier127_multiples((c, ps) in weier127::tests::curve_and_points(1)) {
            let p = ps[0];
            let b = weier_invariants(&c);
            let mut q = p;
            for n in 2..12 {
                q = c.add(&q, &p);
                prop_assert_eq!(x_multiple(&b, p.x, n), q.x, "n = {}", n);
            }
            degrees(&b);
        }
    }

    // Witnesses pass `select`'s rejection checks. l = 3 and 5 divide #E
    // for about a third and a fifth of curves, so most runs see several.
    proptest! {
        #![proptest_config(ProptestConfig::with_cases(64))]

        #[test]
        fn gf2_127_witnesses_have_order_l(c in binary127::tests::curve(), l in prop::sample::select(vec![3u32, 5])) {
            if let Some(enc) = gf2_127_torsion(&c, l) {
                let p = c.decode(enc).unwrap();
                prop_assert!(!p.is_identity());
                prop_assert!(crate::field::gf2_127::is_zero(c.mul(&c.from_affine(&p), l as u128).x));
            }
        }

        #[test]
        fn edwards127_witnesses_have_order_l(c in edwards127::tests::curve(), l in prop::sample::select(vec![3u32, 5])) {
            if let Some(enc) = edwards_torsion(&c, l) {
                let p = c.decode(enc).unwrap();
                let o = edwards127::Point::IDENTITY;
                prop_assert!(!p.is_identity());
                prop_assert!(c.mul(&c.from_affine(&p), l as u128).equals(&o));
            }
        }

        #[test]
        fn weier127_witnesses_have_order_l(c in weier127::tests::curve(), l in prop::sample::select(vec![2u32, 3, 5])) {
            if let Some(enc) = weier_torsion(&c, l) {
                let p = c.decode(enc).unwrap();
                prop_assert!(!p.is_identity());
                if l == 2 {
                    prop_assert!(p.y.is_zero());
                } else {
                    prop_assert!(c.mul(&c.from_affine(&p), l as u128).is_identity());
                }
            }
        }

        /// A curve built with a root x0 of x^3 - 3x + b has a point of
        /// order 2, which the sieve must find.
        #[test]
        fn weier127_finds_a_constructed_2_torsion_point(x0 in crate::field::fp127::tests::fp()) {
            let three = crate::field::fp127::Fp::small(3);
            let Some(c) = weier127::Curve::new(three * x0 - x0 * x0.square()) else {
                return Ok(()); // b = 0, or a singular curve
            };
            let p = c.decode(weier_torsion(&c, 2).expect("a root")).unwrap();
            prop_assert!(p.y.is_zero() && c.is_on_curve(&p));
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(48))]

        /// The closed form against the roots of f_8 = psi_8 / psi_2, whose
        /// rational points of order 8 are what it claims to find.
        #[test]
        fn order8_matches_division_polynomial(c in edwards127::tests::curve()) {
            let o = edwards127::Point::IDENTITY;
            let order8 = |x: Fp| {
                let p = c.decode(x.value().to_le_bytes())?;
                let p4 = c.mul(&c.from_affine(&p), 4);
                (!p4.equals(&o)).then(|| p.encode())
            };
            let f8 = division_polynomial(&edwards_invariants(&c), 8);
            let want = torsion(&f8, order8);
            let got = edwards_order8(&c);
            prop_assert_eq!(got.is_some(), want.is_some());
            prop_assert_eq!(got.is_some(), (Fp::ONE - c.d).is_square());
            if let Some(enc) = got {
                let p = c.from_affine(&c.decode(enc).unwrap());
                prop_assert!(!c.mul(&p, 4).equals(&o) && c.mul(&p, 8).equals(&o));
            }
        }
    }

    #[test]
    fn short_weierstrass_psi3() {
        // 3x^4 + 6a x^2 + 12b x - a^2 for y^2 = x^3 + ax + b, a = -3
        let b = Fp::new(12345);
        let c = weier127::Curve::new(b).unwrap();
        let want = Poly::new(vec![
            -Fp::new(9),
            Fp::new(12) * b,
            -Fp::new(18),
            Fp::ZERO,
            Fp::new(3),
        ]);
        assert_eq!(division_polynomial(&weier_invariants(&c), 3), want);
    }

    /// Ops to try l on a candidate: build psi_l, then gcd(psi_l, x^q - x).
    /// Splitting the gcd when it is not 1 is left out: its moduli have
    /// degree about l, against d = (l^2 - 1)/2 for the Frobenius.
    fn cost<F: Field>(b: &BInvariants<F>, l: u32) -> Ops {
        let b = b.map(Counted);
        count(|| field_roots(&division_polynomial(&b, l as usize))).1
    }

    /// (l, binary M, binary S, prime M, prime S) of `cost`, each with 2 I;
    /// the table in the module doc.
    const OPS: [(u32, u64, u64, u64, u64); 8] = [
        (3, 1063, 509, 2285, 508),
        (5, 9406, 1504, 24878, 1506),
        (7, 37365, 2985, 103164, 2987),
        (11, 231546, 7387, 654159, 7390),
        (13, 451750, 10302, 1281483, 10305),
        (17, 1312266, 17549, 3762855, 17564),
        (19, 2057406, 21891, 5870303, 21896),
        (23, 4407646, 31982, 12600256, 31987),
    ];

    #[test]
    fn op_counts() {
        let seed = [7; 32];
        let gf2_127 = gf2_127_invariants(&select::gf2_127_candidate(&seed, 0).unwrap());
        let ed = (0..)
            .find_map(|j| select::fp127_candidate(&seed, j))
            .unwrap();
        let we = (0..)
            .find_map(|j| select::weier127_candidate(&seed, j))
            .unwrap();
        let near = |got: Ops, m: u64, s: u64| {
            got.i == 2 && got.m.abs_diff(m) * 100 <= m && got.s.abs_diff(s) * 100 <= s
        };
        for (l, bm, bs, pm, ps) in OPS {
            assert!(near(cost(&gf2_127, l), bm, bs), "binary, l = {l}");
            assert!(
                near(cost(&edwards_invariants(&ed), l), pm, ps),
                "Edwards, l = {l}"
            );
            assert!(
                near(cost(&weier_invariants(&we), l), pm, ps),
                "Weierstrass, l = {l}"
            );
        }
    }

    /// ns per M, S, I over GF(2^127): criterion throughput baselines.
    const GF2_127_NS: [f64; 3] = [0.98, 0.51, 330.0];

    /// The cost of each l over GF(2^127), from its op counts.
    fn gf2_127_ns() -> impl Iterator<Item = (u32, f64)> {
        OPS.iter().map(|&(l, m, s, _, _)| {
            let [m_ns, s_ns, i_ns] = GF2_127_NS;
            (l, m as f64 * m_ns + s as f64 * s_ns + 2.0 * i_ns)
        })
    }

    /// The cost of each l over F_p, p = 2^127 - 1, in us: the module doc's
    /// timings.
    const FP127_US: [(u32, f64); 8] = [
        (3, 14.0),
        (5, 88.0),
        (7, 310.0),
        (11, 1820.0),
        (13, 3530.0),
        (17, 10300.0),
        (19, 15900.0),
        (23, 35000.0),
    ];

    /// The largest L such that every odd prime l <= L is worthwhile
    /// against the t(l) a rejection by l saves: cost(l) < Pr[l | #E] t(l).
    fn best(q: u128, costs: impl IntoIterator<Item = (u32, f64)>, t: impl Fn(u32) -> f64) -> u32 {
        let mut last = 1;
        for (l, c) in costs {
            let lf = l as f64;
            let pr = if q % l as u128 == 1 {
                lf / (lf * lf - 1.0)
            } else {
                1.0 / (lf - 1.0)
            };
            if c >= pr * t(l) {
                break;
            }
            last = l;
        }
        last
    }

    /// A family's bounds from its costs per l, its full count t and its
    /// early aborts a(3), a(5), a(7), a(11), in one unit.
    fn tuned(q: u128, costs: &[(u32, f64)], t: f64, a: [f64; 4]) -> Bounds {
        let abort = |l| match l {
            3 => a[0],
            5 => a[1],
            7 => a[2],
            11 => a[3],
            _ => 0.0,
        };
        Bounds {
            find: best(q, costs.iter().copied(), abort),
            prove: best(q, costs.iter().copied(), |_| t),
        }
    }

    #[test]
    fn default_bounds_follow_from_the_model() {
        let q2 = 1u128 << 127;
        let gf2 = |t| best(q2, gf2_127_ns(), |_| t);
        // agm's 0.43 ms binary point count; 7 holds from 0.24 to 2.2 ms
        assert_eq!(gf2(0.43e6), GF2_127_L_MAX);
        assert_eq!(gf2(0.24e6), GF2_127_L_MAX);
        assert_eq!(gf2(0.23e6), 5);
        assert_eq!(gf2(2.2e6), GF2_127_L_MAX);
        // PARI's 3.5 ms count; 11 holds from 2.4 to 5.4 ms
        assert_eq!(gf2(2.4e6), 11);
        assert_eq!(gf2(3.5e6), 11);
        assert_eq!(gf2(5.5e6), 13);
        // PARI's counts and early aborts over F_p, in us
        let p = (1u128 << 127) - 1;
        let edwards = tuned(p, &FP127_US, 54e3, [3.2e3, 740.0, 1.5e3, 3.2e3]);
        assert_eq!(edwards, EDWARDS127);
        let weier = tuned(p, &FP127_US, 56e3, [320.0, 860.0, 1.7e3, 1.8e3]);
        assert_eq!(weier, WEIER127);
    }
}
