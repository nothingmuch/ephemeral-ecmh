//! Small-l torsion sieve: reject a candidate curve before counting points.
//!
//! l | #E(F_q) iff E(F_q) has a point of order l, iff the l-division
//! polynomial psi_l has a root x0 in F_q whose y is in F_q too (and not
//! only on the quadratic twist). The point is exactly a certificate
//! rejection (`select`), so the sieve's witness needs no further work.
//!
//! - g = gcd(psi_l, x^q - x) collects the roots in F_q (`poly::field_roots`):
//!   log2 q squarings mod psi_l, which dominate the cost (127 here).
//!   Usually g = 1. Otherwise all roots of g are tried: a root on the
//!   quadratic twist fails to decode. Binary decoding requires
//!   Tr(x) = 1 and Tr(b/x) = 0. At q = 2^127 and l = 3, Frobenius
//!   may have eigenvalues 1 and -1, so roots can occur on both twists.

//!
//! The binary sieve (`gf2_127`) returns an encoded small-order witness
//! for the certificate verifier.

//!
//! # Cost model
//!
//! Trying l costs one gcd(psi_l, x^q - x), for any curve: psi_l has degree
//! d = (l^2 - 1)/2, and 127 squarings mod psi_l take 127 d^2/2 M over
//! GF(2^127) (a table of x^(2i) mod psi_l makes squaring linear) and about
//! 127 (3/2) d^2 M over F_p. Building psi_l and the gcd itself add O(d^2).
//! Splitting the gcd, when it is not 1, works mod polynomials of degree
//! about l and is left out. Counted with `poly::tests::Counted`, and
//! converted at 0.98 / 0.51 / 330 ns per M / S / I over GF(2^127).
//! These are design-time measurements outside the benchmark runs.
//!
//! | l  | d   | GF(2^127) M | S   | ~time   |
//! |----|-----|-------------|-----|---------|
//! | 3  | 4   | 1063        | 509 | 2.0 us  |
//! | 5  | 12  | 9406        | 1.5k| 11 us   |
//! | 7  | 24  | 37.4k       | 3.0k| 39 us   |
//! | 11 | 60  | 232k        | 7.4k| 231 us  |
//! | 13 | 84  | 452k        | 10k | 449 us  |
//! | 17 | 144 | 1.31M       | 18k | 1.30 ms |
//! | 19 | 180 | 2.06M       | 22k | 2.03 ms |
//! | 23 | 264 | 4.41M       | 32k | 4.34 ms |

//!
//! Trying l, last in the sieve, saves a point count of time T on the fraction
//! Pr[l | #E] of the candidates that reach it, so it is worthwhile iff
//! cost(l) < Pr[l | #E] T. With Frobenius a random element of GL_2(F_l) of
//! determinant q, Pr[l | #E] = l/(l^2 - 1) if q = 1 mod l, else 1/(l - 1),
//! which `sage sage/sieve.sage 2000` checks on binary candidates: 0.50, 0.26,
//! 0.16, 0.09 and 0.08 for l = 3..13. No sieved odd prime divides q - 1
//! for q = 2^127.

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

use crate::curve::binary127;
use crate::curvegen::poly::{Field, Poly, field_roots, find_root};
use crate::curvegen::select;
use crate::field::gf2_127::{Gf, to_u128};
use std::collections::BTreeMap;

/// b-invariants (b2, b4, b6, b8) of a Weierstrass model: its division
/// polynomials depend on nothing else.
pub type BInvariants<F> = [F; 4];

/// y^2 + xy = x^3 + x^2 + B: a1 = a2 = 1, a6 = B.
pub fn gf2_127_invariants(c: &binary127::Curve) -> BInvariants<Gf> {
    [Gf::ONE, Gf::ZERO, Gf::ZERO, c.big_b]
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curvegen::poly::tests::{Counted, Ops, count};
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
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(48))]
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
        let near = |got: Ops, m: u64, s: u64| {
            got.i == 2 && got.m.abs_diff(m) * 100 <= m && got.s.abs_diff(s) * 100 <= s
        };
        for (l, bm, bs, _, _) in OPS {
            assert!(near(cost(&gf2_127, l), bm, bs), "binary, l = {l}");
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
    }
}
