//! Deterministic curve selection and order-certificate generation.
//!
//! `Family` adds witness construction to a certificate `Policy`. The
//! candidate and order conditions belong to `criteria::Criteria`.
//! A rejection label is an odd factor or a cofactor condition:
//! label 8 means order eight on a full Edwards curve, or order four in its
//! two-torsion quotient; label 2 means rational two-torsion on Weierstrass.
//! Witness points are deterministically hashed under `TAG_WITNESS` and
//! multiplied into the required subgroup. A verifier accepts any valid
//! witness, so a sieve may supply a different point for the same rejection.
//!
//! `Count` supplies the full elliptic-curve order, `Factor` resolves odd
//! composites that trial division does not settle, and `Sieve` may reject
//! candidates before counting. The optional `pari` backend supplies counting
//! and factoring. The binary families also support native AGM counting
//! through `agm::Agm` and full-order rejection
//! witnesses, avoiding a factoring dependency.
//!
//! `criteria::find` may count with early abort: a small prime outside the
//! allowed cofactor suffices to reject without obtaining the full order.
//! `prove` needs a witness for every retained candidate before acceptance,
//! so it counts in full unless the family or sieve already supplies one.
//! Structural exclusions follow the family's certificate-entry convention.

use super::criteria::OrderError;
use crate::curve::binary127;
use crate::curvegen::select::{self, Policy, is_prime};
use crate::hash::Salted;

pub use crate::curvegen::select::Binary127;

pub const TAG_WITNESS: &[u8] = b"ephemeral-ecmh/prove/witness";

/// A certificate entry: a rejection label and its encoded group witness.
pub type Rejection = (u128, [u8; 16]);

/// Factoring for the rare order that trial division doesn't settle.
pub trait Factor {
    /// The smallest prime factor of m, an odd composite with no prime
    /// factor below 2^10; or `None` to reject by order instead, in a family
    /// that can (`Family::order_witness`).
    fn smallest_prime_factor(&mut self, m: u128) -> Option<u128>;
}

/// Factors nothing: composites are rejected by order.
pub struct NoFactor;

impl Factor for NoFactor {
    fn smallest_prime_factor(&mut self, _: u128) -> Option<u128> {
        None
    }
}

/// How to witness a counted order's rejection, or certify its acceptance.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Verdict {
    /// #E = cofactor * r with r passing the probable-prime policy.
    Accept(u128),
    /// The order fails criteria for which this format has no rejection witness.
    Inadmissible(OrderError),
    /// Rejected by an odd-order or family-specific cofactor witness.
    Reject(u128),
    /// Rejected by a point of order the smallest prime factor of this odd
    /// composite (which has none below 2^10).
    Composite(u128),
}

/// A certificate policy as the prover sees it: candidates, bounds,
/// cofactor and acceptance come from `Policy`; the prover adds how to
/// witness a rejection.
pub trait Family: Policy<16> + Sized {
    /// What the full order n says about its candidate, as the policy
    /// judges r = n / cofactor: a cofactor with more 2-power than the
    /// family's is the marker label, an odd factor is a rejection, and
    /// an accepted r must pass the order criteria. MOV and anomalous
    /// refusals have no rejection witness in this format; the prover
    /// handles them separately from a composite order.
    fn verdict(n: u128) -> Verdict {
        let cofactor = Self::R.cofactor;
        assert!(
            Self::R.hasse_contains(n) && n.is_multiple_of(cofactor),
            "invalid curve order"
        );
        if n.is_multiple_of(2 * cofactor) {
            return Verdict::Reject(2 * cofactor);
        }
        match odd_verdict(n / cofactor) {
            Verdict::Accept(r) => {
                if let Err(reason) = Self::check_order(n) {
                    return Verdict::Inadmissible(reason);
                }
                Verdict::Accept(r)
            }
            v => v,
        }
    }
    /// The acceptance check on the accepted candidate, as the verifier
    /// applies it. A failed check cannot be skipped without a rejection
    /// certificate, so the prover panics on it.
    fn accept(seed: &[u8; 32], j: u32, r: u128) -> Result<Self::Curve, select::Error> {
        select::accept::<16, Self>(seed, j, r)
    }
    /// An encoded rejection witness for label l and candidate j, with
    /// full elliptic-curve order n. Its group order is family-specific.
    fn witness(c: &Self::Curve, h: &Salted, j: u32, n: u128, l: u128) -> [u8; 16];
    /// The rejection `verdict` would give, if it takes no point count.
    fn quick_reject(c: &Self::Curve) -> Option<Rejection> {
        let _ = c;
        None
    }
    /// The point of a rejection by the order n itself, if the family's
    /// verifier takes one.
    fn order_witness(c: &Self::Curve, h: &Salted, j: u32, n: u128) -> Option<[u8; 16]> {
        let _ = (c, h, j, n);
        None
    }
}

const TRIAL: u128 = 1 << 10;

/// Verdict on an odd m that has to be prime: trial division, then
/// primality exactly as `select` tests it.
fn odd_verdict(m: u128) -> Verdict {
    if let Some(l) = (3..TRIAL).step_by(2).find(|&l| m.is_multiple_of(l)) {
        return Verdict::Reject(l);
    }
    if is_prime(m) {
        Verdict::Accept(m)
    } else {
        Verdict::Composite(m)
    }
}

/// Message of the i-th hashed try for candidate j: j || i, both u32 LE.
fn msg(j: u32, i: u32) -> [u8; 8] {
    (u64::from(i) << 32 | u64::from(j)).to_le_bytes()
}

/// A point of order exactly l (a prime, or 8), from hashed points of a
/// group of order m: multiply by m with the b-part removed (b the prime
/// under l), which leaves order b^e, then by b until l kills it. Retries
/// when that leaves O, or order below 8. Plain (m/l)P would not do: when
/// E\[l\] is Z/l x Z/l and l^2 | m, it is always O.
fn torsion_point<T: Copy>(
    m: u128,
    l: u128,
    hash: impl Fn(u32) -> T,
    mul: impl Fn(&T, u128) -> T,
    is_o: impl Fn(&T) -> bool,
) -> T {
    prime_power_point(m, if l == 8 { 2 } else { l }, l, hash, mul, is_o)
}

/// Extract order `target`, a power of `base`, after removing the entire
/// base-primary component from the multiplier. The final check excludes
/// lower orders, including when the group has noncyclic torsion.
fn prime_power_point<T: Copy>(
    m: u128,
    base: u128,
    target: u128,
    hash: impl Fn(u32) -> T,
    mul: impl Fn(&T, u128) -> T,
    is_o: impl Fn(&T) -> bool,
) -> T {
    assert!(m > 0 && base >= 2 && target >= base && m.is_multiple_of(target));
    let mut power = target;
    while power.is_multiple_of(base) {
        power /= base;
    }
    assert_eq!(power, 1, "target must be a power of the prime base");
    let mut k = m;
    while k.is_multiple_of(base) {
        k /= base;
    }
    (0..=u32::MAX)
        .map(|i| mul(&hash(i), k))
        .filter(|p| !is_o(p))
        .find_map(|mut p| {
            while !is_o(&mul(&p, target)) {
                p = mul(&p, base);
            }
            (!is_o(&mul(&p, target / base))).then_some(p)
        })
        .expect("a point of the requested prime-power order exists")
}

impl Family for Binary127 {
    fn witness(c: &Self::Curve, h: &Salted, j: u32, n: u128, l: u128) -> [u8; 16] {
        // hashed points have Tr(x) = 1: the odd part, of order n/2
        torsion_point(
            n / 2,
            l,
            |i| c.hash_to_curve(h, &msg(j, i)),
            |p, k| c.to_affine(&c.mul(&c.from_affine(p), k)),
            binary127::Affine::is_identity,
        )
        .encode()
    }

    /// Any hashed point: it has Tr(x) = 1, so it is not O and its order is
    /// odd, which is all `gf2_127_rejects` asks.
    fn order_witness(c: &Self::Curve, h: &Salted, j: u32, _: u128) -> Option<[u8; 16]> {
        Some(c.hash_to_curve(h, &msg(j, 0)).encode())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    /// Z/a x Z/b, the non-cyclic case included.
    fn torsion_in_product(a: u128, b: u128, l: u128) -> (u128, u128) {
        let hash = |i: u32| {
            let h = u128::from(i) * 0x9e37_79b9_7f4a_7c15 + 0x7f4a;
            (h % a, (h >> 7) % b)
        };
        let mul = |p: &(u128, u128), k: u128| (p.0 * (k % a) % a, p.1 * (k % b) % b);
        torsion_point(a * b, l, hash, mul, |p| *p == (0, 0))
    }

    proptest! {
        #[test]
        fn torsion_points_have_order_l(
            l in prop::sample::select(vec![2u128, 3, 5, 7, 8, 11]),
            ea in 0u32..3, eb in 0u32..3, sa in 1u128..500, sb in 1u128..500,
        ) {
            let base = if l == 8 { 2 } else { l };
            let (a, b) = (l.pow(ea) * sa, l.pow(eb) * sb);
            // some element has order l: l divides a or b (for 8, 8 does)
            prop_assume!(a % l == 0 || b % l == 0);
            let p = torsion_in_product(a, b, l);
            let mul = |k: u128| (p.0 * k % a, p.1 * k % b);
            prop_assert_ne!(p, (0, 0));
            prop_assert_eq!(mul(l), (0, 0));
            prop_assert_ne!(mul(l / base), (0, 0));
        }

        #[test]
        fn odd_verdicts(m in (1u128 << 64)..(1 << 100)) {
            let m = m | 1;
            match odd_verdict(m) {
                Verdict::Reject(l) => {
                    prop_assert!(m % l == 0 && (3..l).step_by(2).all(|d| m % d != 0));
                }
                Verdict::Accept(r) => prop_assert!(r == m && is_prime(m)),
                Verdict::Inadmissible(_) => unreachable!("odd_verdict only classifies factors"),
                Verdict::Composite(c) => prop_assert!(c == m && !is_prime(m)),
            }
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(64))]
    }

    /// In a quotient by rational two-torsion, the marker 8 witnesses order
    /// 4. A target-order parameter alone cannot infer its prime base.
    #[test]
    fn quotient_witness_has_exact_order_four() {
        for m in [4u128, 8, 12, 24] {
            let p = prime_power_point(
                m,
                2,
                4,
                |i| if i == 0 { m / 2 } else { 1 },
                |p, k| p * k % m,
                |p| *p == 0,
            );
            assert_eq!(4 * p % m, 0);
            assert_ne!(2 * p % m, 0, "group order {m}");
        }
    }

    /// The prime base and target are separate for higher prime powers too.
    #[test]
    fn prime_power_witness_does_not_return_lower_order() {
        for (base, target, m) in [(2u128, 8u128, 32u128), (3, 9, 81), (5, 25, 125)] {
            let p = prime_power_point(
                m,
                base,
                target,
                |i| if i == 0 { m / base } else { 1 },
                |p, k| p * k % m,
                |p| *p == 0,
            );
            assert_eq!(target * p % m, 0);
            assert_ne!(target / base * p % m, 0);
        }
    }

    #[test]
    fn witness_messages_are_distinct() {
        assert_eq!(msg(1, 2), [1, 0, 0, 0, 2, 0, 0, 0]);
    }
}
