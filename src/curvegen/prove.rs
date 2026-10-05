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

use super::agm::Agm;
use super::criteria::{Count, Criteria, Filter, Found, OrderError, find};
use crate::curve::{binary127, edwards127, weier127};
use crate::curvegen::select::{self, Certificate, Policy, affine_mul, is_prime};
use crate::hash::Salted;

mod gf2;
pub use crate::curvegen::select::{Binary127, Edwards127, Weier127};
pub use gf2::{Binary109, certificate109};

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

/// Early abort with a witness: a rejection entry that `select` accepts for
/// this candidate, or `None` to have it counted. It must never reject the
/// candidate `select` would accept.
pub trait Sieve<C> {
    fn reject(&mut self, c: &C) -> Option<Rejection>;
}

// The same witness sieve must skip exactly the same candidates in either path.
impl<C, S: Sieve<C>> Filter<C> for S {
    fn rejects(&mut self, c: &C) -> bool {
        self.reject(c).is_some()
    }
}

/// Supplies no rejection witnesses; remaining candidates need a point count.
#[derive(Clone, Copy, Debug)]
pub struct NoSieve;

impl<C> Sieve<C> for NoSieve {
    fn reject(&mut self, _: &C) -> Option<Rejection> {
        None
    }
}

/// The small-l torsion sieve (`sieve`) over odd l up to the bound; the
/// families' own 8 and 2 come first, as in `select`, unless
/// `Family::quick_reject` has settled them already.
#[derive(Clone, Copy, Debug)]
pub struct SmallL(pub u32);

impl Sieve<binary127::Curve> for SmallL {
    fn reject(&mut self, c: &binary127::Curve) -> Option<Rejection> {
        crate::curvegen::sieve::gf2_127(c, self.0)
    }
}

/// `Edwards127::quick_reject` has tried 8.
impl Sieve<edwards127::Curve> for SmallL {
    fn reject(&mut self, c: &edwards127::Curve) -> Option<Rejection> {
        crate::curvegen::sieve::edwards_odd(c, self.0)
    }
}

impl Sieve<weier127::Curve> for SmallL {
    fn reject(&mut self, c: &weier127::Curve) -> Option<Rejection> {
        crate::curvegen::sieve::weier(c, self.0)
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

impl Family for Edwards127 {
    fn witness(c: &Self::Curve, h: &Salted, j: u32, n: u128, l: u128) -> [u8; 16] {
        // complete: the Edwards law has no exceptions on any of E(F_p)
        torsion_point(
            n,
            l,
            |i| c.hash_to_curve(h, &msg(j, i)),
            |p, k| c.to_affine(&c.mul(&c.from_affine(p), k)),
            edwards127::Affine::is_identity,
        )
        .encode()
    }

    /// 8 | #E iff 1 - d is a square: two square roots (`sieve::edwards_order8`).
    fn quick_reject(c: &Self::Curve) -> Option<Rejection> {
        crate::curvegen::sieve::edwards_order8(c).map(|p| (8, p))
    }
}

impl Family for Weier127 {
    fn witness(c: &Self::Curve, h: &Salted, j: u32, n: u128, l: u128) -> [u8; 16] {
        let hash = |i| c.hash_to_curve(h, &msg(j, i));
        let is_o = weier127::Affine::is_identity;
        let p = if !n.is_multiple_of(2) {
            // no 2-torsion, so RCB is complete
            let mul = |p: &weier127::Affine, k| c.to_affine(&c.mul(&c.from_affine(p), k));
            torsion_point(n, l, hash, mul, is_o)
        } else {
            // the affine chord law is complete on every candidate
            torsion_point(n, l, hash, |p, k| affine_mul(c, p, k), is_o)
        };
        p.encode()
    }
}

/// The first order-admissible candidate and its certificate, if this
/// format can witness the preceding rejections and its acceptance point
/// passes. Failure panics rather than silently choosing a later candidate;
/// direct selection can succeed even when certification cannot.
pub fn prove<F: Family>(
    seed: &[u8; 32],
    count: &mut impl Count<F::Curve>,
    factor: &mut impl Factor,
    sieve: &mut impl Sieve<F::Curve>,
) -> (Certificate, F::Curve) {
    let h = Salted::new(TAG_WITNESS, seed);
    let mut rejections = Vec::new();
    for j in 0..=u32::MAX {
        let Some(c) = F::candidate(seed, j) else {
            if F::ENTRY_PER_INDEX {
                rejections.push((0, [0; 16]));
            }
            continue;
        };
        if let Some(w) = F::quick_reject(&c).or_else(|| sieve.reject(&c)) {
            rejections.push(w);
            continue;
        }
        let n = count.order(&c);
        let l = match F::verdict(n) {
            Verdict::Accept(r) => {
                F::accept(seed, j, r).expect("accepted candidate fails verification");
                let cert = Certificate {
                    index: j,
                    r,
                    rejections,
                };
                return (cert, c);
            }
            Verdict::Inadmissible(reason) => {
                panic!("order rejected ({reason:?}): this seed has no certificate");
            }
            Verdict::Reject(l) => l,
            Verdict::Composite(m) => match factor.smallest_prime_factor(m) {
                Some(l) => l,
                None => {
                    let p = F::order_witness(&c, &h, j, n)
                        .expect("this family rejects composites only by a factor");
                    rejections.push((n, p));
                    continue;
                }
            },
        };
        rejections.push((l, F::witness(&c, &h, j, n, l)));
    }
    unreachable!("2^32 candidates rejected")
}

/// The binary family's curve for a seed, in Rust alone: the sieve at its
/// default bound, then AGM counts.
pub fn find_gf2_127(seed: &[u8; 32]) -> Found<binary127::Curve> {
    find::<Binary127>(
        seed,
        &mut Agm,
        &mut SmallL(crate::curvegen::sieve::GF2_127_L_MAX),
    )
}

/// `find_gf2_127` with its certificate, in Rust alone: a composite r with no
/// small factor is rejected by its order, not factored.
pub fn certify_gf2_127(seed: &[u8; 32]) -> (Certificate, binary127::Curve) {
    prove::<Binary127>(
        seed,
        &mut Agm,
        &mut NoFactor,
        &mut SmallL(crate::curvegen::sieve::GF2_127_L_MAX),
    )
}

pub fn prove_gf2_127(
    seed: &[u8; 32],
    count: &mut impl Count<binary127::Curve>,
    factor: &mut impl Factor,
    sieve: &mut impl Sieve<binary127::Curve>,
) -> (Certificate, binary127::Curve) {
    prove::<Binary127>(seed, count, factor, sieve)
}

pub fn prove_fp127(
    seed: &[u8; 32],
    count: &mut impl Count<edwards127::Curve>,
    factor: &mut impl Factor,
    sieve: &mut impl Sieve<edwards127::Curve>,
) -> (Certificate, edwards127::Curve) {
    prove::<Edwards127>(seed, count, factor, sieve)
}

pub fn prove_weier127(
    seed: &[u8; 32],
    count: &mut impl Count<weier127::Curve>,
    factor: &mut impl Factor,
    sieve: &mut impl Sieve<weier127::Curve>,
) -> (Certificate, weier127::Curve) {
    prove::<Weier127>(seed, count, factor, sieve)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curvegen::criteria::AdmissibleR;
    use crate::field::fp127::Fp;
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

        /// A point of order 8 when 1 - d is a square; that there is none
        /// otherwise is checked against Sage's orders (tests/prove.rs).
        #[test]
        fn edwards127_quick_reject_has_order_8(c in edwards127::tests::curve()) {
            let got = Edwards127::quick_reject(&c);
            prop_assert_eq!(got.is_some(), (Fp::ONE - c.d).is_square());
            prop_assert_eq!(Edwards127::reject_without_count(&c), got.is_some());
            if let Some((l, enc)) = got {
                let p = c.from_affine(&c.decode(enc).unwrap());
                let o = edwards127::Point::IDENTITY;
                prop_assert_eq!(l, 8);
                prop_assert!(!c.mul(&p, 4).equals(&o) && c.mul(&p, 8).equals(&o));
            }
        }
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

    /// The prime r = 2^20 + 7, whose window is about q = r + 5. The order of
    /// 5 mod r is r - 1 (5^((r - 1) / d) != 1 for d = 2, 29, 101, 179, the
    /// prime factors of r - 1 = 2 * 29 * 101 * 179), above 2^20, so q passes
    /// the embedding-degree test, and r != q.
    const TOY_R: u128 = (1 << 20) + 7;

    /// Z/r as a curve of order `TOY_R`: a point is its residue, candidate
    /// j is the integer j, and the hashed point is 1 under a zero seed and
    /// O otherwise, so acceptance fails under any other seed.
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    struct ModR(u32);

    impl select::Arithmetic<16> for ModR {
        type Affine = u32;
        type Point = u32;
        fn decode(&self, enc: [u8; 16]) -> Option<u32> {
            Some((u128::from(enc[0]) % TOY_R) as u32)
        }
        fn hash(&self, salt: &Salted, _: &[u8]) -> u32 {
            let zero = Salted::new(select::TAG_CERT_POINT, &[0; 32]);
            u32::from(salt.digest(&[], 0) == zero.digest(&[], 0))
        }
        fn lift(&self, p: &u32) -> Option<u32> {
            (*p != 0).then_some(*p)
        }
        fn mul(&self, p: &u32, k: u128) -> u32 {
            (u128::from(*p) * k % TOY_R) as u32
        }
        fn is_identity(&self, p: &u32) -> bool {
            *p == 0
        }
    }

    struct SkippingFamily;

    impl Criteria for SkippingFamily {
        type Curve = ModR;
        const TAG: &'static [u8] = b"test";
        const R: AdmissibleR = AdmissibleR::new(TOY_R + 5, 6, 1);
        fn candidate(_: &[u8; 32], j: u32) -> Option<ModR> {
            (j != 0 && j != 2).then_some(ModR(j))
        }
    }

    impl Policy<16> for SkippingFamily {
        const CLEAR: u128 = 1;
        const ENTRY_PER_INDEX: bool = false;
        fn marker(_: &ModR, _: &u32, _: u128) -> Option<bool> {
            None
        }
    }

    impl Family for SkippingFamily {
        fn verdict(n: u128) -> Verdict {
            if n == TOY_R {
                Verdict::Accept(TOY_R)
            } else {
                Verdict::Reject(3)
            }
        }
        fn witness(c: &ModR, _: &Salted, _: u32, _: u128, _: u128) -> [u8; 16] {
            [c.0 as u8; 16]
        }
    }

    struct SkippingCount(Vec<u32>);
    impl Count<ModR> for SkippingCount {
        fn order(&mut self, c: &ModR) -> u128 {
            self.0.push(c.0);
            if c.0 == 4 { TOY_R } else { 9 }
        }
    }

    #[test]
    fn find_and_prove_share_structural_skip_accounting() {
        let seed = [0; 32];
        let mut counted = SkippingCount(vec![]);
        let found = find::<SkippingFamily>(&seed, &mut counted, &mut NoSieve);
        assert_eq!((found.index, found.r, found.curve), (4, TOY_R, ModR(4)));
        assert_eq!(counted.0, [1, 3, 4]);
        let mut counted = SkippingCount(vec![]);
        let (cert, curve) =
            prove::<SkippingFamily>(&seed, &mut counted, &mut NoFactor, &mut NoSieve);
        assert_eq!((cert.index, cert.r, curve), (4, TOY_R, ModR(4)));
        assert_eq!(counted.0, [1, 3, 4]);
        assert_eq!(cert.rejections, [(3, [1; 16]), (3, [3; 16])]);
    }

    #[test]
    fn find_accepts_an_order_even_when_the_certificate_point_fails() {
        let found = find::<SkippingFamily>(&[1; 32], &mut SkippingCount(vec![]), &mut NoSieve);
        assert_eq!((found.index, found.r), (4, TOY_R));
        assert_eq!(
            SkippingFamily::accept(&[1; 32], 4, TOY_R),
            Err(select::Error::OrderCheck)
        );
    }

    #[test]
    #[should_panic(expected = "accepted candidate fails verification")]
    fn prove_rejects_an_unverifiable_acceptance() {
        prove::<SkippingFamily>(
            &[1; 32],
            &mut SkippingCount(vec![]),
            &mut NoFactor,
            &mut NoSieve,
        );
    }

    #[test]
    fn witness_messages_are_distinct() {
        assert_eq!(msg(1, 2), [1, 0, 0, 0, 2, 0, 0, 0]);
    }
}
