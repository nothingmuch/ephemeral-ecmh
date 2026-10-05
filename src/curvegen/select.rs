//! Seed-derived curve selection, verified without point counting.
//!
//! [`Criteria`] names a family's candidates and order conditions; [`Policy`]
//! adds certificate arithmetic and witness conventions. A
//! certificate names the selected index, an order parameter r, and one
//! rejection witness for each earlier candidate. [`verify`], [`accept`] and
//! [`rejects`] are generic over the policy, and the curve arithmetic they
//! need is the [`Arithmetic`] trait, implemented once per curve model.
//!
//! # Acceptance
//!
//! The candidate at `index` must exist, and r must pass the fixed-base
//! probable-prime test [`is_prime`], with `cofactor * r` in the Hasse
//! interval ([`AdmissibleR`]), and r = q is refused (anomalous). Only a
//! cofactor-1 window holds q: 2q and 4q lie outside the Hasse interval, so
//! the refusal is vacuous elsewhere. A point P is hashed to the curve under
//! [`TAG_CERT_POINT`] and the seed, with the index as message;
//! Q = `CLEAR * P` must be nonidentity with rQ = O. The embedding degree
//! of r must exceed [`EMBEDDING_MIN`].
//!
//! If r is prime, Q has order r, so r divides #E, and `cofactor * r` is the
//! only multiple of r in the Hasse interval (its lower end exceeds its
//! width), so #E = `cofactor * r`. The argument is conditional on that
//! primality: an accepted certificate shows r to be a probable prime.
//!
//! # Rejection witnesses
//!
//! A witness is a label l and an encoded point. The point must decode to a
//! nonidentity P. Under the odd rule, l is odd, 3 <= l < the least
//! admissible r, and lP = O: P then has odd order m > 1 dividing l, so
//! m <= l is below every admissible r, while any curve of order
//! `cofactor * r` with r prime has odd part r; the candidate's order is not
//! of that form. The label need not be prime, nor P's exact order. Each
//! family admits one further label, its marker, through [`Policy::marker`]:
//!
//! - Edwards and quotient families, l = 8: `CLEAR * P` is nonidentity and
//!   `2 * CLEAR * P` = O, so the represented point has order 2·CLEAR and
//!   8 divides #E, which `4 * r` with r odd does not.
//! - Weierstrass, l = 2: y = 0, so P has order 2 and #E is even.
//! - Binary, an even l: l lies in the Hasse interval, l/2 fails `is_prime`
//!   and lP = O. Decoded binary points have odd order (Tr(x) = 1 places
//!   them in 2E), so if #E = 2r with r prime in the admitted range, P has
//!   order r, the only multiple of r in the interval is 2r, and 2r/2 passes
//!   `is_prime`: no such curve has a witness. The branch checks a witness,
//!   not that l = #E.
//!
//! # Entries
//!
//! Binary verifiers spend one entry per earlier index and ignore the entry
//! at an index that is no candidate (`ENTRY_PER_INDEX`); the odd families
//! spend entries only at indices their curve constructors admit. Missing
//! or extra entries are rejected. Verdicts come in this order:
//! `RejectionCount` for a binary certificate whose entry count is not its
//! index; `BadRejection(j)` for the first failing witness; `RejectionCount`
//! for an entry shortfall or surplus; then `RNotPrime`, `Hasse`,
//! `Anomalous`, `OrderCheck` (also for a missing accepted candidate) and
//! `Mov`. The rejection formats do not encode every acceptance failure.
//!
//! # Weil descent
//!
//! No descent check is made for GF(2^127). The basic GHS descent to F_2
//! reaches 36 Frobenius classes of curves, none of order 2r with r prime
//! (sage/ghs.sage). Hess's generalization (EUROCRYPT 2003, Section 5.1)
//! reaches 4537 classes, 65 of them of admissible order: a heuristic 2^-52
//! of the admissible isogeny classes, each with a genus 127 or 128 descent.
//! The weakness is accepted: under this heuristic the probability is about
//! 2^-52 per epoch, and by the union bound at most T * 2^-52 across T
//! epochs; a descent compromises only the epoch whose curve it reaches.
//! [`super::select109`] and [`super::select122`] discuss GF(2^109) and
//! GF(2^122).

use super::criteria::{AdmissibleR, Criteria, OrderError};
pub use super::number::{EMBEDDING_MIN, embedding_degree_ok, is_prime};

use crate::curve::binary::{self, Model};
use crate::curve::encoding::Signed;
use crate::curve::{binary127, edwards, edwards127, weier, weier127};
use crate::field::OddField;
use crate::field::fp127::{Fp, P};
use crate::field::gf2_127::{MASK127, from_u128};
use crate::hash::{Salted, halves};

pub const TAG_GF2_127: &[u8] = b"ephemeral-ecmh/curve/gf2";
pub const TAG_FP127: &[u8] = b"ephemeral-ecmh/curve/fp127";
pub const TAG_WEIER127: &[u8] = b"ephemeral-ecmh/curve/weier127";
pub const TAG_CERT_POINT: &[u8] = b"ephemeral-ecmh/cert-point";

/// floor(2 sqrt q), the same for q = 2^127 and p = 2^127 - 1.
const HASSE: u128 = 26087635650665564424;

/// `N` is the point encoding width in bytes.
#[derive(Clone, Debug)]
pub struct Certificate<const N: usize = 16> {
    pub index: u32,
    pub r: u128,
    /// Rejection labels and encoded points for candidates before `index`,
    /// one per earlier index or per earlier candidate (`Policy::ENTRY_PER_INDEX`).
    pub rejections: Vec<(u128, [u8; N])>,
}

#[derive(Debug, PartialEq, Eq)]
pub enum Error {
    RejectionCount,
    BadRejection(u32),
    Order(OrderError),
    OrderCheck,
}

/// What a verifier does on a curve, in the group its certificate speaks
/// of: E, or E/⟨T⟩ for the twisted codecs. `Affine` is what encodings and
/// the hash produce; `Point` is what scalars act on.
pub trait Arithmetic<const N: usize> {
    type Affine: Copy;
    type Point;

    fn decode(&self, enc: [u8; N]) -> Option<Self::Affine>;
    fn hash(&self, salt: &Salted, msg: &[u8]) -> Self::Affine;
    /// None for the identity, which the certificate never admits.
    fn lift(&self, p: &Self::Affine) -> Option<Self::Point>;
    fn mul(&self, p: &Self::Point, k: u128) -> Self::Point;
    fn is_identity(&self, p: &Self::Point) -> bool;
}

/// Certificate arithmetic and witness conventions for a candidate family.
///
/// A certificate names the selected index j, an order parameter r and one
/// rejection witness per earlier candidate (per earlier index under
/// `ENTRY_PER_INDEX`). Verified, it shows that candidate j exists; that r
/// passes the probable-prime test, is admitted by `R` and differs from q;
/// that a nonidentity point of the curve has order dividing r; and that
/// the embedding degree of r exceeds [`EMBEDDING_MIN`]. If r is prime,
/// #E = `R.cofactor` r follows. Each witness shows an earlier candidate's
/// order not to be `R.cofactor` times an admissible prime. The argument
/// is [the module doc](self)'s; the constants here are what varies
/// between families.
pub trait Policy<const N: usize>: Criteria<Curve: Arithmetic<N>> {
    /// Takes the hashed point into the subgroup of order r: the cofactor
    /// as the represented group sees it (1 where decoded points already
    /// have odd order).
    const CLEAR: u128;
    /// One entry per earlier index, ignored at a non-candidate, rather
    /// than one per earlier candidate. Each is a certificate format, fixed
    /// by the family's known-answer certificates; one scheme for all would
    /// change the format.
    const ENTRY_PER_INDEX: bool;

    /// `Some(verdict)` when l is the family's marker label, `None` to
    /// apply the odd rule.
    fn marker(c: &Self::Curve, p: &Affine<N, Self>, l: u128) -> Option<bool>;
}

pub type Affine<const N: usize, P> = <<P as Criteria>::Curve as Arithmetic<N>>::Affine;

/// The first 128-bit half of the tagged digest of j: the raw material of
/// every candidate.
pub(crate) fn half(tag: &[u8], seed: &[u8; 32], j: u32) -> u128 {
    halves(&Salted::new(tag, seed).digest(&[], j))[0]
}

pub(crate) fn cert_salt(seed: &[u8; 32]) -> Salted {
    Salted::new(TAG_CERT_POINT, seed)
}

pub fn verify<const N: usize, P: Policy<N>>(
    seed: &[u8; 32],
    cert: &Certificate<N>,
) -> Result<P::Curve, Error> {
    if P::ENTRY_PER_INDEX && cert.rejections.len() != cert.index as usize {
        return Err(Error::RejectionCount);
    }
    let mut entries = cert.rejections.iter();
    for j in 0..cert.index {
        let Some(c) = P::candidate(seed, j) else {
            if P::ENTRY_PER_INDEX {
                entries.next();
            }
            continue;
        };
        let &(l, enc) = entries.next().ok_or(Error::RejectionCount)?;
        if !rejects::<N, P>(&c, l, enc) {
            return Err(Error::BadRejection(j));
        }
    }
    if entries.next().is_some() {
        return Err(Error::RejectionCount);
    }
    accept::<N, P>(seed, cert.index, cert.r)
}

/// Whether (l, enc) is a rejection witness for c under the policy.
pub fn rejects<const N: usize, P: Policy<N>>(c: &P::Curve, l: u128, enc: [u8; N]) -> bool {
    let Some(p) = c.decode(enc) else {
        return false;
    };
    let Some(q) = c.lift(&p) else {
        return false;
    };
    match P::marker(c, &p, l) {
        Some(verdict) => verdict,
        None => l & 1 == 1 && (3..P::R.lo).contains(&l) && c.is_identity(&c.mul(&q, l)),
    }
}

/// The acceptance policy at `index`, independently of earlier candidates.
pub fn accept<const N: usize, P: Policy<N>>(
    seed: &[u8; 32],
    index: u32,
    r: u128,
) -> Result<P::Curve, Error> {
    let c = P::candidate(seed, index).ok_or(Error::OrderCheck)?;
    P::R.check_r_before_embedding(r).map_err(Error::Order)?;
    let p = c.hash(&cert_salt(seed), &index.to_le_bytes());
    let q = c.lift(&p).map(|p| c.mul(&p, P::CLEAR));
    match q {
        Some(q) if !c.is_identity(&q) && c.is_identity(&c.mul(&q, r)) => {}
        _ => return Err(Error::OrderCheck),
    }
    P::R.check_embedding(r).map_err(Error::Order)?;
    Ok(c)
}

/// The binary marker: an even l in the Hasse interval whose half fails
/// `is_prime`, annihilating P.
pub(crate) fn even_order<const N: usize, C: Arithmetic<N>>(
    r: &AdmissibleR,
    c: &C,
    p: &C::Affine,
    l: u128,
) -> Option<bool> {
    (l & 1 == 0).then(|| {
        r.hasse_contains(l)
            && !is_prime(l / 2)
            && c.lift(p).is_some_and(|q| c.is_identity(&c.mul(&q, l)))
    })
}

/// The Edwards marker: l = 8 with `clear * P` of order 2 in the
/// represented group.
pub(crate) fn eight_torsion<const N: usize, C: Arithmetic<N>>(
    clear: u128,
    c: &C,
    p: &C::Affine,
    l: u128,
) -> Option<bool> {
    (l == 8).then(|| {
        c.lift(p).is_some_and(|q| {
            let q = c.mul(&q, clear);
            !c.is_identity(&q) && c.is_identity(&c.mul(&q, 2))
        })
    })
}

/// The Weierstrass marker: l = 2 with y = 0. The decoded P is nonidentity.
pub(crate) fn two_torsion<F: OddField>(p: &weier::Affine<F>, l: u128) -> Option<bool> {
    (l == 2).then(|| p.y.is_zero())
}

impl<const N: usize, M: Model<Bytes = [u8; N]>> Arithmetic<N> for binary::Curve<M> {
    type Affine = binary::Affine<M>;
    type Point = binary::Point<M>;

    fn decode(&self, enc: [u8; N]) -> Option<Self::Affine> {
        binary::Curve::decode(self, enc)
    }
    fn hash(&self, salt: &Salted, msg: &[u8]) -> Self::Affine {
        self.hash_to_curve(salt, msg)
    }
    fn lift(&self, p: &Self::Affine) -> Option<Self::Point> {
        (!p.is_identity()).then(|| self.from_affine(p))
    }
    fn mul(&self, p: &Self::Point, k: u128) -> Self::Point {
        binary::Curve::mul(self, p, k)
    }
    fn is_identity(&self, p: &Self::Point) -> bool {
        p.is_identity()
    }
}

impl<const N: usize, F: OddField + Signed<Bytes = [u8; N]>> Arithmetic<N> for edwards::Curve<F> {
    type Affine = edwards::Affine<F>;
    type Point = edwards::Point<F>;

    fn decode(&self, enc: [u8; N]) -> Option<Self::Affine> {
        edwards::Curve::decode(self, enc)
    }
    fn hash(&self, salt: &Salted, msg: &[u8]) -> Self::Affine {
        self.hash_to_curve(salt, msg)
    }
    fn lift(&self, p: &Self::Affine) -> Option<Self::Point> {
        (!p.is_identity()).then(|| self.from_affine(p))
    }
    fn mul(&self, p: &Self::Point, k: u128) -> Self::Point {
        edwards::Curve::mul(self, p, k)
    }
    fn is_identity(&self, p: &Self::Point) -> bool {
        p.is_identity()
    }
}

/// Weierstrass verification arithmetic: the complete projective formulas
/// (Renes–Costello–Batina), at every encoding width. The formulas are
/// complete on a subgroup of odd order; on an exceptional pair (one whose
/// difference has order 2, as O and a 2-torsion point) they give
/// (0 : 0 : 0), which is absorbing and which the strict `is_identity`
/// refuses. A ladder's result is thus the exact point or (0 : 0 : 0), and
/// a verdict of identity is reached only when every intermediate lies in
/// a subgroup of odd order and the exact result is O. The verifier asks
/// for the identity of lP with l odd and of rQ with r odd, and
/// for the nonidentity of Q: each is decided as the chord law decides it
/// (`tests::projective_verdicts_match_the_chord_law`).
impl<const N: usize, F: OddField + Signed<Bytes = [u8; N]>> Arithmetic<N> for weier::Curve<F> {
    type Affine = weier::Affine<F>;
    type Point = weier::Point<F>;

    fn decode(&self, enc: [u8; N]) -> Option<Self::Affine> {
        weier::Curve::decode(self, enc)
    }
    fn hash(&self, salt: &Salted, msg: &[u8]) -> Self::Affine {
        self.hash_to_curve(salt, msg)
    }
    fn lift(&self, p: &Self::Affine) -> Option<Self::Point> {
        (!p.is_identity()).then(|| self.from_affine(p))
    }
    fn mul(&self, p: &Self::Point, k: u128) -> Self::Point {
        weier::Curve::mul(self, p, k)
    }
    fn is_identity(&self, p: &Self::Point) -> bool {
        p.is_identity()
    }
}

/// Affine double-and-add by the chord law.
pub(crate) fn affine_mul<F: OddField>(
    c: &weier::Curve<F>,
    p: &weier::Affine<F>,
    k: u128,
) -> weier::Affine<F> {
    (0..128 - k.leading_zeros())
        .rev()
        .fold(weier::Affine::IDENTITY, |acc, i| {
            let acc = c.add(&acc, &acc);
            if k >> i & 1 == 1 { c.add(&acc, p) } else { acc }
        })
}

/// y^2 + xy = x^3 + x^2 + B over GF(2^127), #E = 2r.
pub struct Binary127;

impl Criteria for Binary127 {
    type Curve = binary127::Curve;
    const TAG: &'static [u8] = TAG_GF2_127;
    const R: AdmissibleR = AdmissibleR::new(1 << 127, HASSE, 2);
    fn candidate(seed: &[u8; 32], j: u32) -> Option<Self::Curve> {
        let v = half(Self::TAG, seed, j) & MASK127;
        (v != 0).then(|| binary127::Curve::new(from_u128(v)))
    }
}

impl Policy<16> for Binary127 {
    const CLEAR: u128 = 1;
    const ENTRY_PER_INDEX: bool = true;

    fn marker(c: &Self::Curve, p: &binary127::Affine, l: u128) -> Option<bool> {
        even_order(&Self::R, c, p, l)
    }
}

/// Edwards x^2 + y^2 = 1 + d x^2 y^2 over F_p, p = 2^127 - 1, #E = 4r.
pub struct Edwards127;

impl Criteria for Edwards127 {
    type Curve = edwards127::Curve;
    const TAG: &'static [u8] = TAG_FP127;
    const R: AdmissibleR = AdmissibleR::new(P, HASSE, 4);
    // 8 | #E iff 1 - d is a square; accepted Edwards orders are 4r.
    fn reject_without_count(c: &Self::Curve) -> bool {
        (Fp::ONE - c.d).is_square()
    }

    fn candidate(seed: &[u8; 32], j: u32) -> Option<Self::Curve> {
        edwards127::Curve::new(Fp::new(half(Self::TAG, seed, j) & MASK127))
    }
}

impl Policy<16> for Edwards127 {
    const CLEAR: u128 = 4;
    const ENTRY_PER_INDEX: bool = false;

    fn marker(c: &Self::Curve, p: &edwards127::Affine, l: u128) -> Option<bool> {
        eight_torsion(Self::CLEAR, c, p, l)
    }
}

/// y^2 = x^3 - 3x + b over F_p, p = 2^127 - 1, #E = r.
pub struct Weier127;

impl Criteria for Weier127 {
    type Curve = weier127::Curve;
    const TAG: &'static [u8] = TAG_WEIER127;
    const R: AdmissibleR = AdmissibleR::new(P, HASSE, 1);
    fn candidate(seed: &[u8; 32], j: u32) -> Option<Self::Curve> {
        weier127::Curve::new(Fp::new(half(Self::TAG, seed, j) & MASK127))
    }
}

impl Policy<16> for Weier127 {
    const CLEAR: u128 = 1;
    const ENTRY_PER_INDEX: bool = false;

    fn marker(_: &Self::Curve, p: &weier127::Affine, l: u128) -> Option<bool> {
        two_torsion(p, l)
    }
}

pub fn gf2_127_candidate(seed: &[u8; 32], j: u32) -> Option<binary127::Curve> {
    Binary127::candidate(seed, j)
}

pub fn fp127_candidate(seed: &[u8; 32], j: u32) -> Option<edwards127::Curve> {
    Edwards127::candidate(seed, j)
}

pub fn weier127_candidate(seed: &[u8; 32], j: u32) -> Option<weier127::Curve> {
    Weier127::candidate(seed, j)
}

pub fn verify_gf2_127(seed: &[u8; 32], cert: &Certificate) -> Result<binary127::Curve, Error> {
    verify::<16, Binary127>(seed, cert)
}

pub fn gf2_127_rejects(c: &binary127::Curve, l: u128, enc: [u8; 16]) -> bool {
    rejects::<16, Binary127>(c, l, enc)
}

pub fn accept_gf2_127(seed: &[u8; 32], index: u32, r: u128) -> Result<binary127::Curve, Error> {
    accept::<16, Binary127>(seed, index, r)
}

pub fn verify_fp127(seed: &[u8; 32], cert: &Certificate) -> Result<edwards127::Curve, Error> {
    verify::<16, Edwards127>(seed, cert)
}

pub fn accept_fp127(seed: &[u8; 32], index: u32, r: u128) -> Result<edwards127::Curve, Error> {
    accept::<16, Edwards127>(seed, index, r)
}

pub fn verify_weier127(seed: &[u8; 32], cert: &Certificate) -> Result<weier127::Curve, Error> {
    verify::<16, Weier127>(seed, cert)
}

pub fn accept_weier127(seed: &[u8; 32], index: u32, r: u128) -> Result<weier127::Curve, Error> {
    accept::<16, Weier127>(seed, index, r)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curvegen::number::mul_wide;
    use crate::curvegen::sieve;
    use proptest::prelude::*;

    /// hasse = floor(2 sqrt q): hasse^2 <= 4q < (hasse + 1)^2, in 256 bits.
    pub(crate) fn hasse_is_exact(r: &AdmissibleR) {
        let wide = |a: u128| {
            let (lo, hi) = mul_wide(a, a);
            (hi, lo)
        };
        let four_q = (r.q >> 126, r.q << 2);
        assert!(wide(r.hasse) <= four_q, "{r:?}");
        assert!(four_q < wide(r.hasse + 1), "{r:?}");
    }

    #[test]
    fn admissible_r_matches_the_direct_bounds() {
        for r in [Binary127::R, Edwards127::R, Weier127::R] {
            hasse_is_exact(&r);
            assert!(
                r.hasse_contains(r.cofactor * r.lo) && !r.hasse_contains(r.cofactor * (r.lo - 1))
            );
            assert!(
                r.hasse_contains(r.cofactor * r.hi) && !r.hasse_contains(r.cofactor * (r.hi + 1))
            );
        }
        // Exact lower endpoints of the admissible 127-bit intervals.
        // The binary bound's floor is even and one below the least admissible r;
        // it cannot be an odd rejection label.
        assert_eq!(Binary127::R.lo, 85070591730234615852799834032609270652 + 1);
        assert_eq!(Edwards127::R.lo, 42535295865117307926399917016304635326);
        assert_eq!(Weier127::R.lo, P + 1 - HASSE);
    }

    proptest! {
        #[test]
        fn verifiers_reject_any_accept_r(r in any::<u128>(), hi in (1u128 << 127)..) {
            let seed = [7u8; 32];
            for r in [r, hi] {
                let cert = Certificate { index: 0, r, rejections: vec![] };
                prop_assert!(verify_gf2_127(&seed, &cert).is_err());
                prop_assert!(verify_fp127(&seed, &cert).is_err());
                prop_assert!(verify_weier127(&seed, &cert).is_err());
            }
        }

        /// The extended Edwards identity test agrees with `equals(&IDENTITY)`.
        #[test]
        fn edwards_identity_test_matches_equals(seed in any::<[u8; 32]>(), k in 1u128..64) {
            let c = (0..).find_map(|j| fp127_candidate(&seed, j)).unwrap();
            let p = c.from_affine(&c.hash_to_curve(&cert_salt(&seed), b"identity"));
            for k in [k, 4 * k, 8 * k] {
                let q = c.mul(&p, k);
                prop_assert_eq!(q.is_identity(), q.equals(&edwards127::Point::IDENTITY));
            }
        }
    }

    #[test]
    fn gf2_127_rejects_prime_r_above_2_127() {
        let r = ((1u128 << 127) + 1..)
            .step_by(2)
            .find(|&n| is_prime(n))
            .unwrap();
        let cert = Certificate {
            index: 0,
            r,
            rejections: vec![],
        };
        assert_eq!(
            verify_gf2_127(&[7; 32], &cert).unwrap_err(),
            Error::Order(OrderError::Hasse)
        );
    }

    /// A one-point "curve" whose scalar multiplication reaches O only if
    /// KILLS, and then only for scalars above 1, so that clearing keeps
    /// the point and the order check passes or fails as chosen.
    #[derive(Clone, Copy)]
    struct Toy<const KILLS: bool>;

    impl<const KILLS: bool> Arithmetic<16> for Toy<KILLS> {
        type Affine = u8;
        type Point = u8;
        fn decode(&self, _: [u8; 16]) -> Option<u8> {
            Some(1)
        }
        fn hash(&self, _: &Salted, _: &[u8]) -> u8 {
            1
        }
        fn lift(&self, p: &u8) -> Option<u8> {
            (*p != 0).then_some(*p)
        }
        fn mul(&self, p: &u8, k: u128) -> u8 {
            if KILLS && k > 1 { 0 } else { *p }
        }
        fn is_identity(&self, p: &u8) -> bool {
            *p == 0
        }
    }

    /// q = 2 with r = 3 the only admissible order: 2^2 = 1 mod 3, so the
    /// embedding degree is 2 and every acceptance fails the MOV check.
    struct ToyPolicy<const KILLS: bool>;

    impl<const KILLS: bool> Criteria for ToyPolicy<KILLS> {
        type Curve = Toy<KILLS>;
        const TAG: &'static [u8] = b"test";
        const R: AdmissibleR = AdmissibleR::new(2, 0, 1);
        fn candidate(_: &[u8; 32], _: u32) -> Option<Toy<KILLS>> {
            Some(Toy)
        }
    }

    impl<const KILLS: bool> Policy<16> for ToyPolicy<KILLS> {
        const CLEAR: u128 = 1;
        const ENTRY_PER_INDEX: bool = false;
        fn marker(_: &Toy<KILLS>, _: &u8, _: u128) -> Option<bool> {
            None
        }
    }

    /// The verdict order on a certificate failing two checks: the hashed
    /// point's order before the embedding degree.
    #[test]
    fn order_check_precedes_the_embedding_degree() {
        assert_eq!(
            accept::<16, ToyPolicy<false>>(&[0; 32], 0, 3).err(),
            Some(Error::OrderCheck)
        );
        assert_eq!(
            accept::<16, ToyPolicy<true>>(&[0; 32], 0, 3).err(),
            Some(Error::Order(OrderError::Mov))
        );
        assert_eq!(
            accept::<16, ToyPolicy<true>>(&[0; 32], 0, 4).err(),
            Some(Error::Order(OrderError::RNotPrime))
        );
    }

    /// The odd rule's labels: odd, from 3 up to but excluding the least
    /// admissible r. A toy policy with q = Q and Hasse width H at cofactor
    /// 1 admits r in Q + 1 - H..=Q + 1 + H.
    struct ToyBound<const Q: u128, const H: u128 = 0>;

    impl<const Q: u128, const H: u128> Criteria for ToyBound<Q, H> {
        type Curve = Toy<true>;
        const TAG: &'static [u8] = b"test";
        const R: AdmissibleR = AdmissibleR::new(Q, H, 1);
        fn candidate(_: &[u8; 32], _: u32) -> Option<Toy<true>> {
            Some(Toy)
        }
    }

    impl<const Q: u128, const H: u128> Policy<16> for ToyBound<Q, H> {
        const CLEAR: u128 = 1;
        const ENTRY_PER_INDEX: bool = false;
        fn marker(_: &Toy<true>, _: &u8, _: u128) -> Option<bool> {
            None
        }
    }

    /// At an odd bound, 5, and at an even one, 6, which admits the odd
    /// label immediately below it (Edwards107's case).
    #[test]
    fn odd_labels_stop_below_the_least_admissible_r() {
        assert_eq!(ToyBound::<4>::R.lo, 5);
        let admitted = |l| rejects::<16, ToyBound<4>>(&Toy, l, [0; 16]);
        assert!(admitted(3));
        assert!(!admitted(1) && !admitted(4) && !admitted(5) && !admitted(7));

        assert_eq!(ToyBound::<5>::R.lo, 6);
        let admitted = |l| rejects::<16, ToyBound<5>>(&Toy, l, [0; 16]);
        assert!(admitted(3) && admitted(5));
        assert!(!admitted(1) && !admitted(6) && !admitted(7));
    }

    /// r = q is refused where the window holds it, after the Hasse check
    /// and before the hashed point is examined.
    #[test]
    fn r_equal_to_q_is_anomalous() {
        assert_eq!((ToyBound::<5, 1>::R.lo, ToyBound::<5, 1>::R.hi), (5, 7));
        let at = |r| accept::<16, ToyBound<5, 1>>(&[0; 32], 0, r).err();
        assert_eq!(at(5), Some(Error::Order(OrderError::Anomalous)));
        assert_eq!(at(11), Some(Error::Order(OrderError::Hasse)));
        assert_eq!(at(9), Some(Error::Order(OrderError::RNotPrime)));
        // the other prime in the window is not anomalous
        assert_ne!(at(7), Some(Error::Order(OrderError::Anomalous)));
    }

    /// On candidates with rational 2-torsion, the projective arithmetic
    /// decides as the chord law does: lP = O for odd l, and the acceptance
    /// of Q = `clear` P (nonidentity, rQ = O) for the policies' clearings
    /// and odd r. Where it gives a point, it is the exact one.
    #[test]
    fn projective_verdicts_match_the_chord_law() {
        fn check<F: OddField + sieve::Signed16>(
            candidate: impl Fn(u32) -> Option<weier::Curve<F>>,
            seed: &[u8; 32],
        ) -> usize {
            let (c, t) = (0..)
                .filter_map(candidate)
                .find_map(|c| sieve::weier_torsion(&c, 2).map(|t| (c, c.decode(t).unwrap())))
                .expect("a candidate with 2-torsion");
            assert!(t.y.is_zero() && !t.is_identity());
            let p = c.hash_to_curve(&Salted::new(b"test", seed), b"p".as_slice());
            let odd = [
                1,
                3,
                5,
                7,
                9,
                15,
                17,
                255,
                257,
                (1 << 20) + 1,
                (1 << 61) - 1,
            ];
            let degenerate = |p: &weier::Point<F>| p.x.is_zero() && p.y.is_zero() && p.z.is_zero();
            let mul =
                |p: &weier::Point<F>, k: u128| <weier::Curve<F> as Arithmetic<16>>::mul(&c, p, k);
            let is_identity =
                |p: &weier::Point<F>| <weier::Curve<F> as Arithmetic<16>>::is_identity(&c, p);
            let mut degenerates = 0;
            for base in [p, c.add(&p, &t), t] {
                let lifted = <weier::Curve<F> as Arithmetic<16>>::lift(&c, &base).unwrap();
                for l in odd {
                    let exact = affine_mul(&c, &base, l);
                    let proj = mul(&lifted, l);
                    assert_eq!(is_identity(&proj), exact.is_identity(), "{l}");
                    if degenerate(&proj) {
                        degenerates += 1;
                    } else {
                        assert_eq!(c.to_affine(&proj), exact, "{l}");
                    }
                }
                for clear in [1, 2, 4] {
                    let (q, exact_q) = (mul(&lifted, clear), affine_mul(&c, &base, clear));
                    for r in odd {
                        let accepted = !is_identity(&q) && is_identity(&mul(&q, r));
                        let exact =
                            !exact_q.is_identity() && affine_mul(&c, &exact_q, r).is_identity();
                        assert_eq!(accepted, exact, "{clear} {r}");
                    }
                }
            }
            degenerates
        }
        let seed = &[3; 32];
        let n = check(|j| Weier127::candidate(seed, j), seed);
        assert!(n > 0, "no (0 : 0 : 0) was produced");
    }
}
