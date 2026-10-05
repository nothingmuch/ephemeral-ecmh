//! Direct curve selection from trusted point counts, independent of certificates.
//!
//! The verifier uses the same order-parameter predicates, but binds r to
//! #E through a nonidentity point and the Hasse bound, conditional on r
//! actually being prime. A direct check needs neither that point nor a
//! witness for each rejected candidate.

use super::number::{embedding_degree_ok, is_prime};

/// The order parameters a family admits: r with `cofactor * r` in the
/// Hasse interval `q + 1 ± hasse`, as the range `lo..=hi` of r itself, so
/// that an untrusted r is bounded without forming `cofactor * r`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AdmissibleR {
    /// The field size.
    pub q: u128,
    /// floor(2 sqrt q).
    pub hasse: u128,
    /// #E / r for an accepted curve.
    pub cofactor: u128,
    /// ceil((q + 1 - hasse) / cofactor).
    pub lo: u128,
    /// floor((q + 1 + hasse) / cofactor).
    pub hi: u128,
}

impl AdmissibleR {
    /// `q + 1 + hasse` may exceed u128 (p = 2^128 - 275 does); the
    /// quotient is formed from the parts. The least admissible r exceeds
    /// the interval's width, so a claimed multiple of r in the interval
    /// is the only one.
    pub const fn new(q: u128, hasse: u128, cofactor: u128) -> Self {
        let q1 = q + 1;
        let lo = (q1 - hasse).div_ceil(cofactor);
        let hi = q1 / cofactor + hasse / cofactor + (q1 % cofactor + hasse % cofactor) / cofactor;
        assert!(lo > 2 * hasse);
        Self {
            q,
            hasse,
            cofactor,
            lo,
            hi,
        }
    }

    /// Check a trusted full point count; primality remains a fixed-base
    /// probable-prime test, just as it is for certificate verification.
    pub fn check_order(&self, n: u128) -> Result<u128, OrderError> {
        if !self.hasse_contains(n) {
            return Err(OrderError::Hasse);
        }
        if !n.is_multiple_of(self.cofactor) {
            return Err(OrderError::Cofactor);
        }
        let r = n / self.cofactor;
        self.check_r(r)?;
        Ok(r)
    }

    /// The shared order-parameter conditions. Certificate verification
    /// puts its point check between these phases, so an invalid point
    /// does not incur the embedding-degree computation.
    pub fn check_r(&self, r: u128) -> Result<(), OrderError> {
        self.check_r_before_embedding(r)?;
        self.check_embedding(r)
    }

    pub(crate) fn check_r_before_embedding(&self, r: u128) -> Result<(), OrderError> {
        if !is_prime(r) {
            return Err(OrderError::RNotPrime);
        }
        if !self.admits(r) {
            return Err(OrderError::Hasse);
        }
        if r == self.q {
            return Err(OrderError::Anomalous);
        }
        Ok(())
    }

    pub(crate) fn check_embedding(&self, r: u128) -> Result<(), OrderError> {
        if !embedding_degree_ok(self.q % r, r) {
            return Err(OrderError::Mov);
        }
        Ok(())
    }

    pub fn admits(&self, r: u128) -> bool {
        (self.lo..=self.hi).contains(&r)
    }

    /// Whether n is a possible full order, |q + 1 - n| <= hasse.
    pub fn hasse_contains(&self, n: u128) -> bool {
        n.abs_diff(self.q + 1) <= self.hasse
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OrderError {
    Hasse,
    Cofactor,
    RNotPrime,
    Anomalous,
    Mov,
}

/// A deterministic candidate family and the order conditions it must meet.
/// No point representation or certificate format is needed to check a count.
pub trait Criteria {
    type Curve;
    /// Domain separation for deriving candidates from the seed.
    const TAG: &'static [u8];
    /// Full curve order is cofactor times an admitted probable prime.
    const R: AdmissibleR;
    /// Primes permitted in the cofactor, for a counter's early abort.
    const TORS: u32 = if Self::R.cofactor == 1 { 1 } else { 2 };

    fn check_order(n: u128) -> Result<u128, OrderError> {
        Self::R.check_order(n)
    }

    /// A structural rejection that needs neither counting nor a witness.
    fn reject_without_count(_c: &Self::Curve) -> bool {
        false
    }

    fn candidate(seed: &[u8; 32], j: u32) -> Option<Self::Curve>;
}

/// Full elliptic-curve order of a supported candidate, fitting in u128.
/// A quotient codec does not change the meaning of this count.
pub trait Count<C> {
    fn order(&mut self, c: &C) -> u128;

    /// #E, or `None` once a prime not dividing `tors` is known to divide it
    /// (PARI's `ellsea(E, tors)`). The default always counts in full.
    fn order_early_abort(&mut self, c: &C, tors: u32) -> Option<u128> {
        let _ = tors;
        Some(self.order(c))
    }
}

/// A true result must imply failure of the family's order criteria.
/// Unlike a certificate sieve, a filter need not construct evidence.
pub trait Filter<C> {
    fn rejects(&mut self, c: &C) -> bool;
}

/// The first order-admissible candidate, without a certificate.
#[derive(Clone, Copy, Debug)]
pub struct Found<C> {
    pub index: u32,
    pub r: u128,
    pub curve: C,
}

/// Uses trusted point counts and order criteria only. This can succeed
/// even when the certificate format cannot witness an earlier rejection,
/// or its fixed acceptance point is the identity.
pub fn find<F: Criteria>(
    seed: &[u8; 32],
    count: &mut impl Count<F::Curve>,
    filter: &mut impl Filter<F::Curve>,
) -> Found<F::Curve> {
    for j in 0..=u32::MAX {
        let Some(c) = F::candidate(seed, j) else {
            continue;
        };
        if F::reject_without_count(&c) || filter.rejects(&c) {
            continue;
        }
        let Some(n) = count.order_early_abort(&c, F::TORS) else {
            continue;
        };
        if let Ok(r) = F::check_order(n) {
            return Found {
                index: j,
                r,
                curve: c,
            };
        }
    }
    unreachable!("2^32 candidates rejected")
}
