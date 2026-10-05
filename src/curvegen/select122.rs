//! The two binary policies over GF(2^122): `select`'s scheme with a tag,
//! a candidate stream of their own.
//!
//! # Search orders
//!
//! Candidate j of a seed is derived from the first 128-bit half of
//! `Salted::new(tag, seed).digest(&[], j)`, under a tag per family, so the
//! two families never share a candidate stream:
//! - dense ([`TAG_DENSE`]): B = the half's bits under `MASK122`, as
//!   `gf2_122::from_u128` reads them (a0 = bits 0..61, a1 = bits 64..125).
//!   B in GF(4) (a0, a1 in {0, 1}) is no candidate (see below).
//! - GLS ([`TAG_GLS`]): beta = the half's low 61 bits, as
//!   `gf2_61::from_u64` reads them, and B = beta^4. beta in F_2 is no
//!   candidate.
//!
//! A j that is no candidate still takes a certificate entry, which is
//! ignored (probability 2^-120 and 2^-60).
//!
//! # What differs from GF(2^127)
//!
//! Nothing in the checks: the policy constants, and `select`'s acceptance,
//! whose embedding-degree bound the KAT curves exceed (their degrees are
//! about r itself; sage/kat122.sage prints them). GLS adds no check:
//! #E = (q - 1)^2 + t^2 over q = 2^61 (the twist of the subfield curve
//! y^2 + xy = x^3 + beta^4) is what the count gives, and psi, the GLS
//! endomorphism, reduces the rho work factor by sqrt 2 (2^59.8 against the
//! dense family's 2^60.3), which is a property of the family, not a
//! per-curve check.
//!
//! # Subfield curves and Weil descent (sage/ghs122.sage)
//!
//! GF(2^122)'s subfields are GF(2^61), GF(4) and F_2. B in GF(4) makes
//! E a curve over GF(4) (a = u is there too), a Koblitz curve whose
//! Frobenius reduces the rho work factor by sqrt(61); B = 1 has
//! #E = 2 * prime.
//! Those B are no candidates, nor is a GLS beta in F_2 (B = beta^4 = 1).
//!
//! The GHS magic number m for a descent to GF(2^l) is the F2-dimension of
//! the orbit of sqrt(B) under x -> x^(2^l), plus one. x^122 - 1 =
//! (x + 1)^2 Phi_61(x)^2 over F_2, Phi_61 irreducible (ord_61(2) = 60),
//! so:
//! - to F_2 (n = 122) and GF(4) (n = 61): m <= 2 exactly for sqrt(B) in
//!   GF(4), where the genus (<= 2, or 4 allowing Hess's extra dimension
//!   for a = u) leaves a Jacobian too small to hold E\[r\]; else
//!   m >= 60, genus >= 2^58. Dense and GLS alike (a GLS sqrt(B) outside
//!   F_2 has an orbit of dimension 60), so no curve is GHS-weak here.
//! - to GF(2^61) (n = 2): m <= 2 (1 for GLS), genus <= 2 (4) for every
//!   curve: index calculus over GF(2^61) costs about q = 2^61 for genus
//!   2 (more above), no better than rho. No per-curve check.
//!
//! The KAT certificates come from sage/kat122.sage (binary109's from
//! sage/kat109.sage); the `Dense122` and `Gls122` provers reproduce them
//! byte for byte.

use crate::curve::binary::Affine;
use crate::curve::binary122::{Dense, Gls, M122, M122Gls};
use crate::curvegen::criteria::{AdmissibleR, Criteria};
use crate::curvegen::select::{Certificate, Error, Policy, even_order, half};
use crate::field::gf2_122::gf2_61::{self, MASK61};
use crate::field::gf2_122::{MASK122, from_u128};

pub const TAG_DENSE: &[u8] = b"ephemeral-ecmh/curve/gf2-122";
pub const TAG_GLS: &[u8] = b"ephemeral-ecmh/curve/gf2-122-gls";

const Q: u128 = 1 << 122;
/// floor(2 sqrt q) = 2^62 exactly
const HASSE: u128 = 1 << 62;
const _: () = assert!(HASSE * HASSE == 4 * Q);

const R: AdmissibleR = AdmissibleR::new(Q, HASSE, 2);

/// y^2 + xy = x^3 + u x^2 + B over GF(2^122) with B dense, #E = 2r.
pub struct Dense122;

impl Criteria for Dense122 {
    type Curve = Dense;
    const TAG: &'static [u8] = TAG_DENSE;
    const R: AdmissibleR = R;
    fn candidate(seed: &[u8; 32], j: u32) -> Option<Dense> {
        let v = half(Self::TAG, seed, j) & MASK122;
        (v & !(1 | 1 << 64) != 0).then(|| Dense::new(from_u128(v)))
    }
}

impl Policy<16> for Dense122 {
    const CLEAR: u128 = 1;
    const ENTRY_PER_INDEX: bool = true;

    fn marker(c: &Dense, p: &Affine<M122>, l: u128) -> Option<bool> {
        even_order(&Self::R, c, p, l)
    }
}

/// The GLS family: B = beta^4 for beta in GF(2^61) \ F_2, #E = 2r.
pub struct Gls122;

impl Criteria for Gls122 {
    type Curve = Gls;
    const TAG: &'static [u8] = TAG_GLS;
    const R: AdmissibleR = R;
    fn candidate(seed: &[u8; 32], j: u32) -> Option<Gls> {
        let v = half(Self::TAG, seed, j) as u64 & MASK61;
        (v > 1).then(|| Gls::new(gf2_61::from_u64(v).xsquare(2)))
    }
}

impl Policy<16> for Gls122 {
    const CLEAR: u128 = 1;
    const ENTRY_PER_INDEX: bool = true;

    fn marker(c: &Gls, p: &Affine<M122Gls>, l: u128) -> Option<bool> {
        even_order(&Self::R, c, p, l)
    }
}

pub fn dense_candidate(seed: &[u8; 32], j: u32) -> Option<Dense> {
    Dense122::candidate(seed, j)
}

pub fn gls_candidate(seed: &[u8; 32], j: u32) -> Option<Gls> {
    Gls122::candidate(seed, j)
}

pub fn verify_dense(seed: &[u8; 32], cert: &Certificate) -> Result<Dense, Error> {
    crate::curvegen::select::verify::<16, Dense122>(seed, cert)
}

pub fn verify_gls(seed: &[u8; 32], cert: &Certificate) -> Result<Gls, Error> {
    crate::curvegen::select::verify::<16, Gls122>(seed, cert)
}

pub fn accept_dense(seed: &[u8; 32], index: u32, r: u128) -> Result<Dense, Error> {
    crate::curvegen::select::accept::<16, Dense122>(seed, index, r)
}

pub fn accept_gls(seed: &[u8; 32], index: u32, r: u128) -> Result<Gls, Error> {
    crate::curvegen::select::accept::<16, Gls122>(seed, index, r)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hash::{Salted, halves};
    use proptest::prelude::*;

    #[test]
    fn bounds_match_the_direct_constants() {
        // q + 1 - HASSE is odd, so the ceiling of its half is the floor plus one.
        assert_eq!(R.lo, (Q + 1 - HASSE) / 2 + 1);
        // The floor is even, so the odd rejection labels end below it.
        const { assert!(((Q + 1 - HASSE) / 2) & 1 == 0) };
        assert_eq!(R.hi, (Q + 1 + HASSE) / 2);
        assert!(R.admits(Q / 2) && !R.admits(1 << 126));
        // r = q is out of the window: 2q exceeds q + 1 + HASSE
        const { assert!(R.hi < Q) };
    }

    proptest! {
        #[test]
        fn candidates_follow_the_documented_order(seed in any::<[u8; 32]>(), j in 0u32..1000) {
            let h = |tag| halves(&Salted::new(tag, &seed).digest(&[], j))[0];
            if let Some(c) = dense_candidate(&seed, j) {
                prop_assert_eq!(crate::field::gf2_122::to_u128(c.big_b), h(TAG_DENSE) & MASK122);
            }
            if let Some(c) = gls_candidate(&seed, j) {
                let beta = gf2_61::from_u64(h(TAG_GLS) as u64);
                prop_assert_eq!(gf2_61::to_u64(c.beta), gf2_61::to_u64(beta));
            }
        }
    }
}
