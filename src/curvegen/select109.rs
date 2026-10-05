//! The binary policy over GF(2^109): `select`'s scheme with its own tag,
//! the bounds for q = 2^109 and 14-byte point encodings.
//!
//! GHS needs no check: ord_109(2) = 36, so every B outside F_2 descends
//! to genus at least 2^35 (sage/gf2_109.sage).

use crate::curve::binary109::{self, ENCODED_LEN};
use crate::curvegen::criteria::{AdmissibleR, Criteria};
use crate::curvegen::select::{Error, Policy, even_order, half};
use crate::field::gf2_109::{MASK109, from_u128};

pub const TAG_GF2_109: &[u8] = b"ephemeral-ecmh/curve/gf2-109";

const Q: u128 = 1 << 109;
/// floor(2 sqrt q)
const HASSE: u128 = 50952413380206180;
const _: () = assert!(HASSE * HASSE <= 4 * Q && (HASSE + 1) * (HASSE + 1) > 4 * Q);

/// Indexed rejection witnesses with 14-byte point encodings.
pub type Certificate = crate::curvegen::select::Certificate<ENCODED_LEN>;

/// y^2 + xy = x^3 + x^2 + B over GF(2^109), #E = 2r.
///
/// The marker is the binary one, as for [`super::select::Binary127`]: an
/// even label, #E of a candidate whose odd part is composite.
pub struct Binary109;

impl Criteria for Binary109 {
    type Curve = binary109::Curve;
    const TAG: &'static [u8] = TAG_GF2_109;
    const R: AdmissibleR = AdmissibleR::new(Q, HASSE, 2);
    fn candidate(seed: &[u8; 32], j: u32) -> Option<Self::Curve> {
        let v = half(Self::TAG, seed, j) & MASK109;
        (v != 0).then(|| binary109::Curve::new(from_u128(v)))
    }
}

impl Policy<ENCODED_LEN> for Binary109 {
    const CLEAR: u128 = 1;
    const ENTRY_PER_INDEX: bool = true;

    fn marker(c: &Self::Curve, p: &binary109::Affine, l: u128) -> Option<bool> {
        even_order::<ENCODED_LEN, _>(&Self::R, c, p, l)
    }
}

pub fn candidate(seed: &[u8; 32], j: u32) -> Option<binary109::Curve> {
    Binary109::candidate(seed, j)
}

pub fn verify(seed: &[u8; 32], cert: &Certificate) -> Result<binary109::Curve, Error> {
    crate::curvegen::select::verify::<ENCODED_LEN, Binary109>(seed, cert)
}

pub fn accept(seed: &[u8; 32], index: u32, r: u128) -> Result<binary109::Curve, Error> {
    crate::curvegen::select::accept::<ENCODED_LEN, Binary109>(seed, index, r)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bounds_match_the_direct_constants() {
        // The floor is even, so the odd rejection labels end below it.
        const R_MIN: u128 = 324518553658426701306949330473166;
        assert_eq!(Binary109::R.lo, R_MIN + 1);
        const { assert!(R_MIN & 1 == 0) };
        assert!(Binary109::R.admits(Q / 2) && !Binary109::R.admits(1 << 126));
        // r = q is out of the window: 2q exceeds q + 1 + HASSE
        const { assert!(Binary109::R.hi < Q) };
    }
}
