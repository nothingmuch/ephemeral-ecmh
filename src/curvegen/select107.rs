//! The Edwards and Weierstrass policies over p = 2^107 - 1: `select`'s
//! scheme with their own tags, the bounds for this p and 14-byte points.

use crate::curve::{edwards107, weier107};
use crate::curvegen::criteria::{AdmissibleR, Criteria};
use crate::curvegen::select::{Error, Policy, eight_torsion, half, two_torsion};
use crate::field::fp107::{BYTES, Fp, MASK107, P};

pub const TAG_FP107: &[u8] = b"ephemeral-ecmh/curve/fp107";
pub const TAG_WEIER107: &[u8] = b"ephemeral-ecmh/curve/weier107";

/// floor(2 sqrt p)
const HASSE: u128 = 25476206690103090;

pub type Certificate = crate::curvegen::select::Certificate<BYTES>;

/// Edwards x^2 + y^2 = 1 + d x^2 y^2 over F_p, #E = 4r.
pub struct Edwards107;

impl Criteria for Edwards107 {
    type Curve = edwards107::Curve;
    const TAG: &'static [u8] = TAG_FP107;
    const R: AdmissibleR = AdmissibleR::new(P, HASSE, 4);
    fn candidate(seed: &[u8; 32], j: u32) -> Option<Self::Curve> {
        edwards107::Curve::new(Fp::new(half(Self::TAG, seed, j) & MASK107))
    }
}

impl Policy<BYTES> for Edwards107 {
    const CLEAR: u128 = 4;
    const ENTRY_PER_INDEX: bool = false;

    fn marker(c: &Self::Curve, p: &edwards107::Affine, l: u128) -> Option<bool> {
        eight_torsion(Self::CLEAR, c, p, l)
    }
}

/// y^2 = x^3 - 3x + b over F_p, #E = r.
pub struct Weier107;

impl Criteria for Weier107 {
    type Curve = weier107::Curve;
    const TAG: &'static [u8] = TAG_WEIER107;
    const R: AdmissibleR = AdmissibleR::new(P, HASSE, 1);
    fn candidate(seed: &[u8; 32], j: u32) -> Option<Self::Curve> {
        weier107::Curve::new(Fp::new(half(Self::TAG, seed, j) & MASK107))
    }
}

impl Policy<BYTES> for Weier107 {
    const CLEAR: u128 = 1;
    const ENTRY_PER_INDEX: bool = false;

    fn marker(_: &Self::Curve, p: &weier107::Affine, l: u128) -> Option<bool> {
        two_torsion(p, l)
    }
}

pub fn fp107_candidate(seed: &[u8; 32], j: u32) -> Option<edwards107::Curve> {
    Edwards107::candidate(seed, j)
}

pub fn weier107_candidate(seed: &[u8; 32], j: u32) -> Option<weier107::Curve> {
    Weier107::candidate(seed, j)
}

pub fn verify_fp107(seed: &[u8; 32], cert: &Certificate) -> Result<edwards107::Curve, Error> {
    crate::curvegen::select::verify::<BYTES, Edwards107>(seed, cert)
}

pub fn accept_fp107(seed: &[u8; 32], index: u32, r: u128) -> Result<edwards107::Curve, Error> {
    crate::curvegen::select::accept::<BYTES, Edwards107>(seed, index, r)
}

pub fn verify_weier107(seed: &[u8; 32], cert: &Certificate) -> Result<weier107::Curve, Error> {
    crate::curvegen::select::verify::<BYTES, Weier107>(seed, cert)
}

pub fn accept_weier107(seed: &[u8; 32], index: u32, r: u128) -> Result<weier107::Curve, Error> {
    crate::curvegen::select::accept::<BYTES, Weier107>(seed, index, r)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curvegen::select::is_prime;

    #[test]
    fn p_is_prime() {
        assert!(is_prime(P));
        assert!(!is_prime(P + 2));
    }

    #[test]
    fn bounds_match_the_direct_constants() {
        // The floor is odd and below every admissible r, so it is a valid
        // odd rejection label.
        const R_MIN: u128 = 40564819207303334478842830046259;
        assert_eq!(Edwards107::R.lo, R_MIN + 1);
        assert!(R_MIN & 1 == 1 && (3..Edwards107::R.lo).contains(&R_MIN));
        assert_eq!(Weier107::R.lo, P + 1 - HASSE);
        // r = q is reachable only at cofactor 1
        assert!(Edwards107::R.hi < P && Weier107::R.admits(P));
        assert!(Edwards107::R.admits(P / 4) && Weier107::R.admits(P));
    }
}
