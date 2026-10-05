//! The policy for `twisted128`: `select`'s scheme with the witnesses read
//! in the quotient group G = E/⟨T⟩ and the bounds for p = 2^128 - 275.

use crate::curve::twisted128::{self, BYTES};
use crate::curvegen::criteria::{AdmissibleR, Criteria};
use crate::curvegen::select::{Error, Policy, eight_torsion, half};
use crate::field::fp128::{Fp, P};

pub const TAG_FP128: &[u8] = b"ephemeral-ecmh/curve/fp128";

/// floor(2 sqrt p) = 2^65 - 1: (2^65 - 1)^2 = 2^130 - 2^66 + 1 <= 4p =
/// 2^130 - 1100 < 2^130. (Squaring it would overflow u128.)
pub const HASSE: u128 = (1 << 65) - 1;
/// #E = 4r lies in p + 1 ± HASSE = 2^128 - 274 ± HASSE, which straddles
/// 2^128; `AdmissibleR` bounds r without forming 4r.
const R: AdmissibleR = AdmissibleR::new(P, HASSE, 4);
pub const R_LO: u128 = R.lo;
const _: () = assert!(P == u128::MAX - 274);

pub type Certificate = crate::curvegen::select::Certificate<BYTES>;

/// -x^2 + y^2 = 1 + d x^2 y^2 over F_p, #E = 4r, represented in G of
/// order 2r: the hashed point is cleared by 2, and the l = 8 marker is a
/// point of order 4 in G (E's points of order 4 all double to T).
pub struct Twisted128;

impl Criteria for Twisted128 {
    type Curve = twisted128::Curve;
    const TAG: &'static [u8] = TAG_FP128;
    const R: AdmissibleR = R;
    fn candidate(seed: &[u8; 32], j: u32) -> Option<Self::Curve> {
        twisted128::Curve::new(Fp::new(half(Self::TAG, seed, j)))
    }
}

impl Policy<BYTES> for Twisted128 {
    const CLEAR: u128 = 2;
    const ENTRY_PER_INDEX: bool = false;

    fn marker(c: &Self::Curve, p: &twisted128::Affine, l: u128) -> Option<bool> {
        eight_torsion(Self::CLEAR, c, p, l)
    }
}

pub fn fp128_candidate(seed: &[u8; 32], j: u32) -> Option<twisted128::Curve> {
    Twisted128::candidate(seed, j)
}

pub fn verify_fp128(seed: &[u8; 32], cert: &Certificate) -> Result<twisted128::Curve, Error> {
    crate::curvegen::select::verify::<BYTES, Twisted128>(seed, cert)
}

pub fn accept_fp128(seed: &[u8; 32], index: u32, r: u128) -> Result<twisted128::Curve, Error> {
    crate::curvegen::select::accept::<BYTES, Twisted128>(seed, index, r)
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
    fn r_bounds_are_the_hasse_interval() {
        // Express the endpoints around 2^126 to avoid overflowing p + 1 + HASSE.
        assert_eq!(R.lo, (1 << 126) - (HASSE + 274) / 4);
        assert_eq!(R.hi, (1 << 126) + (HASSE - 274) / 4);
        // r = p is out of the window: 4p exceeds p + 1 + HASSE
        const { assert!(R.hi < P) };
        // 4 R.lo >= p + 1 - HASSE > 4 (R.lo - 1), and likewise at the top,
        // with 4r - 2^128 kept in range by wrapping
        let lo = P + 1 - HASSE;
        assert!(R.lo.wrapping_mul(4) >= lo && (R.lo - 1).wrapping_mul(4) < lo);
        let hi = HASSE - 274; // p + 1 + HASSE - 2^128
        assert!(R.hi.wrapping_mul(4) <= hi && (R.hi + 1).wrapping_mul(4) > hi);
    }
}
