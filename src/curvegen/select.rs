use super::criteria::{AdmissibleR, Criteria};
pub use super::number::{EMBEDDING_MIN, embedding_degree_ok, is_prime};
use crate::curve::binary127;
use crate::field::gf2_127::{MASK127, from_u128};
use crate::hash::{Salted, halves};
pub const TAG_GF2_127: &[u8] = b"ephemeral-ecmh/curve/gf2";
/// floor(2 sqrt q), the same for q = 2^127 and p = 2^127 - 1.
const HASSE: u128 = 26087635650665564424;
/// The first 128-bit half of the tagged digest of j: the raw material of
/// every candidate.
pub(crate) fn half(tag: &[u8], seed: &[u8; 32], j: u32) -> u128 {
    halves(&Salted::new(tag, seed).digest(&[], j))[0]
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
pub fn gf2_127_candidate(seed: &[u8; 32], j: u32) -> Option<binary127::Curve> {
    Binary127::candidate(seed, j)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::curvegen::number::mul_wide;
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
        {
            let r = Binary127::R;
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
    }
}
