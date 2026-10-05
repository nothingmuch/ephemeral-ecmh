//! `binary`'s curves over GF(2^109), the 14-byte family.
//!
//! - Field: `field::gf2_109`, F_2\[z\]/(z^109 + z^5 + z^4 + z^2 + 1),
//!   with PMULL, PCLMULQDQ and portable backends. 109 is odd, so Tr(1) = 1
//!   and the doubles are still exactly the Tr(x) = 1 points
//!   (sage/gf2_109.sage).
//! - Encoding: 109 bits of x and Tr(y) at bit 109, in 14 bytes.
//! - `map`: the constants of `binary::Pornin`'s map for this modulus:
//!   crrl's with two bit parities in place of two bits.

use crate::curve::binary;
use crate::field::gf2_109::Gf;

mod map;

/// 110 bits: x, then Tr(y).
pub const ENCODED_LEN: usize = 14;

/// a = 1 over GF(2^109), in 14 bytes with the top two bits zero.
#[derive(Clone, Copy, Debug)]
pub struct M109;

impl binary::Model for M109 {
    type F = Gf;
    type K = Gf;
    const A: Gf = Gf::ONE;
    #[inline(always)]
    fn mul_a(x: Gf) -> Gf {
        x
    }
    #[inline(always)]
    fn mul_a2(x: Gf) -> Gf {
        x
    }
    #[inline(always)]
    fn mul_1a(_: Gf) -> Gf {
        Gf::ZERO
    }
    const SIGN: u32 = 109;
    type Bytes = [u8; ENCODED_LEN];
}

pub type Curve = binary::Curve<M109>;
pub type Affine = binary::Affine<M109>;
pub type Point = binary::Point<M109>;
pub use binary::{add_batch, sum_batch};

pub fn candidate(c: u128) -> (Gf, u32) {
    binary::candidate::<M109>(c)
}

#[cfg(test)]
pub(crate) mod tests {
    crate::curve::binary::tests::suite!(gf2_109, crate::curve::binary109::M109);
}
