//! `binary`'s curves over GF(2^127), the 16-byte family.
//!
//! - Field: `field::gf2_127`, F_2\[z\]/(z^127 + z^63 + 1), crrl's.
//! - Encoding: 127 bits of x and Tr(y) at bit 127, every bit of 16 bytes.
//! - `map`: the constants of `binary::Pornin`'s map for this modulus,
//!   which are crrl's.

use crate::curve::binary;
use crate::field::gf2_127::Gf;

mod map;

/// a = 1 over GF(2^127), with every bit of 16 bytes in use.
#[derive(Clone, Copy, Debug)]
pub struct M127;

impl binary::Model for M127 {
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
    const SIGN: u32 = 127;
    type Bytes = [u8; 16];
}

pub type Curve = binary::Curve<M127>;
pub type Affine = binary::Affine<M127>;
pub use binary::{add_batch, sum_batch};
pub type MapState = binary::MapState<Gf>;

pub fn candidate(c: u128) -> (Gf, u32) {
    binary::candidate::<M127>(c)
}

#[cfg(test)]
pub(crate) mod tests {
    crate::curve::binary::tests::suite!(gf2_127, crate::curve::binary127::M127);

    proptest! {
        #[test]
        fn candidate_sets_bit_0(c in any::<u128>()) {
            // Tr(x) = x_0 here, so fixing the trace is setting bit 0
            let (x, sign) = candidate::<M>(c);
            prop_assert!(Field::equals(x, Gf::new(c | 1)));
            prop_assert_eq!(sign as u128, c >> 127);
        }
    }
}
