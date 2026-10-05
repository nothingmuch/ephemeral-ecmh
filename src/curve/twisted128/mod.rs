//! `curve::twisted` over F_p, p = 2^128 - 275 (`field::fp128`; p = 5 mod
//! 8, so -1 is a square), in 16-byte encodings of the quotient: u alone.

use crate::curve::twisted;
use crate::field::fp128::Fp;

pub use twisted::BYTES;
pub type Curve = twisted::Curve<Fp>;
pub type Affine = twisted::Affine<Fp>;
pub type Point = twisted::Point<Fp>;
pub type Cached = twisted::Cached<Fp>;

#[cfg(test)]
pub(crate) mod tests {
    crate::curve::twisted::tests::suite!(fp128, Fp, fp);

    #[test]
    fn encodings_from_p_up_are_rejected() {
        use crate::field::fp128::P;
        let c = Curve::new(F::ONE + F::ONE).unwrap();
        for e in P..=u128::MAX {
            assert!(c.decode_u128(e).is_none(), "{e}");
        }
    }
}
