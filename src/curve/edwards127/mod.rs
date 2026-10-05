//! `curve::edwards` over F_p, p = 2^127 - 1 (`field::fp127`), in
//! 16-byte encodings: x and the parity of y, decoded by a square root of
//! 125 squarings.

use crate::curve::edwards;
use crate::field::fp127::Fp;

pub type Curve = edwards::Curve<Fp>;
pub type Affine = edwards::Affine<Fp>;
pub type Point = edwards::Point<Fp>;
pub type Cached = edwards::Cached<Fp>;

#[cfg(test)]
pub(crate) mod tests {
    crate::curve::edwards::tests::suite!(fp127, Fp, fp);
}
