//! `curve::edwards` over F_p, p = 2^107 - 1 (`field::fp107`), in
//! 14-byte encodings: x and the parity of y, decoded by a square root of
//! 105 squarings.

use crate::curve::edwards;
use crate::field::fp107::Fp;

pub type Curve = edwards::Curve<Fp>;
pub type Affine = edwards::Affine<Fp>;
pub type Point = edwards::Point<Fp>;
pub type Cached = edwards::Cached<Fp>;

#[cfg(test)]
pub(crate) mod tests {
    crate::curve::edwards::tests::suite!(fp107, Fp, fp);
}
