//! `curve::weier` over F_p, p = 2^127 - 1 (`field::fp127`), in
//! 16-byte encodings: x and the parity of y, decoded by a square root of
//! 125 squarings.

use crate::curve::weier;
use crate::field::fp127::Fp;

pub type Curve = weier::Curve<Fp>;
pub type OddCurve = weier::OddCurve<Fp>;
pub type Affine = weier::Affine<Fp>;
pub type Point = weier::Point<Fp>;

#[cfg(test)]
pub(crate) mod tests {
    crate::curve::weier::tests::suite!(fp127, Fp, fp);
}
