//! `curve::weier` over GF(p^2), p = 2^61 - 1 (`field::fp61x2`), in
//! 16-byte encodings: x = a + b i as a | b << 61, and sgn0(y).

use crate::curve::weier;
use crate::field::fp61x2::Fq;

pub type Curve = weier::Curve<Fq>;
pub type OddCurve = weier::OddCurve<Fq>;
pub type Affine = weier::Affine<Fq>;
pub type Point = weier::Point<Fq>;

#[cfg(test)]
pub(crate) mod tests {
    crate::curve::weier::tests::suite!(fp61x2, Fq, fq);
}
