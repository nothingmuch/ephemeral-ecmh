//! `curve::edwards` over GF(p^2), p = 2^61 - 1 (`field::fp61x2`), in
//! 16-byte encodings: x = a + b i as a | b << 61, and sgn0(y).

use crate::curve::edwards;
use crate::field::fp61x2::Fq;

pub type Curve = edwards::Curve<Fq>;
pub type Affine = edwards::Affine<Fq>;
pub type Point = edwards::Point<Fq>;
pub type Cached = edwards::Cached<Fq>;

#[cfg(test)]
pub(crate) mod tests {
    crate::curve::edwards::tests::suite!(fp61x2, Fq, fq);
}
