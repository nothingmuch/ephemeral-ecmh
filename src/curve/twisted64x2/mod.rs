//! `curve::twisted` over GF(p^2), p = 2^64 - 59 (`field::fp64x2`), in 16-byte
//! encodings of the quotient: u = a + b i as a | b << 64.

use crate::curve::twisted;
use crate::field::fp64x2::Fq;

pub type Curve = twisted::Curve<Fq>;
pub type Affine = twisted::Affine<Fq>;
pub type Point = twisted::Point<Fq>;
pub type Cached = twisted::Cached<Fq>;

#[cfg(test)]
pub(crate) mod tests {
    crate::curve::twisted::tests::suite!(fp64x2, Fq, fq);
}
