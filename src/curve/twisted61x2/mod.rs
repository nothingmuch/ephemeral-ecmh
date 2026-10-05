//! `curve::twisted` over GF(p^2), p = 2^61 - 1 (`field::fp61x2`), in 16-byte
//! encodings of the quotient: u = a + b i as a | b << 61.

use crate::curve::twisted;
use crate::field::fp61x2::Fq;

pub type Curve = twisted::Curve<Fq>;
pub type Affine = twisted::Affine<Fq>;
pub type Point = twisted::Point<Fq>;
pub type Cached = twisted::Cached<Fq>;

#[cfg(test)]
pub(crate) mod tests {
    crate::curve::twisted::tests::suite!(fp61x2, Fq, fq);
}
