//! Hashing to the curve: try-and-increment on the Montgomery x (what the
//! group traits use) or on the Edwards u. Elligator 2 is `h2c`'s, over
//! the Montgomery model these points already travel in.

use super::{Affine, Cached, Curve};
use crate::curve::encoding::{self, Signed, sign};
use crate::curve::h2c::{Lift, Montgomery, try_and_increment};
use crate::field::OddField;
use crate::hash::Salted;

/// A digest half as an encoding.
impl<F: OddField + Signed> Lift<Affine<F>> for Curve<F> {
    type Candidate = u128;
    fn candidate(&self, c: u128) -> Option<u128> {
        encoding::candidate::<F>(c)
    }
    fn lift(&self, c: u128) -> Option<Affine<F>> {
        self.decode_int(c)
    }
}

/// A digest half as an Edwards u and sgn0(v).
impl<F: OddField + Signed> Lift<Cached<F>> for Curve<F> {
    type Candidate = u128;
    fn candidate(&self, c: u128) -> Option<u128> {
        encoding::candidate::<F>(c)
    }
    fn lift(&self, c: u128) -> Option<Cached<F>> {
        self.edwards_from_u(c)
    }
}

impl<F: OddField + Signed> Curve<F> {
    /// Try-and-increment on x: about 2 square roots, no inversion.
    pub fn hash_to_curve(&self, h: &Salted, msg: &[u8]) -> Affine<F> {
        try_and_increment(self, h, msg)
    }

    /// Try-and-increment on the Edwards u instead: v^2 = (1 - u^2)/(1 - d u^2),
    /// and 1 - d u^2 != 0 since d is a non-square. One `sqrt_ratio` both
    /// tests and takes the square root of the ratio, and the result is
    /// already the `Cached` form, so there's no Montgomery-to-Edwards
    /// conversion (and its inversion) before 8M additions.
    pub fn hash_to_edwards(&self, h: &Salted, msg: &[u8]) -> Cached<F> {
        try_and_increment(self, h, msg)
    }

    /// The point with u = the low BITS bits of c reduced into the field.
    /// Bit BITS selects sgn0(v) when v != 0; higher bits are ignored.
    pub fn edwards_from_u(&self, c: u128) -> Option<Cached<F>> {
        let u: F = encoding::x(c);
        let u2 = u.square();
        let (n, den) = (F::ONE - u2, F::ONE - self.d * u2);
        let s = F::sqrt_ratio(n, den)?;
        let v = if s.sgn0() == sign::<F>(c) { s } else { -s };
        Some(Cached {
            u,
            v,
            duv: self.d * u * v,
        })
    }
}

impl<F: OddField + Signed> Montgomery for Curve<F> {
    type F = F;
    type Affine = Affine<F>;
    /// a2 = (1 + d)/2 is 0 at d = -1, complete where -1 is a non-square.
    fn montgomery(&self) -> (F, F) {
        (self.a2, self.a4)
    }
    fn rational_map(&self, x: F, y: F) -> Affine<F> {
        Affine { x, y }
    }
}
