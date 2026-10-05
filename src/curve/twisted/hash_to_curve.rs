//! Hashing to the curve: try-and-increment on u, and the Montgomery model
//! that `h2c`'s Elligator 2 maps through.

use super::{Affine, Curve};
use crate::curve::h2c::{Lift, Montgomery, try_and_increment};
use crate::field::Packed;
use crate::hash::Salted;

/// A digest half's low BITS bits as an encoding; `decode_u128` rejects
/// those that don't unpack.
impl<F: Packed> Lift<Affine<F>> for Curve<F> {
    type Candidate = u128;
    fn candidate(&self, c: u128) -> Option<u128> {
        Some(c & u128::MAX >> (128 - F::BITS))
    }
    fn lift(&self, c: u128) -> Option<Affine<F>> {
        self.decode_u128(c)
    }
}

impl<F: Packed> Curve<F> {
    /// Each digest half is a candidate encoding; about half decode, one
    /// `sqrt_ratio` each. For uniform candidates this is exactly uniform
    /// on G: skipping those that don't unpack leaves them uniform on F,
    /// and each class has exactly one encoding.
    pub fn hash_to_curve(&self, h: &Salted, msg: &[u8]) -> Affine<F> {
        try_and_increment(self, h, msg)
    }
}

/// A curve with its Montgomery model y^2 = x^3 + a2 x^2 + a4 x, where
/// a2 = (d - 1)/2 and a4 = (1 + d)^2/16: B y^2 = x^3 + A x^2 + x scaled by
/// B = -4/(1 + d). a2^2 - 4 a4 = -d is a non-square as d is, and a2 = 0
/// only at d = 1, a square: Elligator 2 applies on every curve. B costs an
/// inversion that only the Montgomery maps use, so it is kept here rather
/// than in [`Curve`], which every candidate in curve generation builds.
#[derive(Clone, Copy, Debug)]
pub struct MontgomeryModel<F> {
    curve: Curve<F>,
    b: F,
}

impl<F: Packed> MontgomeryModel<F> {
    pub fn new(curve: Curve<F>) -> Self {
        // 1 + d != 0: -1 is a square and d is not
        let four = (F::ONE + F::ONE).square();
        let b = -four * (F::ONE + curve.d).inv();
        Self { curve, b }
    }
}

impl<F: Packed> Montgomery for MontgomeryModel<F> {
    type F = F;
    type Affine = Affine<F>;
    fn montgomery(&self) -> (F, F) {
        let two = F::ONE + F::ONE;
        let (d, e) = (self.curve.d, F::ONE + self.curve.d);
        ((d - F::ONE) * two.inv(), (e * (two + two).inv()).square())
    }
    /// u = x_M/y_M and v = (x_M - 1)/(x_M + 1) with x_M = B x, y_M = B y,
    /// one inversion for both. (0, 0) is T. x_M = -1 has no rational point:
    /// its points map to E's points at infinity, which -d a non-square
    /// keeps irrational.
    fn rational_map(&self, x: F, y: F) -> Affine<F> {
        if y.is_zero() {
            return Affine {
                u: F::ZERO,
                v: -F::ONE,
            };
        }
        let xm = self.b * x;
        let (p, m) = (xm + F::ONE, xm - F::ONE);
        let inv = (y * p).inv();
        Affine {
            u: x * p * inv,
            v: m * y * inv,
        }
    }
}
