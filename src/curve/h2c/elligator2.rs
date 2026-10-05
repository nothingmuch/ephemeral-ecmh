//! Elligator 2 [RFC 9380, Section 6.7.1] on any curve birational to a
//! Montgomery model, written once for the models that have one.

use super::Map;
use crate::curve::encoding;
use crate::field::batch::Invert;
use crate::field::{Field, OddField, Packed};

/// A curve with a Montgomery model y^2 = x^3 + a2 x^2 + a4 x, and the
/// birational map from that model to the curve's own points. RFC 9380
/// needs a4 != 0, which every model here has, a2^2 - 4 a4 a non-square,
/// which completeness gives every model here, and a2 != 0, which
/// [`Elligator2::new`] checks.
pub trait Montgomery {
    type F: Packed;
    type Affine;
    /// (a2, a4).
    fn montgomery(&self) -> (Self::F, Self::F);
    /// The curve's point at the Montgomery model's (x, y).
    fn rational_map(&self, x: Self::F, y: Self::F) -> Self::Affine;
}

#[derive(Clone, Copy, Debug)]
pub struct Elligator2<C: Montgomery> {
    curve: C,
    a2: C::F,
    a4: C::F,
}

impl<C: Montgomery> Elligator2<C> {
    /// `None` where a2 = 0: there x1 = 0 for every input.
    pub fn new(curve: C) -> Option<Self> {
        let (a2, a4) = curve.montgomery();
        (!a2.is_zero()).then_some(Self { curve, a2, a4 })
    }

    fn rhs(&self, x: C::F) -> C::F {
        x * (x * (x + self.a2) + self.a4)
    }
}

impl<C: Montgomery> Map<C::Affine> for Elligator2<C> {
    /// With Z = `NON_SQUARE`: x1 = -a2/(1 + Z r^2), else x2 = -x1 - a2,
    /// since f(x2) = Z r^2 f(x1) for f(x) = x(x^2 + a2 x + a4). r is the low
    /// BITS bits of c reduced into the field, and sgn0(y) = sgn0(r), as
    /// RFC 9380 sets it. 1 I and 1 or 2 square roots, before the map to the
    /// curve's own model. Not uniform on its own.
    fn map(&self, c: u128) -> C::Affine {
        let r: C::F = encoding::x(c);
        let den = C::F::ONE + C::F::NON_SQUARE * r.square();
        // r^2 = -1/Z: take x1 = -a2, whose partner is x2 = 0.
        let x1 = if den.is_zero() {
            -self.a2
        } else {
            -self.a2 * den.inv()
        };
        let (x, y) = match self.rhs(x1).sqrt() {
            Some(y) => (x1, y),
            None => {
                let x2 = -x1 - self.a2;
                let y = self.rhs(x2).sqrt();
                (x2, y.expect("Elligator 2: one of x1, x2 is on the curve"))
            }
        };
        let y = if y.sgn0() != r.sgn0() { -y } else { y };
        self.curve.rational_map(x, y)
    }
}
