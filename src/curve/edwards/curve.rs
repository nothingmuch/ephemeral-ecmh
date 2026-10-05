//! The curve for one namespace's d, and its Montgomery model's constants.

use super::Affine;
use crate::field::OddField;

#[derive(Clone, Copy, Debug)]
pub struct Curve<F> {
    pub d: F,
    pub a2: F,
    pub a4: F,
    /// B = 4/(1-d): Montgomery x_M = B x maps to v = (x_M - 1)/(x_M + 1)
    pub bm: F,
}

impl<F: OddField> Curve<F> {
    /// `None` unless d is a non-square, i.e. the Edwards model is complete.
    pub fn new(d: F) -> Option<Self> {
        if d.sqrt().is_some() {
            return None;
        }
        let e = F::ONE - d; // nonzero: d = 1 is a square
        let two = F::ONE + F::ONE;
        let four = two + two;
        Some(Self {
            d,
            a2: (F::ONE + d) * two.inv(),
            a4: (e * four.inv()).square(),
            bm: four * e.inv(),
        })
    }

    pub fn rhs(&self, x: F) -> F {
        x * (x * (x + self.a2) + self.a4)
    }

    /// Whether x is a Montgomery x-coordinate.
    pub fn x_on_curve(&self, x: F) -> bool {
        self.rhs(x).is_square()
    }

    /// Is u an Edwards u-coordinate? (1 - u^2)/(1 - d u^2) is a square iff
    /// its numerator times denominator is.
    pub fn u_on_curve(&self, u: F) -> bool {
        let u2 = u.square();
        ((F::ONE - u2) * (F::ONE - self.d * u2)).is_square()
    }

    pub fn is_on_curve(&self, p: &Affine<F>) -> bool {
        p.is_identity() || p.y.square() == self.rhs(p.x)
    }
}
