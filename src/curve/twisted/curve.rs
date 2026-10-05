//! The curve for one namespace's d.

use super::Affine;
use crate::field::OddField;

#[derive(Clone, Copy, Debug)]
pub struct Curve<F> {
    pub d: F,
    /// 2d, for the extended and cached additions
    pub k: F,
}

impl<F: OddField> Curve<F> {
    /// `None` unless the curve is complete: -1 a square (q = 1 mod 4) and
    /// d a non-square.
    pub fn new(d: F) -> Option<Self> {
        if d.is_square() || !(-F::ONE).is_square() {
            return None;
        }
        Some(Self { d, k: d + d })
    }

    pub fn is_on_curve(&self, p: &Affine<F>) -> bool {
        let (u2, v2) = (p.u.square(), p.v.square());
        v2 - u2 == F::ONE + self.d * u2 * v2
    }
}
