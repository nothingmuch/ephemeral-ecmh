//! Hashing to the curve by try-and-increment on x.

use super::{Affine, Curve};
use crate::curve::encoding::{self, Signed};
use crate::curve::h2c::{Lift, try_and_increment};
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

impl<F: OddField + Signed> Curve<F> {
    /// Try-and-increment on x. See `Sswu` for a map with per-curve setup.
    pub fn hash_to_curve(&self, h: &Salted, msg: &[u8]) -> Affine<F> {
        try_and_increment(self, h, msg)
    }
}
