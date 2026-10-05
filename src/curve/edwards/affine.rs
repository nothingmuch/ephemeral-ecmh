//! Montgomery affine points, and their encoding.

use super::Curve;
use crate::curve::encoding::{self, Signed};
use crate::field::OddField;

/// Montgomery affine point; O is the off-curve sentinel (0, 1).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Affine<F> {
    pub x: F,
    pub y: F,
}

impl<F: OddField> Affine<F> {
    pub const IDENTITY: Self = Self {
        x: F::ZERO,
        y: F::ONE,
    };

    pub fn is_identity(&self) -> bool {
        *self == Self::IDENTITY
    }

    pub fn neg(&self) -> Self {
        if self.is_identity() {
            *self
        } else {
            Self {
                x: self.x,
                y: -self.y,
            }
        }
    }
}

impl<F: OddField + Signed> Affine<F> {
    /// x and sgn0(y) (`curve::encoding`).
    pub fn encode(&self) -> F::Bytes {
        encoding::encode((!self.is_identity()).then_some((self.x, self.y)))
    }
}

impl<F: OddField + Signed> Curve<F> {
    /// One square root, no inversion.
    pub fn decode(&self, c: F::Bytes) -> Option<Affine<F>> {
        self.decode_int(F::from_bytes(c))
    }

    pub(super) fn decode_int(&self, c: u128) -> Option<Affine<F>> {
        let p = encoding::decode(c, |x| self.rhs(x))?;
        Some(p.map_or(Affine::IDENTITY, |(x, y)| Affine { x, y }))
    }
}
