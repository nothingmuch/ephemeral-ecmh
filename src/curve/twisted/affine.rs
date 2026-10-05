//! Edwards affine points and the 16-byte encoding of their class in G.

use super::{BYTES, Curve};
use crate::field::{OddField, Packed};

/// Edwards affine (u, v), any representative of its class in G.
#[derive(Clone, Copy, Debug)]
pub struct Affine<F> {
    pub u: F,
    pub v: F,
}

impl<F: OddField> Affine<F> {
    pub const IDENTITY: Self = Self {
        u: F::ZERO,
        v: F::ONE,
    };

    /// Equal in G: (u, v) and (-u, -v) are one class.
    pub fn equals(&self, o: &Self) -> bool {
        (self.u == o.u && self.v == o.v) || (self.u == -o.u && self.v == -o.v)
    }

    pub fn is_identity(&self) -> bool {
        self.u.is_zero()
    }

    pub fn neg(&self) -> Self {
        Self {
            u: -self.u,
            v: self.v,
        }
    }
}

impl<F: Packed> Affine<F> {
    /// The u of the representative with sgn0(v) = 0, or with sgn0(u) = 0
    /// when v = 0 (the order-4 points (±i, 0), which T swaps).
    pub fn encode(&self) -> [u8; BYTES] {
        let flip = if self.v.is_zero() {
            self.u.sgn0()
        } else {
            self.v.sgn0()
        };
        let u = if flip { -self.u } else { self.u };
        u.pack().to_le_bytes()
    }
}

impl<F: Packed> Curve<F> {
    /// One `sqrt_ratio`; its cost is the field's.
    pub fn decode(&self, c: [u8; BYTES]) -> Option<Affine<F>> {
        self.decode_u128(u128::from_le_bytes(c))
    }

    pub(crate) fn decode_u128(&self, c: u128) -> Option<Affine<F>> {
        let u = F::unpack(c)?;
        let u2 = u.square();
        // v^2 = (1 + u^2)/(1 - d u^2); d is a non-square, so 1 - d u^2 != 0
        let s = F::sqrt_ratio(F::ONE + u2, F::ONE - self.d * u2)?;
        let v = if s.sgn0() { -s } else { s };
        if v.is_zero() && u.sgn0() {
            return None;
        }
        Some(Affine { u, v })
    }
}
