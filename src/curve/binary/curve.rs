//! What a family fixes (`Model`), and the curve for one namespace's B with
//! the constants derived from it.

use super::Affine;
use crate::field::{Binary, Field};
use core::fmt::Debug;

/// A family of curves y^2 + xy = x^3 + a x^2 + B over one field: a, the
/// type B is drawn from, and the encoding's layout.
///
/// Tr(a) = 1 makes #E = 2 * odd possible and puts the doubles, E\[r\] = 2E,
/// exactly on the hyperplane Tr(x) = 1. For odd m that is a = 1, which the
/// formulas then drop; for even m, Tr(1) = 0 and a must be some other
/// trace-1 element, chosen so that products by a^2 and 1 + a take as few
/// operations as the field allows.
pub trait Model: Copy + Debug + 'static {
    type F: Binary;
    /// B, b = sqrt(B) and beta = sqrt(b).
    type K: Constant<Self::F>;
    const A: Self::F;
    /// a x, a^2 x and (1 + a) x.
    fn mul_a(x: Self::F) -> Self::F;
    fn mul_a2(x: Self::F) -> Self::F;
    fn mul_1a(x: Self::F) -> Self::F;
    /// The bit of the encoding's integer that carries Tr(y), outside
    /// `F::MASK`; every other bit outside the mask must be 0.
    const SIGN: u32;
    /// At most 16 bytes, little-endian.
    type Bytes: Copy + Eq + Debug + Default + AsRef<[u8]> + AsMut<[u8]>;

    /// v's low bytes, as many as fit.
    fn to_bytes(v: u128) -> Self::Bytes {
        let mut b = Self::Bytes::default();
        let n = b.as_ref().len();
        b.as_mut().copy_from_slice(&v.to_le_bytes()[..n]);
        b
    }

    fn from_bytes(b: Self::Bytes) -> u128 {
        let mut v = [0; 16];
        v[..b.as_ref().len()].copy_from_slice(b.as_ref());
        u128::from_le_bytes(v)
    }
}

/// A curve constant: B, b or beta, as the type the family draws it from.
pub trait Constant<F>: Copy + Debug {
    /// The constant as a field element.
    fn gf(self) -> F;
    fn mul(self, x: F) -> F;
    fn sqrt(self) -> Self;
    /// 1/self, for nonzero self.
    fn recip(self) -> Self;
}

/// A dense constant: each m_beta is a full multiplication.
impl<F: Binary> Constant<F> for F {
    #[inline(always)]
    fn gf(self) -> F {
        self
    }
    #[inline(always)]
    fn mul(self, x: F) -> F {
        self * x
    }
    fn sqrt(self) -> Self {
        Binary::sqrt(self)
    }
    fn recip(self) -> Self {
        crate::field::batch::Invert::inv(self)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Curve<M: Model> {
    pub big_b: M::K,
    /// b = sqrt(B)
    pub b: M::K,
    /// beta = sqrt(b)
    pub beta: M::K,
}

impl<M: Model> Curve<M> {
    /// B must be nonzero.
    pub fn new(big_b: M::K) -> Self {
        let b = big_b.sqrt();
        Self {
            big_b,
            b,
            beta: b.sqrt(),
        }
    }

    pub fn is_on_curve(&self, p: &Affine<M>) -> bool {
        if p.is_identity() {
            return true;
        }
        let (x, y) = (p.x, p.y);
        let x2 = x.square();
        (y.square() + x * y).equals(x2 * x + M::mul_a(x2) + self.big_b.gf())
    }
}
