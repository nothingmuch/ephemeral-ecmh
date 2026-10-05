//! Affine points of E\[r\], their group law and their encoding.

use super::{Constant, Curve, Model, add_batch};
use crate::field::batch::Invert;
use crate::field::{Binary, Field};

/// A point of E\[r\]; x = 0 stands for O (the 2-torsion point (0, sqrt B)
/// is not in E\[r\], so x = 0 is otherwise unused).
#[derive(Clone, Copy, Debug)]
pub struct Affine<M: Model> {
    pub x: M::F,
    pub y: M::F,
}

impl<M: Model> PartialEq for Affine<M> {
    fn eq(&self, o: &Self) -> bool {
        self.x.equals(o.x) && (self.is_identity() || self.y.equals(o.y))
    }
}
impl<M: Model> Eq for Affine<M> {}

impl<M: Model> Affine<M> {
    pub const IDENTITY: Self = Self {
        x: M::F::ZERO,
        y: M::F::ZERO,
    };

    pub fn is_identity(&self) -> bool {
        self.x.is_zero()
    }

    pub fn neg(&self) -> Self {
        Self {
            x: self.x,
            y: self.y + self.x,
        }
    }

    /// Single addition (one inversion). Use `add_batch` in bulk.
    pub fn add(&self, o: &Self) -> Self {
        let mut out = [Self::IDENTITY];
        add_batch(
            core::slice::from_ref(self),
            core::slice::from_ref(o),
            &mut out,
        );
        out[0]
    }

    /// `x | Tr(y) << SIGN`, O -> 0. Any other point has Tr(x) = 1, so x != 0.
    pub fn encode(&self) -> M::Bytes {
        if self.is_identity() {
            return M::to_bytes(0);
        }
        M::to_bytes(self.x.value() | (self.y.trace() as u128) << M::SIGN)
    }
}

impl<M: Model> Curve<M> {
    /// Is x (nonzero) the x-coordinate of a point of E\[r\]? Tr(x) = 1 and
    /// Tr(b/x) = 0: an inversion, an m_b and two parities.
    pub fn x_on_curve(&self, x: M::F) -> bool {
        x.trace() == 1 && self.b.mul(x.inv()).trace() == 0
    }

    /// Inverse of `Affine::encode`; rejects non-canonical and off-group inputs.
    pub fn decode(&self, c: M::Bytes) -> Option<Affine<M>> {
        match parse::<M>(c) {
            Parsed::Identity => Some(Affine::IDENTITY),
            Parsed::Invalid => None,
            Parsed::X(x, sign) => self.decode_with_inverse(x, x.inv(), sign),
        }
    }

    /// `decode` on many encodings, with one shared inversion.
    pub fn decode_batch(&self, cs: &[M::Bytes]) -> Vec<Option<Affine<M>>> {
        let ps: Vec<Parsed<M::F>> = cs.iter().map(|&c| parse::<M>(c)).collect();
        let mut u: Vec<M::F> = ps
            .iter()
            .map(|p| {
                if let Parsed::X(x, _) = *p {
                    x
                } else {
                    M::F::ONE
                }
            })
            .collect();
        crate::field::batch::invert(&mut u);
        ps.into_iter()
            .zip(u)
            .map(|(p, u)| match p {
                Parsed::Identity => Some(Affine::IDENTITY),
                Parsed::Invalid => None,
                Parsed::X(x, sign) => self.decode_with_inverse(x, u, sign),
            })
            .collect()
    }

    /// `decode` given x (Tr(x) = 1) and u = 1/x, so callers can batch the inversion.
    pub fn decode_with_inverse(&self, x: M::F, u: M::F, sign: u32) -> Option<Affine<M>> {
        let mut y = x * self.solve_lambda(x, u)?;
        if y.trace() != sign {
            y += x;
        }
        Some(Affine { x, y })
    }

    /// A root w of w^2 + w = x + a + B/x^2, so y = wx, if one exists.
    ///
    /// y^2 + xy = x^3 + a x^2 + B  <=>  (y/x)^2 + y/x = x + a + B/x^2, which
    /// is solvable iff its trace is 0: Tr(x) = Tr(a) = 1 cancel, and
    /// Tr(B/x^2) = Tr(b/x).
    pub(crate) fn solve_lambda(&self, x: M::F, u: M::F) -> Option<M::F> {
        if self.b.mul(u).trace() != 0 {
            return None;
        }
        Some((x + M::A + self.big_b.mul(u.square())).solve())
    }
}

/// An encoding before its inversion: O, bytes no point encodes, or a
/// nonzero x (Tr(x) = 1) and the sign.
#[derive(Clone, Copy)]
enum Parsed<F> {
    Identity,
    Invalid,
    X(F, u32),
}

fn parse<M: Model>(c: M::Bytes) -> Parsed<M::F> {
    let c = M::from_bytes(c);
    let x = M::F::new(c);
    if c == 0 {
        Parsed::Identity
    } else if c & !(M::F::MASK | 1 << M::SIGN) != 0 || x.trace() != 1 {
        Parsed::Invalid
    } else {
        Parsed::X(x, (c >> M::SIGN) as u32 & 1)
    }
}
