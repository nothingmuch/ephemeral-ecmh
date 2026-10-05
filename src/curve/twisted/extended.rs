//! Extended Edwards (X : Y : Z : T) coordinates, the accumulators, and
//! the cached affine addends.

use super::{Affine, Curve};
use crate::field::{OddField, batch};

/// Extended Edwards (X : Y : Z : T), u = X/Z, v = Y/Z, T = XY/Z.
#[derive(Clone, Copy, Debug)]
pub struct Point<F> {
    pub x: F,
    pub y: F,
    pub z: F,
    pub t: F,
}

/// (v - u, v + u, 2d u v), for 7M mixed additions.
#[derive(Clone, Copy, Debug)]
pub struct Cached<F> {
    pub ym: F,
    pub yp: F,
    pub t2d: F,
}

impl<F: OddField> Cached<F> {
    /// -(u, v) = (-u, v): swap v - u and v + u.
    pub fn neg(&self) -> Self {
        Self {
            ym: self.yp,
            yp: self.ym,
            t2d: -self.t2d,
        }
    }
}

impl<F: OddField> Point<F> {
    pub const IDENTITY: Self = Self {
        x: F::ZERO,
        y: F::ONE,
        z: F::ONE,
        t: F::ZERO,
    };

    pub fn neg(&self) -> Self {
        Self {
            x: -self.x,
            t: -self.t,
            ..*self
        }
    }

    /// In G: u = 0 is O or T.
    pub fn is_identity(&self) -> bool {
        self.x.is_zero()
    }

    /// Equal in G: the same point, or one of them plus T, which negates u
    /// and v. Both lie on one line through the origin, (X1 : Y1) =
    /// (X2 : Y2), and that line meets E only in that class (see
    /// `equals_cached`).
    pub fn equals(&self, o: &Self) -> bool {
        self.x * o.y == o.x * self.y
    }

    /// Equal in G to a cached addend: (X : Y) = (u : v), as
    /// (Y + X : Y - X) = (v + u : v - u). Of the line's points λ(u, v) only
    /// ±(u, v), one class, lie on E: λ^2 = 1, or else d (λ u v)^2 = 1 with
    /// d a non-square.
    pub fn equals_cached(&self, a: &Cached<F>) -> bool {
        (self.y + self.x) * a.ym == (self.y - self.x) * a.yp
    }
}

impl<F: OddField> Curve<F> {
    pub fn from_affine(&self, p: &Affine<F>) -> Point<F> {
        Point {
            x: p.u,
            y: p.v,
            z: F::ONE,
            t: p.u * p.v,
        }
    }

    /// 2M, no inversion: the affine point is the addend.
    pub fn to_cached(&self, p: &Affine<F>) -> Cached<F> {
        Cached {
            ym: p.v - p.u,
            yp: p.v + p.u,
            t2d: self.k * p.u * p.v,
        }
    }

    /// One inversion of Z.
    pub fn to_affine(&self, p: &Point<F>) -> Affine<F> {
        Self::normalize_with(p, batch::Invert::inv(p.z))
    }

    fn normalize_with(p: &Point<F>, zi: F) -> Affine<F> {
        Affine {
            u: p.x * zi,
            v: p.y * zi,
        }
    }

    /// Z is never 0 on a complete curve, so every point takes part.
    pub fn normalize_batch(&self, ps: &[Point<F>]) -> Vec<Affine<F>> {
        let mut zi: Vec<F> = ps.iter().map(|p| p.z).collect();
        batch::invert(&mut zi);
        ps.iter()
            .zip(zi)
            .map(|(p, zi)| Self::normalize_with(p, zi))
            .collect()
    }

    /// Complete a = -1 addition (add-2008-hwcd-3): 8M + 1 m_2d.
    pub fn add_ext(&self, p: &Point<F>, q: &Point<F>) -> Point<F> {
        let a = (p.y - p.x) * (q.y - q.x);
        let b = (p.y + p.x) * (q.y + q.x);
        let c = p.t * self.k * q.t;
        let zz = p.z * q.z;
        Self::finish(a, b, c, zz + zz)
    }

    /// Extended + cached: 7M.
    pub fn add_cached(&self, p: &Point<F>, q: &Cached<F>) -> Point<F> {
        let a = (p.y - p.x) * q.ym;
        let b = (p.y + p.x) * q.yp;
        let c = p.t * q.t2d;
        Self::finish(a, b, c, p.z + p.z)
    }

    fn finish(a: F, b: F, c: F, d: F) -> Point<F> {
        let (e, f, g, h) = (b - a, d - c, d + c, b + a);
        Point {
            x: e * f,
            y: g * h,
            z: f * g,
            t: e * h,
        }
    }

    /// Fold the cached additions and normalise once: 9M per point.
    pub fn sum_batch(&self, ps: &[Affine<F>]) -> Affine<F> {
        let acc = ps.iter().fold(Point::IDENTITY, |acc, p| {
            self.add_cached(&acc, &self.to_cached(p))
        });
        self.to_affine(&acc)
    }

    /// Scalars are public (orders and witnesses), so the loop skips the
    /// leading zeros.
    pub fn mul(&self, p: &Point<F>, k: u128) -> Point<F> {
        let mut acc = Point::IDENTITY;
        for i in (0..128 - k.leading_zeros()).rev() {
            acc = self.add_ext(&acc, &acc);
            if k >> i & 1 == 1 {
                acc = self.add_ext(&acc, p);
            }
        }
        acc
    }
}
