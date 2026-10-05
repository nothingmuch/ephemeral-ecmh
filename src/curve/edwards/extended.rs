//! Extended Edwards (X : Y : Z : T) coordinates, the accumulators, and
//! the cached Edwards affine addends.

use super::{Affine, Curve};
use crate::field::OddField;
use crate::field::batch::invert as batch_invert;

/// Extended Edwards (X : Y : Z : T), u = X/Z, v = Y/Z, T = XY/Z.
#[derive(Clone, Copy, Debug)]
pub struct Point<F> {
    pub x: F,
    pub y: F,
    pub z: F,
    pub t: F,
}

/// Edwards affine with d*u*v precomputed, for 8M mixed additions.
#[derive(Clone, Copy, Debug)]
pub struct Cached<F> {
    pub u: F,
    pub v: F,
    pub duv: F,
}

impl<F: OddField> Cached<F> {
    pub fn neg(&self) -> Self {
        Self {
            u: -self.u,
            duv: -self.duv,
            ..*self
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

    pub fn equals(&self, o: &Self) -> bool {
        self.x * o.z == o.x * self.z && self.y * o.z == o.y * self.z
    }

    /// (0 : Z : Z : 0); (0 : -Z : Z : 0) is the point of order 2.
    pub fn is_identity(&self) -> bool {
        self.x.is_zero() && self.y == self.z
    }
}

impl<F: OddField> Curve<F> {
    /// u = x/y, v = (w-1)/(w+1) with w = B x, without inversion (5M).
    /// Only (0,0) and O are exceptional: other y = 0 points, and w = -1,
    /// would need d to be a square.
    pub fn from_affine(&self, p: &Affine<F>) -> Point<F> {
        if p.is_identity() {
            return Point::IDENTITY;
        }
        if p.y.is_zero() {
            return Point {
                x: F::ZERO,
                y: -F::ONE,
                z: F::ONE,
                t: F::ZERO,
            };
        }
        let w = self.bm * p.x;
        let (wm, wp) = (w - F::ONE, w + F::ONE);
        Point {
            x: p.x * wp,
            y: wm * p.y,
            z: p.y * wp,
            t: p.x * wm,
        }
    }

    /// x = (Z+Y)/(B (Z-Y)), y = x Z/X: one inversion of B (Z-Y) X.
    pub fn to_affine(&self, p: &Point<F>) -> Affine<F> {
        let den = self.bm * (p.z - p.y) * p.x;
        self.normalize_with(p, den.inv())
    }

    fn normalize_with(&self, p: &Point<F>, inv: F) -> Affine<F> {
        if p.x.is_zero() {
            return if p.y == p.z {
                Affine::IDENTITY
            } else {
                Affine {
                    x: F::ZERO,
                    y: F::ZERO,
                }
            };
        }
        let zy = p.z + p.y;
        Affine {
            x: zy * p.x * inv,
            y: zy * p.z * inv,
        }
    }

    pub fn normalize_batch(&self, ps: &[Point<F>]) -> Vec<Affine<F>> {
        let mut den: Vec<F> = ps
            .iter()
            .map(|p| {
                if p.x.is_zero() {
                    F::ONE
                } else {
                    self.bm * (p.z - p.y) * p.x
                }
            })
            .collect();
        batch_invert(&mut den);
        ps.iter()
            .zip(den)
            .map(|(p, i)| self.normalize_with(p, i))
            .collect()
    }

    pub fn to_cached_batch(&self, ps: &[Affine<F>]) -> Vec<Cached<F>> {
        let ext: Vec<Point<F>> = ps.iter().map(|p| self.from_affine(p)).collect();
        let mut zi: Vec<F> = ext.iter().map(|p| p.z).collect();
        batch_invert(&mut zi);
        ext.iter()
            .zip(zi)
            .map(|(p, zi)| {
                let (u, v) = (p.x * zi, p.y * zi);
                Cached {
                    u,
                    v,
                    duv: self.d * u * v,
                }
            })
            .collect()
    }

    /// Complete a = 1 addition (Hisil–Wong–Carter–Dawson): 9M + 1 m_d.
    pub fn add_ext(&self, p: &Point<F>, q: &Point<F>) -> Point<F> {
        let c = self.d * p.t * q.t;
        Self::finish(p, p.x * q.x, p.y * q.y, c, p.z * q.z, q.x + q.y)
    }

    /// Extended + cached Edwards affine: 8M.
    pub fn add_cached(&self, p: &Point<F>, q: &Cached<F>) -> Point<F> {
        Self::finish(p, p.x * q.u, p.y * q.v, p.t * q.duv, p.z, q.u + q.v)
    }

    fn finish(p: &Point<F>, a: F, b: F, c: F, d: F, xy2: F) -> Point<F> {
        let e = (p.x + p.y) * xy2 - (a + b);
        let (f, g, h) = (d - c, d + c, b - a);
        Point {
            x: e * f,
            y: g * h,
            z: f * g,
            t: e * h,
        }
    }

    /// Scalars are public (orders and witnesses), so the loop skips their
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
