//! Projective (X : Y : Z) coordinates, the accumulators, with the
//! Renes–Costello–Batina complete addition.

use super::{Affine, Curve};
use crate::field::OddField;
use crate::field::batch::invert as batch_invert;

/// Projective (X : Y : Z), x = X/Z, y = Y/Z; O = (0 : 1 : 0).
#[derive(Clone, Copy, Debug)]
pub struct Point<F> {
    pub x: F,
    pub y: F,
    pub z: F,
}

impl<F: OddField> Point<F> {
    pub const IDENTITY: Self = Self {
        x: F::ZERO,
        y: F::ONE,
        z: F::ZERO,
    };

    pub fn neg(&self) -> Self {
        Self {
            y: -self.y,
            ..*self
        }
    }

    /// Strict: (0 : 0 : 0), which RCB returns for exceptional inputs on
    /// curves with 2-torsion (never on `OddCurve`), is not O.
    pub fn is_identity(&self) -> bool {
        self.x.is_zero() && self.z.is_zero() && !self.y.is_zero()
    }

    /// Only meaningful on `OddCurve`, where (0 : 0 : 0) cannot occur.
    pub fn equals(&self, o: &Self) -> bool {
        self.x * o.z == o.x * self.z && self.y * o.z == o.y * self.z
    }
}

impl<F: OddField> Curve<F> {
    pub fn from_affine(&self, p: &Affine<F>) -> Point<F> {
        if p.is_identity() {
            Point::IDENTITY
        } else {
            Point {
                x: p.x,
                y: p.y,
                z: F::ONE,
            }
        }
    }

    pub fn to_affine(&self, p: &Point<F>) -> Affine<F> {
        self.normalize_with(p, p.z.inv())
    }

    fn normalize_with(&self, p: &Point<F>, zi: F) -> Affine<F> {
        if p.z.is_zero() {
            Affine::IDENTITY
        } else {
            Affine {
                x: p.x * zi,
                y: p.y * zi,
            }
        }
    }

    pub fn normalize_batch(&self, ps: &[Point<F>]) -> Vec<Affine<F>> {
        let mut zi: Vec<F> = ps
            .iter()
            .map(|p| if p.z.is_zero() { F::ONE } else { p.z })
            .collect();
        batch_invert(&mut zi);
        ps.iter()
            .zip(zi)
            .map(|(p, zi)| self.normalize_with(p, zi))
            .collect()
    }

    /// RCB Alg. 4, transcribed step by step.
    pub fn add_ext(&self, p: &Point<F>, q: &Point<F>) -> Point<F> {
        let b = self.b;
        let (x1, y1, z1, x2, y2, z2) = (p.x, p.y, p.z, q.x, q.y, q.z);
        let mut t0 = x1 * x2;
        let mut t1 = y1 * y2;
        let mut t2 = z1 * z2;
        let mut t3 = x1 + y1;
        let mut t4 = x2 + y2;
        t3 = t3 * t4;
        t4 = t0 + t1;
        t3 = t3 - t4;
        t4 = y1 + z1;
        let mut x3 = y2 + z2;
        t4 = t4 * x3;
        x3 = t1 + t2;
        t4 = t4 - x3;
        x3 = x1 + z1;
        let mut y3 = x2 + z2;
        x3 = x3 * y3;
        y3 = t0 + t2;
        y3 = x3 - y3;
        let mut z3 = b * t2;
        x3 = y3 - z3;
        z3 = x3 + x3;
        x3 = x3 + z3;
        z3 = t1 - x3;
        x3 = t1 + x3;
        y3 = b * y3;
        t1 = t2 + t2;
        t2 = t1 + t2;
        y3 = y3 - t2;
        y3 = y3 - t0;
        t1 = y3 + y3;
        y3 = t1 + y3;
        t1 = t0 + t0;
        t0 = t1 + t0;
        t0 = t0 - t2;
        t1 = t4 * y3;
        t2 = t0 * y3;
        y3 = x3 * z3;
        y3 = y3 + t2;
        x3 = t3 * x3;
        x3 = x3 - t1;
        z3 = t4 * z3;
        t1 = t3 * t0;
        z3 = z3 + t1;
        Point {
            x: x3,
            y: y3,
            z: z3,
        }
    }

    /// RCB Alg. 5, the mixed addition: Alg. 4 with Z2 = 1, so two of its
    /// products are gone, 11M + 2 m_b. O has no Z = 1 form, so an O addend
    /// takes a branch; O as the accumulator needs none.
    pub fn add_affine(&self, p: &Point<F>, q: &Affine<F>) -> Point<F> {
        if q.is_identity() {
            return *p;
        }
        let b = self.b;
        let (x1, y1, z1, x2, y2) = (p.x, p.y, p.z, q.x, q.y);
        let mut t0 = x1 * x2;
        let mut t1 = y1 * y2;
        let mut t3 = x2 + y2;
        let mut t4 = x1 + y1;
        t3 = t3 * t4;
        t4 = t0 + t1;
        t3 = t3 - t4;
        t4 = y2 * z1;
        t4 = t4 + y1;
        let mut y3 = x2 * z1;
        y3 = y3 + x1;
        let mut z3 = b * z1;
        let mut x3 = y3 - z3;
        z3 = x3 + x3;
        x3 = x3 + z3;
        z3 = t1 - x3;
        x3 = t1 + x3;
        y3 = b * y3;
        t1 = z1 + z1;
        let mut t2 = t1 + z1;
        y3 = y3 - t2;
        y3 = y3 - t0;
        t1 = y3 + y3;
        y3 = t1 + y3;
        t1 = t0 + t0;
        t0 = t1 + t0;
        t0 = t0 - t2;
        t1 = t4 * y3;
        t2 = t0 * y3;
        y3 = x3 * z3;
        y3 = y3 + t2;
        x3 = t3 * x3;
        x3 = x3 - t1;
        z3 = t4 * z3;
        t1 = t3 * t0;
        z3 = z3 + t1;
        Point {
            x: x3,
            y: y3,
            z: z3,
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
