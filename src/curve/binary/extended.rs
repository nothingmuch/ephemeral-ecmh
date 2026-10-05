//! Pornin's extended (X : S : Z : T) coordinates: the accumulators, with
//! complete additions for any a of trace 1 (ePrint 2022/1325, §3 and §5.1).
//!
//! Over (x, s) with x = b/x̄ and s = b(ȳ + (1 + a) x̄ + x̄^2)/x̄^2 for an
//! affine (x̄, ȳ), a point is x = beta X/Z, s = beta S/Z^2, T = XZ, and
//! the group is E\[r\] + N with N = (x = 0, s = b) standing for O. a
//! enters as E = a^2 T1 T2 in the addition, a T in the doubling, and
//! (1 + a) x̄ in the lift and its inverse; for a = 1 these are T1 T2, T
//! and 0. Negation (ȳ -> ȳ + x̄ adds b/x̄ = x to s) and equality do not
//! see it.

use super::{Affine, Constant, Curve, Model};
use crate::field::Field;
use crate::field::batch::Invert;
use crate::field::batch::invert as batch_invert;

/// Pornin's (X : S : Z : T). Complete: no exceptional cases.
#[derive(Clone, Copy, Debug)]
pub struct Point<M: Model> {
    pub x: M::F,
    pub s: M::F,
    pub z: M::F,
    pub t: M::F,
}

impl<M: Model> Point<M> {
    pub fn neg(&self) -> Self {
        Self {
            s: self.s + self.t,
            ..*self
        }
    }

    /// w^2 = s/x = S/T, so equal points have S1 T2 = S2 T1 (2M, no inversion).
    pub fn equals(&self, o: &Self) -> bool {
        (self.s * o.t).equals(o.s * self.t)
    }

    /// X = 0 only at N, the group's identity (`normalize_with`).
    pub fn is_identity(&self) -> bool {
        self.x.is_zero()
    }
}

impl<M: Model> Curve<M> {
    /// N = (x = 0, s = b).
    pub fn neutral(&self) -> Point<M> {
        Point {
            x: M::F::ZERO,
            s: self.beta.gf(),
            z: M::F::ONE,
            t: M::F::ZERO,
        }
    }

    /// P -> P + N without inversion:
    /// (beta : beta(y + x^2 + (1 + a) x) : x : beta x), 1S + 2 m_beta.
    pub fn from_affine(&self, p: &Affine<M>) -> Point<M> {
        if p.is_identity() {
            return self.neutral();
        }
        let (s, t) = self.lift(p);
        Point {
            x: self.beta.gf(),
            s,
            z: p.x,
            t,
        }
    }

    /// `from_affine`'s S and T.
    fn lift(&self, p: &Affine<M>) -> (M::F, M::F) {
        let s = self.beta.mul(p.y + p.x.square() + M::mul_1a(p.x));
        (s, self.beta.mul(p.x))
    }

    /// One inversion (of X): x = beta Z/X, y = beta S/X^2 + x^2 + (1 + a) x.
    pub fn to_affine(&self, p: &Point<M>) -> Affine<M> {
        self.normalize_with(p, p.x.inv())
    }

    fn normalize_with(&self, p: &Point<M>, u: M::F) -> Affine<M> {
        if p.x.is_zero() {
            return Affine::IDENTITY;
        }
        let x = self.beta.mul(p.z * u);
        Affine {
            x,
            y: self.beta.mul(p.s * u.square()) + x.square() + M::mul_1a(x),
        }
    }

    pub fn normalize_batch(&self, ps: &[Point<M>]) -> Vec<Affine<M>> {
        let mut u: Vec<M::F> = ps
            .iter()
            .map(|p| if p.x.is_zero() { M::F::ONE } else { p.x })
            .collect();
        batch_invert(&mut u);
        ps.iter()
            .zip(u)
            .map(|(p, u)| self.normalize_with(p, u))
            .collect()
    }

    /// 8M + 2S + 2 m_beta.
    pub fn add(&self, p: &Point<M>, q: &Point<M>) -> Point<M> {
        let (x1x2, s1s2, z1z2, t1t2) = (p.x * q.x, p.s * q.s, p.z * q.z, p.t * q.t);
        self.finish_add(p, q.s + q.t, x1x2, s1s2, z1z2, t1t2)
    }

    /// Extended + affine: 7M + 2S + 3 m_beta, plus 1S + 2 m_beta to lift `a`.
    pub fn add_affine(&self, p: &Point<M>, a: &Affine<M>) -> Point<M> {
        if a.is_identity() {
            return *p;
        }
        let (s2, t2) = self.lift(a);
        let (x1x2, s1s2, z1z2, t1t2) = (self.beta.mul(p.x), p.s * s2, p.z * a.x, p.t * t2);
        self.finish_add(p, s2 + t2, x1x2, s1s2, z1z2, t1t2)
    }

    fn finish_add(
        &self,
        p: &Point<M>,
        st2: M::F,
        x1x2: M::F,
        s1s2: M::F,
        z1z2: M::F,
        t1t2: M::F,
    ) -> Point<M> {
        let d = (p.s + p.t) * st2;
        let e = M::mul_a2(t1t2);
        let (f, g) = (x1x2.square(), z1z2.square());
        let x = d + s1s2;
        let s = self.beta.mul(g * (s1s2 + e) + f * (d + e));
        let z = self.beta.mul(f + g);
        Point { x, s, z, t: x * z }
    }

    /// 2p, 3M + 5S + 2 m_beta and a product by a (ePrint 2022/1325,
    /// Fig. 3), against `add`'s 8M + 2S + 2 m_beta: `add(p, p)` with its
    /// terms collected.
    pub fn double(&self, p: &Point<M>) -> Point<M> {
        let (xx, zz) = (p.x.square(), p.z.square());
        let v = xx + zz;
        let w = v * (p.s + M::mul_a(p.t)) + xx * p.t;
        let x = p.t.square();
        let z = self.beta.mul(v.square());
        Point {
            x,
            s: self.beta.mul(w.square()),
            z,
            t: x * z,
        }
    }

    /// Double-and-add; tests and certificate checks only. Scalars are
    /// public (orders and witnesses), so the loop skips their leading zeros.
    pub fn mul(&self, p: &Point<M>, k: u128) -> Point<M> {
        let mut acc = self.neutral();
        for i in (0..128 - k.leading_zeros()).rev() {
            acc = self.double(&acc);
            if k >> i & 1 == 1 {
                acc = self.add(&acc, p);
            }
        }
        acc
    }
}
