//! λ-projective accumulators for `binary`'s group (Oliveira, López, Aranha,
//! Rodríguez-Henríquez, "Lambda coordinates for binary elliptic curves",
//! CHES 2013 / JCEN 2014).
//!
//! - Curve and field: `binary`'s, with the same hashing, encoding and
//!   affine points.
//! - Accumulators (`Point`): (X : L : Z) with x = X/Z, λ = L/Z and
//!   λ = x + y/x. Adding a λ-affine addend (`Affine`) is 8M + 2S, with no
//!   curve constant; preparing one costs an inversion (1 I + 1M), shared
//!   across a batch by `lift_batch`, or none from
//!   [`super::Curve::hash_to_lambda`].
//!
//! Pornin's complete (X : S : Z : T) addition costs 8M + 2S + 2 m_beta,
//! and beta = B^(1/4) is dense on a random curve, so each m_beta is a full
//! multiplication. The group law itself never mentions B, and neither do
//! these formulas: mixed addition with a λ-affine addend is 8M + 2S, with
//! no curve constant. a appears only in
//! - the curve, (λ^2 + λ + a) x^2 = x^4 + B, or
//!   (L^2 + LZ + a Z^2) X^2 = X^4 + B Z^4: substitute y = x(λ + x) into
//!   y^2 + xy = x^3 + a x^2 + B;
//! - doubling, x(2P) = λ^2 + λ + a ([`Affine::double`]);
//!
//! and not in mixed addition, whose a terms cancel between the two
//! points' curve equations.
//!
//! The formulas are not complete: O as either input and x_P = x_Q (P = Q,
//! which needs a doubling, or P = -Q, whose sum is O and routine in RIBLT,
//! where a cell empties) take explicit branches, which is acceptable on
//! public data. The addend needs λ, one inversion per item, where Pornin's
//! lift is 1S + 2 m_beta.

use super::Model;
use crate::curve::h2c::{BatchLift, Lift, try_and_increment, try_and_increment_batch};
use crate::field::batch::Invert;
use crate::field::{Binary, Field};
use crate::group::{Accumulate, Decode, Encode, Group, HashToCurve, Negate, SumBatch};
use crate::hash::Salted;

/// `binary`'s curve under the λ-projective `Group` impl.
#[derive(Clone, Copy, Debug)]
pub struct Curve<M: Model>(pub super::Curve<M>);

/// (X : L : Z) on (L^2 + LZ + a Z^2) X^2 = X^4 + B Z^4; Z = 0 is O.
#[derive(Clone, Copy, Debug)]
pub struct Point<M: Model> {
    pub x: M::F,
    pub l: M::F,
    pub z: M::F,
}

/// A λ-affine point (x, λ) of E\[r\]; x = 0 stands for O, as in `binary::Affine`.
#[derive(Clone, Copy, Debug)]
pub struct Affine<M: Model> {
    pub x: M::F,
    pub l: M::F,
}

impl<M: Model> PartialEq for Affine<M> {
    fn eq(&self, o: &Self) -> bool {
        self.x.equals(o.x) && (self.is_identity() || self.l.equals(o.l))
    }
}
impl<M: Model> Eq for Affine<M> {}

impl<M: Model> Affine<M> {
    pub const IDENTITY: Self = Self {
        x: M::F::ZERO,
        l: M::F::ZERO,
    };

    pub fn is_identity(&self) -> bool {
        self.x.is_zero()
    }

    /// -(x, y) = (x, y + x), so λ gains y/x's + 1.
    pub fn neg(&self) -> Self {
        Self {
            x: self.x,
            l: self.l + M::F::ONE,
        }
    }

    /// λ = x + y/x: 1 I + 1M.
    pub fn lift(p: &super::Affine<M>) -> Self {
        Self::lift_with(p, p.x.inv())
    }

    fn lift_with(p: &super::Affine<M>, u: M::F) -> Self {
        if p.is_identity() {
            return Self::IDENTITY;
        }
        Self {
            x: p.x,
            l: p.x + p.y * u,
        }
    }

    /// 2P from affine P: 1M + 3S. x(2P) = λ^2 + λ + a is Z, and nonzero:
    /// x(2P) = 0 would make 2P the 2-torsion point (0, sqrt B), and P of
    /// order 4, which E\[r\] lacks. λ(2P) = x_P^2/x(2P) + λ^2 + a + 1, and
    /// λ^2 + a + 1 = Z + λ + 1, so L = X + x_P^2 + (λ + 1) Z with X = Z^2.
    pub fn double(&self) -> Point<M> {
        // O has no λ; the formula would take its placeholder off the curve
        if self.is_identity() {
            return Point::IDENTITY;
        }
        let t = self.l.square() + self.l + M::A;
        let x = t.square();
        Point {
            x,
            l: x + self.x.square() + (self.l + M::F::ONE) * t,
            z: t,
        }
    }
}

/// `Affine::lift` on many points with one shared inversion: 4M each + 1 I.
pub fn lift_batch<M: Model>(ps: &[super::Affine<M>]) -> Vec<Affine<M>> {
    let mut u: Vec<M::F> = ps
        .iter()
        .map(|p| if p.is_identity() { M::F::ONE } else { p.x })
        .collect();
    crate::field::batch::invert(&mut u);
    ps.iter()
        .zip(u)
        .map(|(p, u)| Affine::lift_with(p, u))
        .collect()
}

impl<M: Model> From<Affine<M>> for Point<M> {
    fn from(a: Affine<M>) -> Self {
        if a.is_identity() {
            return Self::IDENTITY;
        }
        Self {
            x: a.x,
            l: a.l,
            z: M::F::ONE,
        }
    }
}

impl<M: Model> Point<M> {
    /// The paper's O; `is_identity` takes any Z = 0 as O.
    pub const IDENTITY: Self = Self {
        x: M::F::ONE,
        l: M::F::ONE,
        z: M::F::ZERO,
    };

    pub fn is_identity(&self) -> bool {
        self.z.is_zero()
    }

    /// Mixed addition, 8M + 2S; P = Q adds `Affine::double`'s 1M + 3S.
    ///
    /// The affine chord law is x3 = x_P x_Q (λ_P + λ_Q) / (x_P + x_Q)^2,
    /// λ3 = x_Q (x3 + x_P)^2 / (x3 x_P) + λ_P + 1. With A = L + λ_Q Z =
    /// Z(λ_P + λ_Q) and B = (X + x_Q Z)^2 = Z^2 (x_P + x_Q)^2, it becomes
    /// X3 = A^2 X x_Q Z, L3 = (A x_Q Z + B)^2 + AB(L + Z), Z3 = ABZ.
    pub fn add_affine(&self, q: &Affine<M>) -> Self {
        if q.is_identity() {
            return *self;
        }
        if self.is_identity() {
            return (*q).into();
        }
        let t = q.x * self.z;
        let e = self.x + t;
        let a = self.l + q.l * self.z;
        // x_P = x_Q: P = Q, where the formulas collapse to (0 : 0 : 0), or
        // P = -Q, where they would give O as (X : X : 0) (B = 0, A = Z);
        // A = 0 iff λ_P = λ_Q tells them apart. Past this, A != 0 as well:
        // with x_P != x_Q it would put P + Q at x3 = 0, the 2-torsion point
        // (0, sqrt B) that E[r] lacks, so Z3 = ABZ != 0.
        if e.is_zero() {
            return if a.is_zero() {
                q.double()
            } else {
                Self::IDENTITY
            };
        }
        let b = e.square();
        let (at, ab) = (a * t, a * b);
        Self {
            x: at * (a * self.x),
            l: (at + b).square() + ab * (self.l + self.z),
            z: ab * self.z,
        }
    }

    /// L1 Z2 = L2 Z1: 2M, no inversion. λ alone fixes a point of
    /// E\[r\] \ {O}: it fixes w, w^2 = λ + 1 + a, and w is injective there
    /// (`wcodec`).
    pub fn equals(&self, o: &Self) -> bool {
        match (self.is_identity(), o.is_identity()) {
            (true, true) => true,
            (false, false) => (self.l * o.z).equals(o.l * self.z),
            _ => false,
        }
    }

    /// `equals` against an affine point: L = l Z, 1M.
    pub fn equals_affine(&self, a: &Affine<M>) -> bool {
        match (self.is_identity(), a.is_identity()) {
            (true, true) => true,
            (false, false) => self.l.equals(a.l * self.z),
            _ => false,
        }
    }

    /// One inversion (of Z) and 3M: x = X/Z, y = x(λ + x).
    pub fn to_affine(&self) -> super::Affine<M> {
        self.normalize_with(self.z.inv())
    }

    fn normalize_with(&self, u: M::F) -> super::Affine<M> {
        if self.is_identity() {
            return super::Affine::IDENTITY;
        }
        let x = self.x * u;
        super::Affine {
            x,
            y: x * (self.l * u + x),
        }
    }
}

/// `Point::to_affine` on many points with one shared inversion: 6M each + 1 I.
pub fn normalize_batch<M: Model>(ps: &[Point<M>]) -> Vec<super::Affine<M>> {
    let mut u: Vec<M::F> = ps
        .iter()
        .map(|p| if p.is_identity() { M::F::ONE } else { p.z })
        .collect();
    crate::field::batch::invert(&mut u);
    ps.iter().zip(u).map(|(p, u)| p.normalize_with(u)).collect()
}

/// `binary::Curve`'s candidates, decoded straight to λ-affine.
impl<M: Model> Lift<Affine<M>> for super::Curve<M> {
    type Candidate = (M::F, u32);
    fn candidate(&self, c: u128) -> Option<(M::F, u32)> {
        Some(super::candidate::<M>(c))
    }
    fn lift(&self, (x, sign): (M::F, u32)) -> Option<Affine<M>> {
        self.decode_lambda(x, x.inv(), sign)
    }
}

impl<M: Model> BatchLift<Affine<M>> for super::Curve<M> {
    type F = M::F;
    fn denominator(&self, (x, _): (M::F, u32)) -> M::F {
        x
    }
    fn lift_inv(&self, (x, sign): (M::F, u32), u: M::F) -> Option<Affine<M>> {
        self.decode_lambda(x, u, sign)
    }
}

impl<M: Model> super::Curve<M> {
    /// `hash_to_curve`'s point as (x, λ = x + y/x), the λ-affine addend,
    /// without the lift's inversion: the root w = y/x is what the decoding
    /// solves for, and 1M still picks the sign.
    pub fn hash_to_lambda(&self, h: &Salted, msg: &[u8]) -> Affine<M> {
        try_and_increment(self, h, msg)
    }

    /// `hash_to_lambda` with shared inversions.
    pub fn hash_to_lambda_batch(&self, h: &Salted, msgs: &[&[u8]]) -> Vec<Affine<M>> {
        try_and_increment_batch(self, h, msgs)
    }

    /// `decode_with_inverse` to (x, λ): y = wx picks the root by Tr(y),
    /// and λ = x + y/x is x + w or x + w + 1.
    fn decode_lambda(&self, x: M::F, u: M::F, sign: u32) -> Option<Affine<M>> {
        let mut w = self.solve_lambda(x, u)?;
        if (x * w).trace() != sign {
            w += M::F::ONE;
        }
        Some(Affine { x, l: w + x })
    }
}

/// `binary`'s points, hashes and encodings, with λ-projective accumulators:
/// the same digests as `binary::Curve`'s.
impl<M: Model> Group for Curve<M> {
    type Affine = super::Affine<M>;
    type Point = Point<M>;
    fn identity(&self) -> Self::Point {
        Point::IDENTITY
    }
    fn is_identity(&self, p: &Self::Point) -> bool {
        p.is_identity()
    }
    fn to_affine(&self, p: &Self::Point) -> Self::Affine {
        p.to_affine()
    }
    fn to_affine_batch(&self, ps: &[Self::Point]) -> Vec<Self::Affine> {
        normalize_batch(ps)
    }
}

impl<M: Model> Accumulate for Curve<M> {
    /// λ-affine: 8M + 2S per add, but an inversion to prepare, so prepare
    /// in batches (or hash straight to it, `Curve::hash_to_lambda`).
    type Addend = Affine<M>;
    fn prepare(&self, a: &Self::Affine) -> Self::Addend {
        Affine::lift(a)
    }
    fn prepare_batch(&self, a: &[Self::Affine]) -> Vec<Self::Addend> {
        lift_batch(a)
    }
    fn add(&self, p: &Self::Point, a: &Self::Addend) -> Self::Point {
        p.add_affine(a)
    }
    /// Lifting to λ takes an inversion per call; bulk callers prepare.
    fn add_affine(&self, p: &Self::Point, a: &Self::Affine) -> Self::Point {
        p.add_affine(&Affine::lift(a))
    }
}

impl<M: Model> Negate for Curve<M> {
    fn neg(&self, a: &Self::Affine) -> Self::Affine {
        a.neg()
    }
    fn neg_addend(&self, a: &Self::Addend) -> Self::Addend {
        a.neg()
    }
    fn equals_addend(&self, p: &Point<M>, a: &Affine<M>) -> bool {
        p.equals_affine(a)
    }
}

impl<M: Model> HashToCurve for Curve<M> {
    fn hash(&self, h: &Salted, msg: &[u8]) -> Self::Affine {
        self.0.hash_to_curve(h, msg)
    }
    fn hash_batch(&self, h: &Salted, msgs: &[&[u8]]) -> Vec<Self::Affine> {
        self.0.hash_to_curve_batch(h, msgs)
    }
}

impl<M: Model> SumBatch for Curve<M> {
    fn sum_batch(&self, a: &[Self::Affine]) -> Self::Affine {
        super::sum_batch(a)
    }
}

impl<M: Model> Encode for Curve<M> {
    type Encoding = M::Bytes;
    fn encode(&self, a: &Self::Affine) -> Self::Encoding {
        a.encode()
    }
}
impl<M: Model> Decode for Curve<M> {
    fn decode(&self, e: &Self::Encoding) -> Option<Self::Affine> {
        self.0.decode(*e)
    }
    fn decode_batch(&self, es: &[Self::Encoding]) -> Vec<Option<Self::Affine>> {
        self.0.decode_batch(es)
    }
}
