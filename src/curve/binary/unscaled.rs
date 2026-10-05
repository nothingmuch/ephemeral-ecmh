//! Pornin's complete addition with no curve constant in it: `extended`'s
//! coordinates unscaled, with the constant moved into the addends, and
//! the group on the wire as `wcodec`'s w.
//!
//! - Accumulators (`Point`): (X : S : Z), x = X/Z and s = S/Z^2 the
//!   (x, s) of P + N (x = b/x̄, s = b(ȳ + (1 + a) x̄ + x̄^2)/x̄^2 for an
//!   affine (x̄, ȳ)), so N = (0 : b : 1). `extended` keeps x/beta and
//!   s/beta, which puts beta on its S3 and Z3: 2 m_beta per add.
//! - Addends (`Affine`): (u, v) = (x/b, s/b) of P + N, that is u = 1/x̄
//!   and v = (λ + 1 + a)/x̄ = w^2 u, and (0, 1) for O. Pornin's mixed
//!   addition (ePrint 2022/1325, §5.1) on these is his complete formula
//!   with beta cancelled: 7M + 2S, a entering only as a^2 (`Model::mul_a2`).
//!   T = XZ, which his formula carries as a coordinate, is computed here
//!   at 1M rather than stored, so a cell is three elements, as `lambda`'s.
//! - Negation: (u, v + u), as s gains x. Equality: w^2 = s/x = S/(XZ) =
//!   v/u, and w fixes the point (`wcodec`).
//! - Wire, hash and decode: `wcodec`'s w, landing on (u, v) directly. The
//!   two roots of x^2 + d x = b multiply to b, so u = (x̄ + d)/b and
//!   v = w^2 u: 1M and 1 m_(1/b) past `wcodec`'s own decode.

use super::{Constant, Model, wcodec};
use crate::curve::h2c::{BatchLift, Lift, try_and_increment, try_and_increment_batch};
use crate::field::batch::Invert;
use crate::field::{Binary, Field};
use crate::group::{Accumulate, Decode, Encode, Group, HashToCurve, Negate, SumBatch};
use crate::hash::Salted;

/// `binary`'s curve with unscaled accumulators and the w codec.
#[derive(Clone, Copy, Debug)]
pub struct Curve<M: Model> {
    pub w: wcodec::Curve<M>,
    /// 1/b
    pub b_inv: M::K,
}

/// (X : S : Z), x = X/Z and s = S/Z^2 the (x, s) of P + N; X = 0 is N.
#[derive(Clone, Copy, Debug)]
pub struct Point<M: Model> {
    pub x: M::F,
    pub s: M::F,
    pub z: M::F,
}

/// (u, v) = (1/x̄, w^2/x̄) for a point of E\[r\] other than O, (0, 1) for O.
#[derive(Clone, Copy, Debug)]
pub struct Affine<M: Model> {
    pub u: M::F,
    pub v: M::F,
}

impl<M: Model> PartialEq for Affine<M> {
    fn eq(&self, o: &Self) -> bool {
        self.u.equals(o.u) && self.v.equals(o.v)
    }
}
impl<M: Model> Eq for Affine<M> {}

impl<M: Model> Affine<M> {
    pub const IDENTITY: Self = Self {
        u: M::F::ZERO,
        v: M::F::ONE,
    };

    pub fn is_identity(&self) -> bool {
        self.u.is_zero()
    }

    pub fn neg(&self) -> Self {
        Self {
            u: self.u,
            v: self.v + self.u,
        }
    }

    /// The addend of the affine point p: 1 I + 2M.
    pub fn lift(p: &super::Affine<M>) -> Self {
        Self::lift_with(p, p.x.inv())
    }

    /// `lift` given 1/x̄: v = (x̄ + ȳ/x̄ + 1 + a)/x̄.
    pub fn lift_with(p: &super::Affine<M>, u: M::F) -> Self {
        if p.is_identity() {
            return Self::IDENTITY;
        }
        Self {
            u,
            v: (p.x + p.y * u + M::F::ONE + M::A) * u,
        }
    }
}

impl<M: Model> Point<M> {
    pub fn is_identity(&self) -> bool {
        self.x.is_zero()
    }

    /// Mixed addition, 7M + 2S: Pornin's, with T1 = X1 Z1 and the addend
    /// (u, v) for his (X2 : S2 : 1 : T2) = (u : v : 1 : u).
    pub fn add(&self, a: &Affine<M>) -> Self {
        let t = self.x * self.z;
        let (xu, sv, tu) = (self.x * a.u, self.s * a.v, t * a.u);
        let d = (self.s + t) * (a.u + a.v);
        let e = M::mul_a2(tu);
        let (f, g) = (xu.square(), self.z.square());
        Self {
            x: d + sv,
            s: g * (sv + e) + f * (d + e),
            z: f + g,
        }
    }

    /// The same point: S1 X2 Z2 = S2 X1 Z1 (4M), as w^2 = S/(XZ) and N
    /// alone has X = 0.
    pub fn equals(&self, o: &Self) -> bool {
        (self.s * (o.x * o.z)).equals(o.s * (self.x * self.z))
    }

    /// The point a is: S u = v X Z (3M), as w^2 = v/u.
    pub fn equals_affine(&self, a: &Affine<M>) -> bool {
        (self.s * a.u).equals(a.v * (self.x * self.z))
    }
}

impl<M: Model> Curve<M> {
    pub fn new(c: super::Curve<M>) -> Self {
        Self {
            w: wcodec::Curve::new(c),
            b_inv: c.b.recip(),
        }
    }

    pub fn neutral(&self) -> Point<M> {
        Point {
            x: M::F::ZERO,
            s: self.w.c.b.gf(),
            z: M::F::ONE,
        }
    }

    /// (u, v) = (X/(bZ), S/(bZ^2)) given 1/Z: 3M + 1 m_(1/b).
    fn normalize_with(&self, p: &Point<M>, zi: M::F) -> Affine<M> {
        let c = self.b_inv.mul(zi);
        Affine {
            u: p.x * c,
            v: p.s * c * zi,
        }
    }

    /// The addend w stands for, given d = d(w) and 1/d^2.
    fn lift_w(&self, w: M::F, d: M::F, v: M::F) -> Option<Affine<M>> {
        let u = self.b_inv.mul(self.w.c.root_w(d, v)? + d);
        // w^2 = d + w + a, as d = d(w).
        Some(Affine {
            u,
            v: (d + w + M::A) * u,
        })
    }

    pub fn decode(&self, c: M::Bytes) -> Option<Affine<M>> {
        self.w
            .decode_as(c, Affine::IDENTITY, |w, d, v| self.lift_w(w, d, v))
    }

    /// `decode` on many encodings, with one shared inversion.
    pub fn decode_batch(&self, cs: &[M::Bytes]) -> Vec<Option<Affine<M>>> {
        self.w
            .decode_batch_as(cs, Affine::IDENTITY, |w, d, v| self.lift_w(w, d, v))
    }

    /// The encoding of w = sqrt(w^2), or O's (0) for None.
    fn encode_w2(&self, w2: Option<M::F>) -> M::Bytes {
        M::to_bytes(w2.map_or(0, |w2| (Binary::sqrt(w2) + self.w.w0).value()))
    }

    /// w^2 = v/u: one inversion and a square root.
    pub fn encode(&self, a: &Affine<M>) -> M::Bytes {
        self.encode_w2((!a.is_identity()).then(|| a.v * a.u.inv()))
    }

    /// w^2 = S/(XZ) without the addend, whose 1/Z would be a second
    /// inversion: 1 I + 2M and a square root.
    pub fn encode_point(&self, p: &Point<M>) -> M::Bytes {
        self.encode_w2((!p.is_identity()).then(|| p.s * (p.x * p.z).inv()))
    }
}

impl<M: Model> Lift<Affine<M>> for Curve<M> {
    type Candidate = (M::F, M::F);
    fn candidate(&self, c: u128) -> Option<(M::F, M::F)> {
        let w = M::F::new(c);
        Some((w, wcodec::d::<M>(w)))
    }
    fn lift(&self, (w, d): (M::F, M::F)) -> Option<Affine<M>> {
        self.lift_w(w, d, d.square().inv())
    }
}

impl<M: Model> BatchLift<Affine<M>> for Curve<M> {
    type F = M::F;
    fn denominator(&self, (_, d): (M::F, M::F)) -> M::F {
        d.square()
    }
    fn lift_inv(&self, (w, d): (M::F, M::F), v: M::F) -> Option<Affine<M>> {
        self.lift_w(w, d, v)
    }
}

impl<M: Model> Group for Curve<M> {
    type Affine = Affine<M>;
    type Point = Point<M>;
    fn identity(&self) -> Point<M> {
        self.neutral()
    }
    fn is_identity(&self, p: &Point<M>) -> bool {
        p.is_identity()
    }
    /// One inversion of Z, never 0 as the addition is complete.
    fn to_affine(&self, p: &Point<M>) -> Affine<M> {
        self.normalize_with(p, p.z.inv())
    }
    fn to_affine_batch(&self, ps: &[Point<M>]) -> Vec<Affine<M>> {
        let mut zi: Vec<M::F> = ps.iter().map(|p| p.z).collect();
        crate::field::batch::invert(&mut zi);
        ps.iter()
            .zip(zi)
            .map(|(p, zi)| self.normalize_with(p, zi))
            .collect()
    }
}

impl<M: Model> Accumulate for Curve<M> {
    /// Hashes and decodes are already addends: nothing to prepare.
    type Addend = Affine<M>;
    fn prepare(&self, a: &Affine<M>) -> Affine<M> {
        *a
    }
    fn prepare_batch(&self, a: &[Affine<M>]) -> Vec<Affine<M>> {
        a.to_vec()
    }
    fn add(&self, p: &Point<M>, a: &Affine<M>) -> Point<M> {
        p.add(a)
    }
    fn add_affine(&self, p: &Point<M>, a: &Affine<M>) -> Point<M> {
        p.add(a)
    }
}

impl<M: Model> Negate for Curve<M> {
    fn neg(&self, a: &Affine<M>) -> Affine<M> {
        a.neg()
    }
    fn neg_addend(&self, a: &Affine<M>) -> Affine<M> {
        a.neg()
    }
    fn equals_addend(&self, p: &Point<M>, a: &Affine<M>) -> bool {
        p.equals_affine(a)
    }
}

impl<M: Model> HashToCurve for Curve<M> {
    /// Try-and-increment on w, as `wcodec::Curve::hash`.
    fn hash(&self, h: &Salted, msg: &[u8]) -> Affine<M> {
        try_and_increment(self, h, msg)
    }
    fn hash_batch(&self, h: &Salted, msgs: &[&[u8]]) -> Vec<Affine<M>> {
        try_and_increment_batch(self, h, msgs)
    }
}

impl<M: Model> SumBatch for Curve<M> {
    /// Mixed additions into one accumulator, then one inversion.
    fn sum_batch(&self, a: &[Affine<M>]) -> Affine<M> {
        let p = a.iter().fold(self.neutral(), |p, q| p.add(q));
        self.to_affine(&p)
    }
}

impl<M: Model> Encode for Curve<M> {
    type Encoding = M::Bytes;
    fn encode(&self, a: &Affine<M>) -> M::Bytes {
        Curve::encode(self, a)
    }
    fn encode_point(&self, p: &Point<M>) -> M::Bytes {
        Curve::encode_point(self, p)
    }
    /// w^2 = S/(XZ) straight from the accumulators, one shared inversion.
    fn encode_batch(&self, ps: &[Point<M>]) -> Vec<M::Bytes> {
        let mut t: Vec<M::F> = ps
            .iter()
            .map(|p| {
                if p.is_identity() {
                    M::F::ONE
                } else {
                    p.x * p.z
                }
            })
            .collect();
        crate::field::batch::invert(&mut t);
        ps.iter()
            .zip(t)
            .map(|(p, t)| self.encode_w2((!p.is_identity()).then(|| p.s * t)))
            .collect()
    }
}

impl<M: Model> Decode for Curve<M> {
    fn decode(&self, e: &M::Bytes) -> Option<Affine<M>> {
        Curve::decode(self, *e)
    }
    fn decode_batch(&self, es: &[M::Bytes]) -> Vec<Option<Affine<M>>> {
        Curve::decode_batch(self, es)
    }
}
