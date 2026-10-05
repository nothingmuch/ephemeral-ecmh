//! `binary`'s group on the wire as one field element w, hashed to and
//! decoded straight into λ-affine addends. w is Pornin's (ePrint 2022/1325,
//! §4.3): his w^2 = s/x of P + N is λ + 1 + a. His decode keeps the root
//! x(P + N) = b/x(P) of the quadratic below, and this one keeps x(P).
//!
//! For P = (x, y) != O in E\[r\], λ = x + y/x = w^2 + 1 + a (w^2 for
//! a = 1), and w determines P. The λ form of the curve,
//! (λ^2 + λ + a) x^2 = x^4 + B, is (x^2 + d x + b)^2 = 0 with
//! d = w^2 + w + a and b = sqrt(B). Its two roots are x(P) and x(P + T),
//! T = (0, b) the 2-torsion point. They differ by d, and Tr(d) = Tr(a) = 1,
//! so exactly one has Tr(x) = 1, which is membership in E\[r\] = 2E.
//!
//! - Wire: w + w0 in the bits of `F::MASK`, 0 for O. w0 is a fixed w that
//!   decodes to nothing. There is no sign bit, since w + 1 is -P: m - 1
//!   bits of group in m, where `binary`'s own encoding spends m + 1.
//! - Decode: d = w^2 + w + a, reject unless Tr(b/d^2) = 0, then x is
//!   d solve(b/d^2) or that + d. One inversion (batchable), 1M + 1 m_b,
//!   3S and a `solve`.
//! - Encode, from λ-projective: λ = L/Z, w = sqrt(λ + 1 + a). One
//!   inversion (batchable), 1M and a square root.
//! - Hash: try-and-increment on w, the bits of a digest half in
//!   `F::MASK`, which is the decode without w0. Each valid w is one point,
//!   so for independent uniform digest halves (the random-oracle
//!   assumption) the hash is uniform on E\[r\] \ {O}, with both signs in
//!   its range.

use super::lambda::{Affine, Point};
use super::{Constant, Model};
use crate::curve::h2c::{BatchLift, Lift, try_and_increment, try_and_increment_batch};
use crate::field::batch::Invert;
use crate::field::{Binary, Field};
use crate::group::{Accumulate, Decode, Encode, Group, HashToCurve, Negate, SumBatch};
use crate::hash::Salted;

/// `binary`'s curve with the w codec and λ-projective accumulators.
#[derive(Clone, Copy, Debug)]
pub struct Curve<M: Model> {
    pub c: super::Curve<M>,
    /// The first of 0, 2, 4, ... that decodes to nothing (w and w + 1
    /// share d, so odd ones add nothing), which leaves 0 on the wire for O.
    pub w0: M::F,
}

pub(super) fn d<M: Model>(w: M::F) -> M::F {
    w.square() + w + M::A
}

/// The w an encoding's integer stands for, if nothing is set outside the mask.
fn w_of<M: Model>(c: u128) -> Option<M::F> {
    (c & !M::F::MASK == 0).then(|| M::F::new(c))
}

impl<M: Model> Curve<M> {
    pub fn new(c: super::Curve<M>) -> Self {
        let w0 = (0..)
            .step_by(2)
            .map(M::F::new)
            .find(|&w| {
                let d = d::<M>(w);
                c.solve_w(w, d, d.square().inv()).is_none()
            })
            .unwrap();
        Self { c, w0 }
    }

    pub fn decode(&self, c: M::Bytes) -> Option<Affine<M>> {
        self.decode_as(c, Affine::IDENTITY, |w, d, v| self.c.solve_w(w, d, v))
    }

    /// `decode` on many encodings, with one shared inversion.
    pub fn decode_batch(&self, cs: &[M::Bytes]) -> Vec<Option<Affine<M>>> {
        self.decode_batch_as(cs, Affine::IDENTITY, |w, d, v| self.c.solve_w(w, d, v))
    }

    /// `decode` into any form of the point: `o` for 0, otherwise
    /// `lift(w, d(w), 1/d(w)^2)`.
    pub(super) fn decode_as<P>(
        &self,
        c: M::Bytes,
        o: P,
        lift: impl Fn(M::F, M::F, M::F) -> Option<P>,
    ) -> Option<P> {
        let c = M::from_bytes(c);
        if c == 0 {
            return Some(o);
        }
        let w = w_of::<M>(c)? + self.w0;
        let d = d::<M>(w);
        lift(w, d, d.square().inv())
    }

    /// `decode_as` on many encodings, with one shared inversion.
    pub(super) fn decode_batch_as<P: Copy>(
        &self,
        cs: &[M::Bytes],
        o: P,
        lift: impl Fn(M::F, M::F, M::F) -> Option<P>,
    ) -> Vec<Option<P>> {
        let ws: Vec<Option<M::F>> = cs
            .iter()
            .map(|&c| {
                let c = M::from_bytes(c);
                (c != 0)
                    .then(|| w_of::<M>(c))
                    .flatten()
                    .map(|w| w + self.w0)
            })
            .collect();
        let ds: Vec<M::F> = ws.iter().map(|w| w.map_or(M::F::ONE, d::<M>)).collect();
        let mut vs: Vec<M::F> = ds.iter().map(|d| d.square()).collect();
        crate::field::batch::invert(&mut vs);
        cs.iter()
            .zip(ws)
            .zip(ds.into_iter().zip(vs))
            .map(|((c, w), (d, v))| match w {
                Some(w) => lift(w, d, v),
                None if M::from_bytes(*c) == 0 => Some(o),
                None => None,
            })
            .collect()
    }

    pub fn encode(&self, a: &Affine<M>) -> M::Bytes {
        if a.is_identity() {
            return M::to_bytes(0);
        }
        M::to_bytes((Binary::sqrt(a.l + M::F::ONE + M::A) + self.w0).value())
    }

    /// Try-and-increment on w, a digest half at a time; each try succeeds
    /// with probability (r - 1)/2^m, about 1/2, for uniform digest halves.
    pub fn hash(&self, h: &Salted, msg: &[u8]) -> Affine<M> {
        try_and_increment(self, h, msg)
    }

    /// `hash` with one shared inversion per round. Same outputs as the
    /// single version.
    pub fn hash_batch(&self, h: &Salted, msgs: &[&[u8]]) -> Vec<Affine<M>> {
        try_and_increment_batch(self, h, msgs)
    }
}

/// A digest half as w (its bits in `F::MASK`) and d(w); Tr(d) = Tr(a) =
/// 1, so d is never zero.
impl<M: Model> Lift<Affine<M>> for Curve<M> {
    type Candidate = (M::F, M::F);
    fn candidate(&self, c: u128) -> Option<(M::F, M::F)> {
        let w = M::F::new(c);
        Some((w, d::<M>(w)))
    }
    fn lift(&self, (w, d): (M::F, M::F)) -> Option<Affine<M>> {
        self.c.solve_w(w, d, d.square().inv())
    }
}

impl<M: Model> BatchLift<Affine<M>> for Curve<M> {
    type F = M::F;
    fn denominator(&self, (_, d): (M::F, M::F)) -> M::F {
        d.square()
    }
    fn lift_inv(&self, (w, d): (M::F, M::F), v: M::F) -> Option<Affine<M>> {
        self.c.solve_w(w, d, v)
    }
}

impl<M: Model> super::Curve<M> {
    /// The point of E\[r\] that w stands for, given d = d(w) and v = 1/d^2,
    /// or None if x^2 + d x = b has no root (Tr(b/d^2) = 1).
    fn solve_w(&self, w: M::F, d: M::F, v: M::F) -> Option<Affine<M>> {
        // λ = w^2 + 1 + a = d + w + 1
        Some(Affine {
            x: self.root_w(d, v)?,
            l: d + w + M::F::ONE,
        })
    }

    /// x(P) for the P that w stands for: the root of x^2 + d x = b with
    /// Tr(x) = 1, given v = 1/d^2. The other root is x + d, and the two
    /// multiply to b.
    pub(super) fn root_w(&self, d: M::F, v: M::F) -> Option<M::F> {
        let e = self.b.mul(v);
        if e.trace() != 0 {
            return None;
        }
        let x = d * e.solve();
        Some(if x.trace() != 1 { x + d } else { x })
    }
}

impl<M: Model> Point<M> {
    /// (X/Z, L/Z): one inversion and 2M.
    pub fn to_lambda(&self) -> Affine<M> {
        self.lambda_with(self.z.inv())
    }

    fn lambda_with(&self, u: M::F) -> Affine<M> {
        if self.is_identity() {
            return Affine::IDENTITY;
        }
        Affine {
            x: self.x * u,
            l: self.l * u,
        }
    }
}

/// `Point::to_lambda` on many points with one shared inversion: 5M each + 1 I.
pub fn to_lambda_batch<M: Model>(ps: &[Point<M>]) -> Vec<Affine<M>> {
    let mut u: Vec<M::F> = ps
        .iter()
        .map(|p| if p.is_identity() { M::F::ONE } else { p.z })
        .collect();
    crate::field::batch::invert(&mut u);
    ps.iter().zip(u).map(|(p, u)| p.lambda_with(u)).collect()
}

impl<M: Model> Group for Curve<M> {
    type Affine = Affine<M>;
    type Point = Point<M>;
    fn identity(&self) -> Point<M> {
        Point::IDENTITY
    }
    fn is_identity(&self, p: &Point<M>) -> bool {
        p.is_identity()
    }
    fn to_affine(&self, p: &Point<M>) -> Affine<M> {
        p.to_lambda()
    }
    fn to_affine_batch(&self, ps: &[Point<M>]) -> Vec<Affine<M>> {
        to_lambda_batch(ps)
    }
}

impl<M: Model> Accumulate for Curve<M> {
    /// Hashes and decodes are already λ-affine: nothing to prepare.
    type Addend = Affine<M>;
    fn prepare(&self, a: &Affine<M>) -> Affine<M> {
        *a
    }
    fn prepare_batch(&self, a: &[Affine<M>]) -> Vec<Affine<M>> {
        a.to_vec()
    }
    fn add(&self, p: &Point<M>, a: &Affine<M>) -> Point<M> {
        p.add_affine(a)
    }
    fn add_affine(&self, p: &Point<M>, a: &Affine<M>) -> Point<M> {
        p.add_affine(a)
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
    fn hash(&self, h: &Salted, msg: &[u8]) -> Affine<M> {
        Curve::hash(self, h, msg)
    }
    fn hash_batch(&self, h: &Salted, msgs: &[&[u8]]) -> Vec<Affine<M>> {
        Curve::hash_batch(self, h, msgs)
    }
}

impl<M: Model> SumBatch for Curve<M> {
    /// Mixed additions into one accumulator, then one inversion.
    fn sum_batch(&self, a: &[Affine<M>]) -> Affine<M> {
        a.iter()
            .fold(Point::IDENTITY, |p, q| p.add_affine(q))
            .to_lambda()
    }
}

impl<M: Model> Encode for Curve<M> {
    type Encoding = M::Bytes;
    fn encode(&self, a: &Affine<M>) -> M::Bytes {
        Curve::encode(self, a)
    }
}
impl<M: Model> Decode for Curve<M> {
    fn decode(&self, e: &Self::Encoding) -> Option<Affine<M>> {
        Curve::decode(self, *e)
    }
    fn decode_batch(&self, es: &[Self::Encoding]) -> Vec<Option<Affine<M>>> {
        Curve::decode_batch(self, es)
    }
}
