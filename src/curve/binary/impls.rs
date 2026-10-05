//! The `group` traits: the extended accumulators, try-and-increment hashing
//! and the x-and-sign encoding, for every `Model`.

use super::{Affine, Curve, Model, Point, sum_batch};
use crate::group::{Accumulate, Decode, Encode, Group, HashToCurve, Negate, SumBatch};
use crate::hash::Salted;

impl<M: Model> Group for Curve<M> {
    type Affine = Affine<M>;
    type Point = Point<M>;
    fn identity(&self) -> Self::Point {
        self.neutral()
    }
    fn is_identity(&self, p: &Self::Point) -> bool {
        p.is_identity()
    }
    fn to_affine(&self, p: &Self::Point) -> Self::Affine {
        Curve::to_affine(self, p)
    }
    fn to_affine_batch(&self, ps: &[Self::Point]) -> Vec<Self::Affine> {
        self.normalize_batch(ps)
    }
}

impl<M: Model> Accumulate for Curve<M> {
    /// Extended: the lift (1S + 2 m_beta) happens once, not per cell.
    type Addend = Point<M>;
    fn prepare(&self, a: &Self::Affine) -> Self::Addend {
        self.from_affine(a)
    }
    fn add(&self, p: &Self::Point, a: &Self::Addend) -> Self::Point {
        Curve::add(self, p, a)
    }
    fn add_affine(&self, p: &Self::Point, a: &Self::Affine) -> Self::Point {
        Curve::add_affine(self, p, a)
    }
}

impl<M: Model> Negate for Curve<M> {
    fn neg(&self, a: &Self::Affine) -> Self::Affine {
        a.neg()
    }
    fn neg_addend(&self, a: &Self::Addend) -> Self::Addend {
        a.neg()
    }
    /// S1 T2 = S2 T1: 2M.
    fn equals_addend(&self, p: &Self::Point, a: &Self::Addend) -> bool {
        p.equals(a)
    }
}

impl<M: Model> HashToCurve for Curve<M> {
    fn hash(&self, h: &Salted, msg: &[u8]) -> Self::Affine {
        self.hash_to_curve(h, msg)
    }
    fn hash_batch(&self, h: &Salted, msgs: &[&[u8]]) -> Vec<Self::Affine> {
        self.hash_to_curve_batch(h, msgs)
    }
}

impl<M: Model> SumBatch for Curve<M> {
    fn sum_batch(&self, a: &[Self::Affine]) -> Self::Affine {
        sum_batch(a)
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
        Curve::decode(self, *e)
    }
    fn decode_batch(&self, es: &[Self::Encoding]) -> Vec<Option<Self::Affine>> {
        Curve::decode_batch(self, es)
    }
}
