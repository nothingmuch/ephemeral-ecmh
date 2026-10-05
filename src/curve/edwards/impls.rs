//! The `group` traits: the laws' extended accumulators and cached addends, and
//! the hash and encoding of `curve::encoding`, for every `Signed` field.

use super::{Affine, Cached, Curve, Point};
use crate::curve::encoding::Signed;
use crate::field::OddField;
use crate::group::{Accumulate, Decode, Encode, Group, HashToCurve, Negate, SumBatch};
use crate::hash::Salted;

impl<F: OddField> Group for Curve<F> {
    type Affine = Affine<F>;
    type Point = Point<F>;
    fn identity(&self) -> Point<F> {
        Point::IDENTITY
    }
    /// X = 0 at O and at the 2-torsion point (0, -1); Y = Z tells them apart.
    fn is_identity(&self, p: &Point<F>) -> bool {
        p.x.is_zero() && p.y == p.z
    }
    fn to_affine(&self, p: &Point<F>) -> Affine<F> {
        Curve::to_affine(self, p)
    }
    fn to_affine_batch(&self, ps: &[Point<F>]) -> Vec<Affine<F>> {
        self.normalize_batch(ps)
    }
}

impl<F: OddField> Accumulate for Curve<F> {
    /// Cached Edwards affine (8M per add, against 9M + 1 m_d). It costs an
    /// inversion, so prepare in batches.
    type Addend = Cached<F>;
    fn prepare(&self, a: &Affine<F>) -> Cached<F> {
        self.to_cached_batch(core::slice::from_ref(a))[0]
    }
    fn prepare_batch(&self, a: &[Affine<F>]) -> Vec<Cached<F>> {
        self.to_cached_batch(a)
    }
    fn add(&self, p: &Point<F>, a: &Cached<F>) -> Point<F> {
        self.add_cached(p, a)
    }
    fn add_affine(&self, p: &Point<F>, a: &Affine<F>) -> Point<F> {
        self.add_ext(p, &self.from_affine(a))
    }
}

impl<F: OddField> Negate for Curve<F> {
    fn neg(&self, a: &Affine<F>) -> Affine<F> {
        a.neg()
    }
    fn neg_addend(&self, a: &Cached<F>) -> Cached<F> {
        a.neg()
    }
    /// X = u Z and Y = v Z: 2M. Z is never 0 on a complete curve.
    fn equals_addend(&self, p: &Point<F>, a: &Cached<F>) -> bool {
        p.x == a.u * p.z && p.y == a.v * p.z
    }
}

impl<F: OddField> SumBatch for Curve<F> {
    fn sum_batch(&self, a: &[Affine<F>]) -> Affine<F> {
        Curve::sum_batch(self, a)
    }
}

impl<F: OddField + Signed> HashToCurve for Curve<F> {
    fn hash(&self, h: &Salted, msg: &[u8]) -> Affine<F> {
        self.hash_to_curve(h, msg)
    }
}

impl<F: OddField + Signed> Encode for Curve<F> {
    type Encoding = F::Bytes;
    fn encode(&self, a: &Affine<F>) -> F::Bytes {
        a.encode()
    }
}

impl<F: OddField + Signed> Decode for Curve<F> {
    fn decode(&self, e: &F::Bytes) -> Option<Affine<F>> {
        Curve::decode(self, *e)
    }
}
