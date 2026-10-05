//! The `group` capability traits: extended accumulators, cached addends.

use super::{Affine, BYTES, Cached, Curve, Point};
use crate::field::{OddField, Packed};
use crate::hash::Salted;

impl<F: OddField> crate::group::Group for Curve<F> {
    type Affine = Affine<F>;
    type Point = Point<F>;
    fn identity(&self) -> Point<F> {
        Point::IDENTITY
    }
    fn is_identity(&self, p: &Point<F>) -> bool {
        p.is_identity()
    }
    fn to_affine(&self, p: &Point<F>) -> Affine<F> {
        Curve::<F>::to_affine(self, p)
    }
    fn to_affine_batch(&self, ps: &[Point<F>]) -> Vec<Affine<F>> {
        self.normalize_batch(ps)
    }
}

impl<F: OddField> crate::group::Accumulate for Curve<F> {
    type Addend = Cached<F>;
    fn prepare(&self, a: &Affine<F>) -> Cached<F> {
        self.to_cached(a)
    }
    fn add(&self, p: &Point<F>, a: &Cached<F>) -> Point<F> {
        self.add_cached(p, a)
    }
    fn add_affine(&self, p: &Point<F>, a: &Affine<F>) -> Point<F> {
        self.add_cached(p, &self.to_cached(a))
    }
}

impl<F: OddField> crate::group::Negate for Curve<F> {
    fn neg(&self, a: &Affine<F>) -> Affine<F> {
        a.neg()
    }
    fn neg_addend(&self, a: &Cached<F>) -> Cached<F> {
        a.neg()
    }
    /// (Y + X)(v - u) = (Y - X)(v + u): 2M.
    fn equals_addend(&self, p: &Point<F>, a: &Cached<F>) -> bool {
        p.equals_cached(a)
    }
}

impl<F: Packed> crate::group::HashToCurve for Curve<F> {
    fn hash(&self, h: &Salted, msg: &[u8]) -> Affine<F> {
        self.hash_to_curve(h, msg)
    }
}

impl<F: OddField> crate::group::SumBatch for Curve<F> {
    fn sum_batch(&self, a: &[Affine<F>]) -> Affine<F> {
        Curve::<F>::sum_batch(self, a)
    }
}

impl<F: Packed> crate::group::Encode for Curve<F> {
    type Encoding = [u8; BYTES];
    fn encode(&self, a: &Affine<F>) -> [u8; BYTES] {
        a.encode()
    }
}
impl<F: Packed> crate::group::Decode for Curve<F> {
    fn decode(&self, e: &Self::Encoding) -> Option<Affine<F>> {
        Curve::<F>::decode(self, *e)
    }
}
