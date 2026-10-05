//! The `group` traits on `OddCurve`: the laws' projective accumulators and
//! affine addends, and the hash and encoding of `curve::encoding`, for every
//! `Signed` field.

use super::{Affine, Curve, OddCurve, Point};
use crate::curve::encoding::Signed;
use crate::field::OddField;
use crate::group::{Accumulate, Decode, Encode, Group, HashToCurve, Negate, SumBatch};
use crate::hash::Salted;

impl<F: OddField> Group for OddCurve<F> {
    type Affine = Affine<F>;
    type Point = Point<F>;
    fn identity(&self) -> Point<F> {
        Point::IDENTITY
    }
    /// `OddCurve` excludes RCB's degenerate (0 : 0 : 0), so Z = 0 suffices.
    fn is_identity(&self, p: &Point<F>) -> bool {
        p.z.is_zero()
    }
    fn to_affine(&self, p: &Point<F>) -> Affine<F> {
        Curve::to_affine(&self.curve(), p)
    }
    fn to_affine_batch(&self, ps: &[Point<F>]) -> Vec<Affine<F>> {
        self.curve().normalize_batch(ps)
    }
}

impl<F: OddField> Accumulate for OddCurve<F> {
    type Addend = Affine<F>;
    fn prepare(&self, a: &Affine<F>) -> Affine<F> {
        *a
    }
    fn add(&self, p: &Point<F>, a: &Affine<F>) -> Point<F> {
        Curve::add_affine(&self.curve(), p, a)
    }
    fn add_affine(&self, p: &Point<F>, a: &Affine<F>) -> Point<F> {
        Curve::add_affine(&self.curve(), p, a)
    }
}

impl<F: OddField> Negate for OddCurve<F> {
    fn neg(&self, a: &Affine<F>) -> Affine<F> {
        a.neg()
    }
    fn neg_addend(&self, a: &Affine<F>) -> Affine<F> {
        a.neg()
    }
    /// X = x Z and Y = y Z: 2M. At Z = 0 (where Y != 0 on odd-order
    /// curves) that fails, so only an O addend needs a case of its own.
    fn equals_addend(&self, p: &Point<F>, a: &Affine<F>) -> bool {
        if a.is_identity() {
            return p.is_identity();
        }
        p.x == a.x * p.z && p.y == a.y * p.z
    }
}

impl<F: OddField> SumBatch for OddCurve<F> {
    fn sum_batch(&self, a: &[Affine<F>]) -> Affine<F> {
        Curve::sum_batch(&self.curve(), a)
    }
}

impl<F: OddField + Signed> HashToCurve for OddCurve<F> {
    fn hash(&self, h: &Salted, msg: &[u8]) -> Affine<F> {
        self.curve().hash_to_curve(h, msg)
    }
}

impl<F: OddField + Signed> Encode for OddCurve<F> {
    type Encoding = F::Bytes;
    fn encode(&self, a: &Affine<F>) -> F::Bytes {
        a.encode()
    }
}

impl<F: OddField + Signed> Decode for OddCurve<F> {
    fn decode(&self, e: &F::Bytes) -> Option<Affine<F>> {
        Curve::decode(&self.curve(), *e)
    }
}
