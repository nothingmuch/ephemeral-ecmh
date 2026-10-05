//! Hashing to E\[r\] by try-and-increment, one at a time or with shared
//! inversions.

use super::{Affine, Curve, Model};
use crate::curve::h2c::{BatchLift, Lift, try_and_increment, try_and_increment_batch};
use crate::field::Binary;
use crate::field::batch::Invert;
use crate::hash::Salted;

/// The trace-1 x and the sign bit that a 128-bit digest half stands for:
/// the bits in `F::MASK`, plus a if that has trace 0, then bit SIGN.
/// Tr(a) = 1, so the correction is two-to-one from m-bit field elements
/// onto the trace-1 hyperplane, and a uniform digest gives a uniform
/// trace-1 x.
pub fn candidate<M: Model>(c: u128) -> (M::F, u32) {
    let x = M::F::new(c);
    // a or 0, without branching on the digest
    let fix = M::A.value() & ((x.trace() ^ 1) as u128).wrapping_neg();
    (x + M::F::new(fix), (c >> M::SIGN) as u32 & 1)
}

/// A digest half as `candidate`'s x and sign. x has trace 1, so it is
/// never zero.
impl<M: Model> Lift<Affine<M>> for Curve<M> {
    type Candidate = (M::F, u32);
    fn candidate(&self, c: u128) -> Option<(M::F, u32)> {
        Some(candidate::<M>(c))
    }
    fn lift(&self, (x, sign): (M::F, u32)) -> Option<Affine<M>> {
        self.decode_with_inverse(x, x.inv(), sign)
    }
}

impl<M: Model> BatchLift<Affine<M>> for Curve<M> {
    type F = M::F;
    fn denominator(&self, (x, _): (M::F, u32)) -> M::F {
        x
    }
    fn lift_inv(&self, (x, sign): (M::F, u32), u: M::F) -> Option<Affine<M>> {
        self.decode_with_inverse(x, u, sign)
    }
}

impl<M: Model> Curve<M> {
    /// Try-and-increment: decode each digest half's `candidate` until one
    /// lands. Uniform on E\[r\] \ {O} for a random-oracle digest; each try
    /// succeeds with probability (r - 1)/2^m, about 1/2, and a failed try
    /// costs I + 1M.
    pub fn hash_to_curve(&self, h: &Salted, msg: &[u8]) -> Affine<M> {
        try_and_increment(self, h, msg)
    }

    /// `hash_to_curve` for many items, with one shared inversion per round.
    /// Each item needs about 2 tries on average; a batch of n items needs
    /// about log2(n) + 1 rounds, and after k rounds some item remains with
    /// probability at most n (1 - p)^k, p = (r - 1)/2^m. Each try then
    /// costs about 4M + trace. Same outputs as the single version.
    pub fn hash_to_curve_batch(&self, h: &Salted, msgs: &[&[u8]]) -> Vec<Affine<M>> {
        try_and_increment_batch(self, h, msgs)
    }
}
