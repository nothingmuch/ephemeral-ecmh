//! Group capabilities used by the checksum implementations.
//!
//! [`Group`] defines an affine representation and an accumulator representation.
//! The remaining traits provide operations on these representations:
//!
//! - [`Accumulate`]: preparation of addends and addition to accumulators.
//! - [`Negate`]: negation and comparison of an accumulator with an addend.
//! - [`HashToCurve`]: deterministic mapping of message bytes to affine points.
//! - [`SumBatch`]: summation of affine points.
//! - [`Encode`]: canonical, fixed-width point encoding.
//! - [`Decode`]: decoding of canonical point encodings.
//!
//! `ecmh::Ecmh` requires hashing, negation, and encoding; `ecmh::digest_batch`
//! requires hashing, batch summation, and encoding. The RIBLT simulator requires
//! hashing and negation and does not serialize points.

use crate::hash::Salted;

/// A group instance with affine and accumulator point representations.
/// Point and addend arguments must belong to the same group instance.
pub trait Group: Copy {
    type Affine: Copy;
    type Point: Copy;
    fn identity(&self) -> Self::Point;
    /// Test whether an accumulator represents the identity without normalizing it.
    /// The RIBLT simulator combines this result with the cell count to test emptiness.
    fn is_identity(&self, p: &Self::Point) -> bool;
    fn to_affine(&self, p: &Self::Point) -> Self::Affine;
    /// Normalize each accumulator in input order. The default normalizes
    /// separately; implementations may override it to share inversions.
    fn to_affine_batch(&self, ps: &[Self::Point]) -> Vec<Self::Affine> {
        ps.iter().map(|p| self.to_affine(p)).collect()
    }
}

/// Adding points into an accumulator.
pub trait Accumulate: Group {
    /// A prepared point representation for repeated addition to accumulators.
    /// Preparation can be reused across every cell assigned to the same item.
    type Addend: Copy;
    fn prepare(&self, a: &Self::Affine) -> Self::Addend;
    /// Prepare each affine point in input order. Implementations may override
    /// the default elementwise operation to share preparation work.
    fn prepare_batch(&self, a: &[Self::Affine]) -> Vec<Self::Addend> {
        a.iter().map(|a| self.prepare(a)).collect()
    }
    fn add(&self, p: &Self::Point, a: &Self::Addend) -> Self::Point;
    fn add_affine(&self, p: &Self::Point, a: &Self::Affine) -> Self::Point;
}

/// Negating points, affine and prepared, so that sums can be subtracted from.
pub trait Negate: Accumulate {
    fn neg(&self, a: &Self::Affine) -> Self::Affine;
    fn neg_addend(&self, a: &Self::Addend) -> Self::Addend;
    /// Test whether `p` represents the same group element as `a`.
    /// The default tests `p - a == O`; implementations may compare coordinates
    /// directly. The RIBLT simulator uses this test with a signed item addend.
    fn equals_addend(&self, p: &Self::Point, a: &Self::Addend) -> bool {
        self.is_identity(&self.add(p, &self.neg_addend(a)))
    }
}

/// Hashing bytes to an affine point.
pub trait HashToCurve: Group {
    fn hash(&self, h: &Salted, msg: &[u8]) -> Self::Affine;
    /// Hash each message in input order, with the same outputs as `hash`.
    /// Implementations may override the default to share work across messages.
    fn hash_batch(&self, h: &Salted, msgs: &[&[u8]]) -> Vec<Self::Affine> {
        msgs.iter().map(|m| self.hash(h, m)).collect()
    }
}

/// Summing many affine points at once, e.g. with batched inversions.
pub trait SumBatch: Group {
    fn sum_batch(&self, a: &[Self::Affine]) -> Self::Affine;
}

/// A point's canonical encoding.
pub trait Encode: Group {
    /// A point's fixed-width encoding: 16 bytes for the 127-bit families.
    type Encoding: Copy + Eq + core::fmt::Debug + AsRef<[u8]> + AsMut<[u8]>;
    fn encode(&self, a: &Self::Affine) -> Self::Encoding;
    /// Encode an accumulator. The default normalizes with `to_affine`;
    /// implementations may override it to skip the affine form.
    fn encode_point(&self, p: &Self::Point) -> Self::Encoding {
        self.encode(&self.to_affine(p))
    }
    /// Encode accumulators in input order after one `to_affine_batch` call.
    /// Shared normalization depends on the group's batch implementation.
    fn encode_batch(&self, ps: &[Self::Point]) -> Vec<Self::Encoding> {
        let a = self.to_affine_batch(ps);
        a.iter().map(|a| self.encode(a)).collect()
    }
}

/// Decode canonical point encodings for this group instance.
pub trait Decode: Encode {
    /// Return the point encoded by `e`, or `None` if no point's canonical
    /// encoding equals `e`.
    fn decode(&self, e: &Self::Encoding) -> Option<Self::Affine>;
    /// Decode each encoding in input order, with the same outputs as `decode`.
    /// Implementations may override the default to share inversions.
    fn decode_batch(&self, es: &[Self::Encoding]) -> Vec<Option<Self::Affine>> {
        es.iter().map(|e| self.decode(e)).collect()
    }
}
