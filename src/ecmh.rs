//! ECMH over a curve implementing the required group capabilities:
//! H(multiset) = sum of H(item).
//!
//! `Ecmh` is the streaming form (one hashed point into an accumulator per
//! insert/remove, as a RIBLT cell does). `digest_batch` is the bulk form:
//! hash and sum through the group's batch operations. Each asks of the
//! curve only the `group` capabilities it uses.

use crate::group::{Encode, Group, HashToCurve, Negate, SumBatch};
use crate::hash::Salted;

pub const TAG_ITEM: &[u8] = b"ephemeral-ecmh/item";

/// Streaming multiset hash under one salt on one curve.
pub struct Ecmh<G: Group> {
    pub group: G,
    pub salt: Salted,
    pub acc: G::Point,
}

impl<G: HashToCurve + Negate + Encode> Ecmh<G> {
    pub fn new(group: G, salt: &[u8; 32]) -> Self {
        Self {
            group,
            salt: Salted::new(TAG_ITEM, salt),
            acc: group.identity(),
        }
    }

    pub fn insert(&mut self, item: &[u8]) {
        let p = self.group.hash(&self.salt, item);
        self.acc = self.group.add_affine(&self.acc, &p);
    }

    pub fn remove(&mut self, item: &[u8]) {
        let p = self.group.neg(&self.group.hash(&self.salt, item));
        self.acc = self.group.add_affine(&self.acc, &p);
    }

    pub fn digest(&self) -> G::Encoding {
        self.group.encode_point(&self.acc)
    }
}

/// Bulk multiset hash using the group's batch operations.
pub fn digest_batch<G: HashToCurve + SumBatch + Encode>(
    group: G,
    salt: &[u8; 32],
    items: &[&[u8]],
) -> G::Encoding {
    let h = Salted::new(TAG_ITEM, salt);
    group.encode(&group.sum_batch(&group.hash_batch(&h, items)))
}
