//! ristretto255 through the group traits, the checksum group of Buchanan's
//! riblt-ecmh (a fork of yangl1996/riblt), so that the RIBLT workloads run
//! over it in the same harness as the random-curve families.

use curve25519_dalek::ristretto::RistrettoPoint;
use curve25519_dalek::traits::{Identity, IsIdentity};
use ephemeral_ecmh::group::{Accumulate, Group, HashToCurve, Negate};
use ephemeral_ecmh::hash::Salted;

/// curve25519-dalek's `RistrettoPoint` serves as affine point, accumulator
/// and addend: dalek has no public affine form, and its cached
/// (projective Niels) form is private, so every addition converts its
/// right operand.
#[derive(Clone, Copy, Debug)]
pub struct Ristretto255;

impl Group for Ristretto255 {
    type Affine = RistrettoPoint;
    type Point = RistrettoPoint;
    fn identity(&self) -> RistrettoPoint {
        RistrettoPoint::identity()
    }
    fn is_identity(&self, p: &RistrettoPoint) -> bool {
        p.is_identity()
    }
    fn to_affine(&self, p: &RistrettoPoint) -> RistrettoPoint {
        *p
    }
}

impl Accumulate for Ristretto255 {
    type Addend = RistrettoPoint;
    fn prepare(&self, a: &RistrettoPoint) -> RistrettoPoint {
        *a
    }
    fn add(&self, p: &RistrettoPoint, a: &RistrettoPoint) -> RistrettoPoint {
        p + a
    }
    fn add_affine(&self, p: &RistrettoPoint, a: &RistrettoPoint) -> RistrettoPoint {
        p + a
    }
}

impl Negate for Ristretto255 {
    fn neg(&self, a: &RistrettoPoint) -> RistrettoPoint {
        -a
    }
    fn neg_addend(&self, a: &RistrettoPoint) -> RistrettoPoint {
        -a
    }
    /// Ristretto equality, as Go's `Element.Equal` in riblt-ecmh's purity
    /// test: two products per side, no subtraction.
    fn equals_addend(&self, p: &RistrettoPoint, a: &RistrettoPoint) -> bool {
        p == a
    }
}

impl HashToCurve for Ristretto255 {
    /// `from_uniform_bytes` of two salted SHA-256 digests, counters 0 and 1.
    /// riblt-ecmh applies it to SHA-512 of the bare item, one compression
    /// for its 64-byte items; the salted digests take two compressions each
    /// for 64-byte items and one for 36-byte ones, and put the salt into
    /// the hash as for the other families.
    fn hash(&self, h: &Salted, msg: &[u8]) -> RistrettoPoint {
        let mut wide = [0u8; 64];
        wide[..32].copy_from_slice(&h.digest(msg, 0));
        wide[32..].copy_from_slice(&h.digest(msg, 1));
        RistrettoPoint::from_uniform_bytes(&wide)
    }
}
