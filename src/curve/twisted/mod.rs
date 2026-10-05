//! Complete a = -1 twisted Edwards curves over an odd field F in which -1
//! is a square, generic over `field::Packed`, with 16-byte encodings of a
//! quotient group.
//!
//! - Curve: E_d : -u^2 + v^2 = 1 + d u^2 v^2, d a random non-square per
//!   namespace. -1 a square (q = 1 mod 4) makes this complete.
//! - Group: G = E(F)/⟨T⟩, T = (0, -1), of order #E/2: 2r where a
//!   certified d has #E = 4r, r prime. Sums are computed on E; equality
//!   and encodings are those of G.
//! - Accumulators (`Point`): extended (X : Y : Z : T), 7M per add of a
//!   `Cached` addend and 9M per add of another extended point
//!   (Hisil–Wong–Carter–Dawson, unified, so they double too).
//! - Addends (`Cached`): (v - u, v + u, 2d u v) from the affine point, 2M
//!   and no inversion.
//!
//! Why the quotient: E's points need u and v's sign, BITS + 1 bits, but
//! (u, v) + T = (-u, -v), so a class of G has one representative with
//! sgn0(v) = 0 (sgn0(u) = 0 when v = 0, the class of (±i, 0)), and its
//! packed u (BITS <= 128 bits) identifies the class. A collision in G
//! projects to one in the order-r part as a collision in E would, so
//! nothing is lost by hashing to E and not clearing the cofactor.
//!
//! 2-torsion: -d is a non-square, so the points at infinity aren't
//! rational and T is E's only point of order 2; E's 2-part is cyclic, and
//! (±i, 0), i^2 = -1, have order 4. So 4 | #E always, and 8 | #E iff some
//! point has order 8.

mod affine;
mod curve;
mod extended;
mod hash_to_curve;
mod impls;

pub use affine::Affine;
pub use curve::Curve;
pub use extended::{Cached, Point};
pub use hash_to_curve::MontgomeryModel;

/// Encodings: the canonical representative's packed u, little-endian.
pub const BYTES: usize = 16;

#[cfg(test)]
pub(crate) mod tests;
