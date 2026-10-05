//! Complete Edwards curves over an odd field F, generic over
//! `field::OddField`: E_d : u^2 + v^2 = 1 + d u^2 v^2, d a random
//! non-square per namespace, which makes the addition law complete.
//!
//! - Group: all of E(F), #E = 4r. A collision there projects to one in
//!   the order-r part, so the cofactor is not cleared.
//! - Points (`Affine`) travel in the Montgomery model M_d (below).
//! - Accumulators (`Point`): extended Edwards (X : Y : Z : T). Adding a
//!   prepared addend (`Cached`: Edwards affine with d u v precomputed) is
//!   8M; adding another extended point is 9M + 1 m_d. Preparing converts
//!   from Montgomery affine, an inversion shared across a batch.
//! - Batch sums: the Montgomery affine chord law, 5M + 1S per addition
//!   and one inversion shared by the batch.
//!
//! - Batch sums use the Montgomery model M_d : y^2 = x^3 + a2 x^2 + a4 x,
//!   with a2 = (1+d)/2 and a4 = (1-d)^2/16. That is B y^2 = x^3 + A x^2 + x,
//!   A = 2(1+d)/(1-d), B = 4/(1-d), with x and y scaled by 1/B so the chord
//!   law needs no multiplication by B: x3 = l^2 - a2 - x1 - x2. The
//!   per-namespace parameter only appears in additions (and a4 in tangents).
//! - Streaming accumulators use the Edwards model. The Montgomery x-only
//!   ladder does not apply: ECMH needs general additions.
//! - Here B = A + 2 = 4/(1 - d), and over F_p, p = 3 mod 4, 8 | #E iff
//!   1 - d is a square (2-descent, `curvegen::sieve::edwards_order8`;
//!   sage/kat_common.sage asserts the rule against #E). The twist class
//!   with B square, equivalently B = 1 with A + 2 square, thus has 8 | #E;
//!   cofactor 4, as on Curve1174, needs B = A + 2 non-square.

mod affine;
mod batch;
mod curve;
mod extended;
mod hash_to_curve;
mod impls;
#[cfg(test)]
pub(crate) mod tests;

pub use affine::Affine;
pub use curve::Curve;
pub use extended::{Cached, Point};
