//! Prime-order short Weierstrass curves over a field of characteristic
//! above 3, generic over `field::OddField`: W_b : y^2 = x^3 - 3x + b, b
//! random per namespace, #E = r prime.
//!
//! - Accumulators (`Point`): projective (X : Y : Z) with the
//!   Renes–Costello–Batina addition. The addend is the `Affine` point
//!   itself, by their mixed addition, 11M + 2 m_b; adding two `Point`s is
//!   12M + 2 m_b.
//! - Batch sums: the affine chord law, 5M + 1S per addition and one
//!   inversion shared by the batch. b never appears; a only enters tangents.
//! - `jacobian`: Jacobian accumulators instead, 7M + 4S per addition but
//!   with branches where RCB is complete.
//!
//! secp256k1's own shape (a = 0) has j = 0: every b gives one of six
//! group orders, and its order-6 automorphism speeds up rho by sqrt 6. So
//! per-namespace curves fix a = -3 and hash b instead.
//!
//! Single additions use the Renes–Costello–Batina complete projective formulas
//! (ePrint 2015/1060, a = -3): Alg. 4, 12M + 2 m_b, and Alg. 5 with an
//! affine addend, 11M + 2 m_b. They are complete only without 2-torsion:
//! on `OddCurve`, the group, but not on the rejected candidates that
//! certificate checking also computes on with the raw `Curve` (see
//! `select`, which uses Alg. 4 alone).

mod affine;
mod batch;
mod curve;
mod hash_to_curve;
mod impls;
pub mod jacobian;
mod projective;
mod sswu;
#[cfg(test)]
pub(crate) mod tests;

pub use affine::Affine;
pub use curve::{Curve, OddCurve};
pub use projective::Point;
pub use sswu::{Sswu, SswuField};
