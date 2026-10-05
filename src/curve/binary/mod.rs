//! Ordinary binary curves E_B : y^2 + xy = x^3 + a x^2 + B over GF(2^m),
//! B random per namespace, generic over the family's `Model`: the field, a
//! with Tr(a) = 1, the type B is drawn from, and the encoding's layout.
//!
//! - Group: E\[r\] = 2E = { Tr(x) = 1 } + O, #E = 2r.
//! - Points (`Affine`): x = 0 stands for O, since the 2-torsion point
//!   (0, sqrt B) is outside E\[r\]. The affine law never mentions B.
//! - Accumulators (`Point`): Pornin's extended (X : S : Z : T), complete.
//!   Adding a prepared addend (another `Point`) is 8M + 2S + 2 m_beta;
//!   adding an `Affine` is 7M + 2S + 3 m_beta. Preparing costs
//!   1S + 2 m_beta, once per item rather than per cell update. beta =
//!   B^(1/4) is dense on a random curve, so each m_beta is a full
//!   multiplication; for a family whose B lies in a subfield it costs less
//!   (2 M61 for `binary122`'s GLS family).
//! - Batch sums: the affine chord law, 5M + 1S per addition and one
//!   inversion shared by the batch.
//! - `lambda`: the same curves with λ-projective (X : L : Z) accumulators,
//!   8M + 2S per add with no curve constant, but not complete.
//! - Encoding: x | Tr(y) << SIGN in `Model::Bytes`, 0 for O; hashing is
//!   try-and-increment on the same x and sign.
//! - `map`: Pornin's deterministic map, for the families whose modulus
//!   gives it constants (`Pornin`).
//!
//! The families are `binary127` and `binary109`, a = 1 over `gf2_127` and
//! `gf2_109`, and `binary122`'s two, a = u over `gf2_122` with B dense or
//! in GF(2^61).

mod affine;
mod batch;
mod curve;
mod extended;
mod hash_to_curve;
mod impls;
pub mod lambda;
mod map;
#[cfg(test)]
pub(crate) mod tests;
pub use affine::Affine;
pub use batch::{add_batch, sum_batch};
pub use curve::{Constant, Curve, Model};
pub use extended::Point;
pub use hash_to_curve::candidate;
pub use map::{MapState, Pornin};
#[cfg(test)]
pub(crate) use map::{check_map_vectors, map_suite};
