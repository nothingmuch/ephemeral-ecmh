//! Ordinary binary curves E_B : y^2 + xy = x^3 + a x^2 + B over GF(2^m),
//! B random per namespace, generic over the family's `Model`: the field, a
//! with Tr(a) = 1, the type B is drawn from, and the encoding's layout.
//!
//! - Points (`Affine`): x = 0 stands for O, since the 2-torsion point
//!   (0, sqrt B) is outside E\[r\]. The affine law never mentions B.
//! - Batch sums: the affine chord law, 5M + 1S per addition and one
//!   inversion shared by the batch.
//! - Encoding: x | Tr(y) << SIGN in `Model::Bytes`, 0 for O; hashing is
//!   try-and-increment on the same x and sign.
//!
//! The families are `binary127` and `binary109`, a = 1 over `gf2_127` and
//! `gf2_109`, and `binary122`'s two, a = u over `gf2_122` with B dense or
//! in GF(2^61).

mod affine;
mod batch;
mod curve;
#[cfg(test)]
pub(crate) mod tests;
pub use affine::Affine;
pub use batch::{add_batch, sum_batch};
pub use curve::{Constant, Curve, Model};
