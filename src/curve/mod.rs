//! Curve models and their concrete field instances.
//!
//! - `binary127`, `binary109`, `binary122`: y^2 + xy = x^3 + a x^2 + B over
//!   GF(2^127) and GF(2^109) with a = 1, and over GF(2^122) =
//!   GF(2^61)\[u\] with a = u and B dense or in GF(2^61) (GLS), sharing
//!   the laws of `binary`, with its λ-projective accumulators
//!   `binary::lambda`.
//! - `edwards127`, `edwards107`: complete Edwards curves over F_p, for
//!   p = 2^127 - 1 and 2^107 - 1, sharing the laws of `edwards`.
//! - `twisted128`: complete a = -1 twisted Edwards curves over F_p,
//!   p = 2^128 - 275, encoding the quotient by their 2-torsion point,
//!   sharing the laws of `twisted`.
//! - `weier127`, `weier107`: short Weierstrass curves over the same two
//!   prime fields, sharing the laws of `weier`. Their selectors require
//!   prime order.
//! - `edwards61x2`, `weier61x2`: both, over GF(p^2), p = 2^61 - 1.
//! - `twisted61x2`, `twisted64x2`, `twisted_goldilocks2`: `twisted128`'s
//!   model and codec over the three quadratic fields, where -1 is always
//!   a square.
//!
//! Constructing a curve does not certify its order. Family-specific
//! `curvegen` verifiers check the stated order policy, subject to their
//! documented primality assumptions.

pub mod binary;
pub mod binary109;
pub mod binary122;
pub mod binary127;
pub mod edwards;
pub mod edwards107;
pub mod edwards127;
pub mod edwards61x2;
pub mod encoding;
pub mod h2c;
pub mod twisted;
pub mod twisted128;
pub mod twisted61x2;
pub mod twisted64x2;
pub mod twisted_goldilocks2;
pub mod weier;
pub mod weier107;
pub mod weier127;
pub mod weier61x2;
