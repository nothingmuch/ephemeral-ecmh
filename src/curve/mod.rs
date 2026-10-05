//! Curve models and their concrete field instances.
//!
//! - `binary127`, `binary109`, `binary122`: y^2 + xy = x^3 + a x^2 + B over
//!   GF(2^127) and GF(2^109) with a = 1, and over GF(2^122) =
//!   GF(2^61)\[u\] with a = u and B dense or in GF(2^61) (GLS), sharing
//!   the laws of `binary`, with its λ-projective accumulators
//!   `binary::lambda`.
//! - `edwards127`, `edwards107`: complete Edwards curves over F_p, for
//!   p = 2^127 - 1 and 2^107 - 1, sharing the laws of `edwards`.
//! - `weier127`, `weier107`: short Weierstrass curves over the same two
//!   prime fields, sharing the laws of `weier`. Their selectors require
//!   prime order.
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
pub mod encoding;
pub mod h2c;
pub mod weier;
pub mod weier107;
pub mod weier127;
