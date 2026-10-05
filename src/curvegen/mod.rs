//! Deterministic curve candidates and selection certificates from a supplied seed.
//!
//! The application chooses the seed and its namespace. Selection uses explicit
//! order, point, and embedding-degree checks. The order policy uses a fixed-base
//! probable-prime test; it does not supply a primality certificate.
//!
//! - `criteria`: candidate/order contracts and direct selection from point counts.

pub mod criteria;
mod number;
pub mod prove;
pub mod select;
