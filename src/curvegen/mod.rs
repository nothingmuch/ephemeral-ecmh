//! Deterministic curve candidates and selection certificates from a supplied seed.
//!
//! The application chooses the seed and its namespace. Selection uses explicit
//! order, point, and embedding-degree checks. The order policy uses a fixed-base
//! probable-prime test; it does not supply a primality certificate.
//!
//! - `criteria`: candidate/order contracts and direct selection from point counts.
//! - `select`, `select107`, `select109`, `select122`, `select128`:
//!   candidates from a seed, and the certificate a verifier checks, for
//!   the 127-, 107-, 109-, 122- and 128-bit families.
//! - `select_fp2`: quadratic-field candidates and order certificates.
//! - `prove`: constructs selection certificates using point counts.
//! - `sieve`: rejects candidates by small torsion, before counting points.
//! - `agm`: counts points on binary curves (canonical lift, with `agm::zq`).
//! - `poly`: the polynomial arithmetic the sieve's division polynomials use.
//! - `pari`: point counting with PARI/GP, behind the `pari` feature.

pub mod agm;
pub mod criteria;
mod number;
#[cfg(feature = "pari")]
pub mod pari;
pub mod poly;
pub mod prove;
pub mod select;
pub mod select107;
pub mod select109;
pub mod select122;
pub mod select128;
pub mod sieve;

pub mod select_fp2;

pub use number::embedding_degree_ok;
