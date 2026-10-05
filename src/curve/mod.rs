//! Curve models and their concrete field instances.
//!
//!
//! Constructing a curve does not certify its order. Family-specific
//! `curvegen` verifiers check the stated order policy, subject to their
//! documented primality assumptions.

pub mod binary;
pub mod binary127;
