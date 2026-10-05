//! Rust point counts against the PARI-derived vectors from
//! sage/agm_vectors.sage.

#[path = "common/agm_vectors.rs"]
mod agm_vectors;

use ephemeral_ecmh::curve::binary127::Curve;
use ephemeral_ecmh::curvegen::agm::{Counter, Modulus, counter};
use ephemeral_ecmh::field::gf2_127::from_u128;

/// Sage's default pentanomial modulus.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Gf13;

impl Modulus for Gf13 {
    const M: usize = 13;
    const LOW: &'static [usize] = &[4, 3, 1, 0];
}

#[test]
fn gf127_orders_match_pari() {
    for &(b, n) in agm_vectors::GF127 {
        assert_eq!(counter().order(&Curve::new(from_u128(b))), n, "B = {b:#x}");
    }
}

/// #E_B = q + 1 + t, for t the trace of the a2 = 0 twist.
fn check_traces<F: Modulus>(vectors: &[(u128, u128)]) {
    let c = Counter::<F>::new();
    let q1 = (1 << F::M) + 1;
    for &(b, n) in vectors {
        assert_eq!(c.trace(b), n as i128 - q1, "B = {b:#x}");
    }
}

#[test]
fn other_moduli_match_pari() {
    check_traces::<Gf13>(agm_vectors::GF13);
}
