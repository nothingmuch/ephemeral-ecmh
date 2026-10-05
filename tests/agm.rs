//! Rust point counts against PARI-derived sage/agm_vectors.sage vectors
//! and the binary certificates of sage/kat.sage.

#[path = "common/agm_vectors.rs"]
mod agm_vectors;
mod common;

use common::kats;
use ephemeral_ecmh::curve::binary127::Curve;
use ephemeral_ecmh::curvegen::agm::{Counter, Modulus, counter};
use ephemeral_ecmh::curvegen::select::{gf2_127_candidate, is_prime};
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

/// The Rust counter finds the certified index: each earlier candidate has
/// l | #E and #E/2 composite, and the certified one has #E = 2r.
#[test]
fn counts_reproduce_gf2_127_certificates() {
    for k in kats::GF2_127_CERTS {
        for j in 0..=k.index {
            let n = counter().order(&gf2_127_candidate(&k.seed, j).unwrap());
            if j < k.index {
                let l = k.rejections[j as usize].0;
                assert!(
                    n.is_multiple_of(l) && !is_prime(n / 2),
                    "seed {:?} j {j}",
                    k.seed
                );
            } else {
                assert_eq!(n, 2 * k.r);
            }
        }
    }
}
