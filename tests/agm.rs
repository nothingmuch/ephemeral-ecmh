//! Rust point counts (src/curvegen/agm/) vs PARI: sage/agm_vectors.sage, the
//! binary certificates of sage/kat.sage, sage/kat109.sage and
//! sage/kat122.sage, and PARI itself (feature `pari`).

#[path = "common/agm_vectors.rs"]
mod agm_vectors;
mod common;
#[path = "common/kats109.rs"]
mod kats109;
#[path = "common/kats122.rs"]
mod kats122;

use common::kats;
use ephemeral_ecmh::curve::binary127::Curve;
use ephemeral_ecmh::curvegen::agm::{
    Counter, Gf109, Modulus, counter, counter61, counter109, counter122,
};
use ephemeral_ecmh::curvegen::select::{gf2_127_candidate, is_prime};
use ephemeral_ecmh::curvegen::{select109, select122};
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
    check_traces::<Gf109>(agm_vectors::GF109);
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

/// As `counts_reproduce_gf2_127_certificates`, for sage/kat109.sage.
#[test]
fn counts_reproduce_gf2_109_certificates() {
    for k in kats109::GF2_109_CERTS {
        for j in 0..=k.index {
            let n = counter109().order(&select109::candidate(&k.seed, j).unwrap());
            if j < k.index {
                let l = k.rejections[j as usize].0;
                assert!(n.is_multiple_of(l) && !is_prime(n / 2), "j {j}");
            } else {
                assert_eq!(n, 2 * k.r);
            }
        }
    }
}

/// Every candidate up to a sage/kat122.sage certificate's index: the
/// accepted one has #E = 2r, each rejected one an order its label accounts
/// for (an odd l dividing #E, or #E itself) with #E/2 not prime, and label
/// 0 marks an index that is no candidate.
fn check_gf2_122_kat(k: &kats122::CertKat, order: impl Fn(u32) -> Option<u128>) {
    for j in 0..=k.index {
        let l = k.rejections.get(j as usize).map(|e| e.0);
        let Some(n) = order(j) else {
            assert_eq!(l, Some(0), "j {j}");
            continue;
        };
        match l {
            None => assert_eq!(n, 2 * k.r),
            Some(l) => {
                assert!(!is_prime(n / 2), "j {j}");
                assert!(if l % 2 == 1 { n % l == 0 } else { n == l }, "j {j}");
            }
        }
    }
}

#[test]
fn counts_reproduce_gf2_122_dense_certificates() {
    for k in kats122::GF2_122_CERTS {
        check_gf2_122_kat(k, |j| {
            select122::dense_candidate(&k.seed, j).map(|c| counter122().order(&c))
        });
    }
}

#[test]
fn counts_reproduce_gf2_122_gls_certificates() {
    for k in kats122::GF2_122_GLS_CERTS {
        check_gf2_122_kat(k, |j| {
            select122::gls_candidate(&k.seed, j).map(|c| counter61().order(&c))
        });
    }
}

/// The counters against PARI on the candidates of a seed no certificate
/// uses: random B, independent of the Sage fixtures.
#[cfg(feature = "pari")]
mod pari {
    use super::*;
    use ephemeral_ecmh::curvegen::criteria::Count;
    use ephemeral_ecmh::curvegen::pari::Pari;
    use ephemeral_ecmh::hash::Salted;

    fn seed() -> [u8; 32] {
        Salted::new(b"agm-vs-pari", &[0; 32]).digest(&[], 0)
    }

    #[test]
    fn gf2_109_matches_pari() {
        for j in 0..32 {
            let c = select109::candidate(&seed(), j).unwrap();
            assert_eq!(counter109().order(&c), Pari.order(&c), "j {j}");
        }
    }

    #[test]
    fn gf2_122_dense_matches_pari() {
        for j in 0..32 {
            let c = select122::dense_candidate(&seed(), j).unwrap();
            assert_eq!(counter122().order(&c), Pari.order(&c), "j {j}");
        }
    }

    /// PARI counts over GF(2^122), so this checks the subfield identity
    /// #E = (q - 1)^2 + t1^2 too.
    #[test]
    fn gf2_122_gls_matches_pari() {
        for j in 0..32 {
            let c = select122::gls_candidate(&seed(), j).unwrap();
            assert_eq!(counter61().order(&c), Pari.order(&c), "j {j}");
        }
    }
}
