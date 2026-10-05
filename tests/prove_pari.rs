//! The prover with PARI counting, against kat.sage and the verifier.
#![cfg(feature = "pari")]

#[path = "common/kats.rs"]
mod kats;
#[path = "common/orders.rs"]
mod orders;
use ephemeral_ecmh::curve::{binary127, edwards127, weier127};
use ephemeral_ecmh::curvegen::criteria::{self, Count, Criteria};
use ephemeral_ecmh::curvegen::pari::{self, Pari};
use ephemeral_ecmh::curvegen::prove::{
    self, Edwards127, Factor, Family, NoSieve, Sieve, SmallL, Verdict, Weier127,
};
use ephemeral_ecmh::curvegen::select::{self, Certificate};
use ephemeral_ecmh::curvegen::sieve;
use ephemeral_ecmh::field::fp127::Fp;
use ephemeral_ecmh::field::gf2_127::from_u128;
use ephemeral_ecmh::hash::Salted;
use kats::CertKat;

fn ls(c: &Certificate) -> Vec<u128> {
    c.rejections.iter().map(|e| e.0).collect()
}

fn kat_ls(k: &CertKat) -> Vec<u128> {
    k.rejections.iter().map(|e| e.0).collect()
}

#[test]
fn seadata_is_found() {
    assert!(pari::seadata().is_some(), "set GP_DATA_DIR");
}

/// Counted curves per family to check PARI against. PARI is the
/// reference; what is under test is the binding, which converts every
/// curve of a family the same way, so a few cover it, prime orders and
/// composite ones alike.
const SAMPLE: usize = 8;

#[test]
fn counts_match_sage() {
    for &(v, n) in &orders::GF2_127.orders[..SAMPLE] {
        assert_eq!(Pari.order(&binary127::Curve::new(from_u128(v))), n);
    }
    for &(d, n) in &orders::FP127.orders[..SAMPLE] {
        assert_eq!(Pari.order(&edwards127::Curve::new(Fp::new(d)).unwrap()), n);
    }
    for &(b, n) in &orders::WEIER127.orders[..SAMPLE] {
        assert_eq!(Pari.order(&weier127::Curve::new(Fp::new(b)).unwrap()), n);
    }
}

/// ellsea(E, tors) counts in full, or aborts on a candidate the prover
/// rejects.
#[test]
fn early_abort_is_sound() {
    fn check<F: Family>(n: u128, got: Option<u128>) -> usize {
        match got {
            Some(m) => assert_eq!(m, n),
            None => assert!(!matches!(F::verdict(n), Verdict::Accept(_))),
        }
        usize::from(got.is_none())
    }
    let mut aborted = 0;
    for &(d, n) in &orders::FP127.orders[..SAMPLE] {
        let c = edwards127::Curve::new(Fp::new(d)).unwrap();
        aborted += check::<Edwards127>(n, Pari.order_early_abort(&c, Edwards127::TORS));
    }
    for &(b, n) in &orders::WEIER127.orders[..SAMPLE] {
        let c = weier127::Curve::new(Fp::new(b)).unwrap();
        aborted += check::<Weier127>(n, Pari.order_early_abort(&c, Weier127::TORS));
    }
    assert!(aborted > 0);
}

#[test]
fn factors_match_sage() {
    for o in [&orders::GF2_127, &orders::FP127, &orders::WEIER127] {
        for &(m, l) in o.factors {
            assert_eq!(Pari.smallest_prime_factor(m), Some(l));
        }
    }
}

fn check_kat<C>(
    k: &CertKat,
    prove: fn(&[u8; 32]) -> Certificate,
    find: fn(&[u8; 32]) -> (u32, u128),
    verify: fn(&[u8; 32], &Certificate) -> Result<C, select::Error>,
) {
    let cert = prove(&k.seed);
    assert_eq!((cert.index, cert.r), (k.index, k.r));
    assert_eq!(ls(&cert), kat_ls(k));
    verify(&k.seed, &cert).unwrap();
    assert_eq!(find(&k.seed), (k.index, k.r));
}

fn prove_gf2_127(seed: &[u8; 32]) -> Certificate {
    prove::prove_gf2_127(seed, &mut Pari, &mut Pari, &mut NoSieve).0
}

fn prove_fp127(seed: &[u8; 32]) -> Certificate {
    prove::prove_fp127(seed, &mut Pari, &mut Pari, &mut NoSieve).0
}

fn prove_weier127(seed: &[u8; 32]) -> Certificate {
    prove::prove_weier127(seed, &mut Pari, &mut Pari, &mut NoSieve).0
}

fn find_gf2_127(seed: &[u8; 32]) -> (u32, u128) {
    let f = criteria::find::<prove::Binary127>(seed, &mut Pari, &mut NoSieve);
    (f.index, f.r)
}

fn find_fp127(seed: &[u8; 32]) -> (u32, u128) {
    let f = criteria::find::<prove::Edwards127>(seed, &mut Pari, &mut NoSieve);
    (f.index, f.r)
}

fn find_weier127(seed: &[u8; 32]) -> (u32, u128) {
    let f = criteria::find::<prove::Weier127>(seed, &mut Pari, &mut NoSieve);
    (f.index, f.r)
}

/// The prover is the same code for every seed, and its counts are PARI's
/// (counts_match_sage), so one KAT per family shows it reproduces kat.sage:
/// the one with the fewest candidates to count. tests/kat.rs checks every
/// KAT certificate with the verifier.
fn cheapest(ks: &[CertKat]) -> &CertKat {
    ks.iter().min_by_key(|k| k.index).unwrap()
}

#[test]
fn gf2_127_kat() {
    let k = cheapest(kats::GF2_127_CERTS);
    check_kat(k, prove_gf2_127, find_gf2_127, select::verify_gf2_127);
}

#[test]
fn fp127_kat() {
    let k = cheapest(kats::FP127_CERTS);
    check_kat(k, prove_fp127, find_fp127, select::verify_fp127);
}

#[test]
fn weier127_kat() {
    let k = cheapest(kats::WEIER127_CERTS);
    check_kat(k, prove_weier127, find_weier127, select::verify_weier127);
}

/// The sieve changes which candidates PARI counts, never the selection:
/// `find` and `prove`, each at its bound, select the curve that `find`
/// selects unsieved, and the certificate verifies, its odd l up to the
/// bound being sieve rejections. Two seeds per family, picked for short
/// searches.
fn check_sieve<F: Family>(
    seeds: [u32; 2],
    bounds: sieve::Bounds,
    verify: fn(&[u8; 32], &Certificate) -> Result<F::Curve, select::Error>,
) where
    Pari: Count<F::Curve>,
    SmallL: Sieve<F::Curve>,
{
    let h = Salted::new(b"sieve-test", &[0; 32]);
    for i in seeds {
        let seed = h.digest(&[], i);
        let plain = criteria::find::<F>(&seed, &mut Pari, &mut NoSieve);
        let found = criteria::find::<F>(&seed, &mut Pari, &mut SmallL(bounds.find));
        assert_eq!((found.index, found.r), (plain.index, plain.r), "seed {i}");
        let mut sieve = SmallL(bounds.prove);
        let (cert, _) = prove::prove::<F>(&seed, &mut Pari, &mut Pari, &mut sieve);
        assert_eq!((cert.index, cert.r), (plain.index, plain.r), "seed {i}");
        let l_max = bounds.prove.into();
        assert!(
            ls(&cert).iter().any(|&l| l % 2 == 1 && l <= l_max),
            "seed {i}"
        );
        verify(&seed, &cert).unwrap();
    }
}

#[test]
fn edwards127_sieve() {
    check_sieve::<Edwards127>([20, 11], sieve::EDWARDS127, select::verify_fp127);
}

#[test]
fn weier127_sieve() {
    check_sieve::<Weier127>([12, 11], sieve::WEIER127, select::verify_weier127);
}
