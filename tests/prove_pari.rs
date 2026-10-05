//! The prover with PARI counting, against kat.sage and the verifier.
#![cfg(feature = "pari")]

#[path = "common/kats.rs"]
mod kats;
#[path = "common/kats109.rs"]
mod kats109;
#[path = "common/kats122.rs"]
mod kats122;
#[path = "common/orders.rs"]
mod orders;
use ephemeral_ecmh::curve::{binary127, edwards127, weier127};
use ephemeral_ecmh::curvegen::agm::Agm;
use ephemeral_ecmh::curvegen::criteria::{self, Count, Criteria};
use ephemeral_ecmh::curvegen::pari::{self, Pari};
use ephemeral_ecmh::curvegen::prove::{
    self, Binary109, Dense122, Edwards127, Factor, Family, Gls122, NoFactor, NoSieve, Sieve,
    SmallL, Verdict, Weier127,
};
use ephemeral_ecmh::curvegen::select::{self, Certificate};
use ephemeral_ecmh::curvegen::{select109, select122, sieve};
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

/// Every candidate up to a Sage certificate's index, counted: the
/// accepted one has #E = 2r, and each rejected one an order its label
/// accounts for (an odd l dividing #E, or #E itself), with #E/2 not prime.
/// Label 0 marks an index that is no candidate.
fn check_binary_kat<C>(
    seed: &[u8; 32],
    index: u32,
    r: u128,
    rejections: &[(u128, u128)],
    candidate: fn(&[u8; 32], u32) -> Option<C>,
) where
    Pari: Count<C>,
{
    for j in 0..=index {
        let Some(c) = candidate(seed, j) else {
            assert_eq!(rejections[j as usize].0, 0, "j {j}");
            continue;
        };
        let n = Pari.order(&c);
        if j == index {
            assert_eq!(n, 2 * r);
            continue;
        }
        let l = rejections[j as usize].0;
        assert!(!select::is_prime(n / 2), "j {j}");
        assert!(if l % 2 == 1 { n % l == 0 } else { n == l }, "j {j}");
    }
}

/// binary109 with PARI's factors, and counts by PARI (on the cheapest
/// KAT) or the AGM (on all): sage/kat109.sage's certificates to the byte,
/// which the verifier takes.
#[test]
fn gf2_109_kats() {
    let check = |k: &kats109::CertKat, cert: Certificate| {
        assert_eq!((cert.index, cert.r), (k.index, k.r));
        let cert = prove::certificate109(cert);
        let kat: Vec<_> = k
            .rejections
            .iter()
            .map(|&(l, p)| (l, p.to_le_bytes()[..14].try_into().unwrap()))
            .collect();
        assert_eq!(cert.rejections, kat);
        select109::verify(&k.seed, &cert).unwrap();
    };
    let k = kats109::GF2_109_CERTS
        .iter()
        .min_by_key(|k| k.index)
        .unwrap();
    check(
        k,
        prove::prove::<Binary109>(&k.seed, &mut Pari, &mut Pari, &mut NoSieve).0,
    );
    let f = criteria::find::<Binary109>(&k.seed, &mut Pari, &mut NoSieve);
    assert_eq!((f.index, f.r), (k.index, k.r));
    for k in kats109::GF2_109_CERTS {
        check(
            k,
            prove::prove::<Binary109>(&k.seed, &mut Agm, &mut Pari, &mut NoSieve).0,
        );
    }
}

/// The 122-bit families with PARI's counts and no factoring.
#[test]
fn gf2_122_kats() {
    let cheapest = |ks: &'static [kats122::CertKat]| ks.iter().min_by_key(|k| k.index).unwrap();
    let k = cheapest(kats122::GF2_122_CERTS);
    let (cert, _) = prove::prove::<Dense122>(&k.seed, &mut Pari, &mut NoFactor, &mut NoSieve);
    assert_eq!((cert.index, cert.r), (k.index, k.r));
    select122::verify_dense(&k.seed, &cert).unwrap();
    let k = cheapest(kats122::GF2_122_GLS_CERTS);
    let (cert, _) = prove::prove::<Gls122>(&k.seed, &mut Pari, &mut NoFactor, &mut NoSieve);
    assert_eq!((cert.index, cert.r), (k.index, k.r));
    select122::verify_gls(&k.seed, &cert).unwrap();
}

/// A certificate is a function of the seed: proving twice, and counting
/// by the AGM or by PARI under the same factor backend, give the same
/// index, r, labels and encoded witnesses.
fn same_certificates<F: Family>(seed: &[u8; 32], factor: &mut impl Factor)
where
    Agm: Count<F::Curve>,
    Pari: Count<F::Curve>,
{
    let parts = |c: Certificate| (c.index, c.r, c.rejections);
    let agm = parts(prove::prove::<F>(seed, &mut Agm, factor, &mut NoSieve).0);
    let again = parts(prove::prove::<F>(seed, &mut Agm, factor, &mut NoSieve).0);
    let pari = parts(prove::prove::<F>(seed, &mut Pari, factor, &mut NoSieve).0);
    assert!(!agm.2.is_empty(), "no witnesses to compare");
    assert_eq!(again, agm);
    assert_eq!(pari, agm);
}

/// On each family's cheapest KAT seed with a rejection.
#[test]
fn gf2_certificates_are_deterministic_and_counter_independent() {
    let k = kats109::GF2_109_CERTS
        .iter()
        .filter(|k| k.index > 0)
        .min_by_key(|k| k.index)
        .unwrap();
    same_certificates::<Binary109>(&k.seed, &mut Pari);
    let cheapest = |ks: &'static [kats122::CertKat]| {
        let k = ks.iter().filter(|k| k.index > 0).min_by_key(|k| k.index);
        k.unwrap().seed
    };
    same_certificates::<Dense122>(&cheapest(kats122::GF2_122_CERTS), &mut NoFactor);
    same_certificates::<Gls122>(&cheapest(kats122::GF2_122_GLS_CERTS), &mut NoFactor);
}

#[test]
fn counts_match_gf2_109_kats() {
    for k in kats109::GF2_109_CERTS {
        let c = select109::candidate;
        check_binary_kat(&k.seed, k.index, k.r, k.rejections, c);
    }
}

#[test]
fn counts_match_gf2_122_kats() {
    for k in kats122::GF2_122_CERTS {
        let c = select122::dense_candidate;
        check_binary_kat(&k.seed, k.index, k.r, k.rejections, c);
    }
    for k in kats122::GF2_122_GLS_CERTS {
        let c = select122::gls_candidate;
        check_binary_kat(&k.seed, k.index, k.r, k.rejections, c);
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
