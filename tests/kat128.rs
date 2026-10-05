//! Rust vs Sage (sage/kat128.sage) for `twisted128`.

#[path = "common/kats128.rs"]
mod kats128;

use ephemeral_ecmh::curve::twisted128::{Affine, BYTES};
use ephemeral_ecmh::curvegen::select128::{Certificate, R_LO, fp128_candidate, verify_fp128};
use ephemeral_ecmh::field::fp128::Fp;
use ephemeral_ecmh::group::{Accumulate, Encode, HashToCurve, SumBatch};
use ephemeral_ecmh::hash::Salted;
use kats128::CertKat;
use proptest::prelude::*;

fn cert(k: &CertKat) -> Certificate {
    Certificate {
        index: k.index,
        r: k.r,
        rejections: k
            .rejections
            .iter()
            .map(|&(l, p)| (l, p.to_le_bytes()))
            .collect(),
    }
}

/// Through the `group` traits, as the checksums use them.
fn check_vectors<G: HashToCurve + Accumulate + SumBatch + Encode<Encoding = [u8; BYTES]>>(
    g: G,
    k: &CertKat,
) {
    let h = Salted::new(b"kat", &k.seed);
    let pts: Vec<G::Affine> = (0..k.hashes.len() as u8)
        .map(|i| g.hash(&h, &[i]))
        .collect();
    for (p, want) in pts.iter().zip(k.hashes) {
        assert_eq!(g.encode(p), want.to_le_bytes());
    }
    assert_eq!(g.encode(&g.sum_batch(&pts)), k.sum.to_le_bytes());
    let acc = pts
        .iter()
        .fold(g.identity(), |acc, p| g.add_affine(&acc, p));
    assert_eq!(g.encode(&g.to_affine(&acc)), k.sum.to_le_bytes());
}

#[test]
fn twisted128_certificates_and_vectors() {
    for k in kats128::FP128_CERTS {
        check_vectors(verify_fp128(&k.seed, &cert(k)).unwrap(), k);
    }
}

#[test]
fn negated_witnesses_verify() {
    for k in kats128::FP128_CERTS {
        let mut c = cert(k);
        let mut w = c.rejections.iter_mut();
        for j in 0..c.index {
            if let Some(g) = fp128_candidate(&k.seed, j) {
                let (_, e) = w.next().unwrap();
                *e = g.decode(*e).unwrap().neg().encode();
            }
        }
        assert!(verify_fp128(&k.seed, &c).is_ok());
    }
}

/// The witness checks' boundaries, on the first certificate's rejections.
#[test]
fn witness_boundaries() {
    let k = &kats128::FP128_CERTS[0];
    let at = |l: u128| k.rejections.iter().position(|w| w.0 == l).unwrap();
    let (eight, three) = (at(8), at(3));
    let with = |j: usize, w: (u128, [u8; BYTES])| {
        let mut c = cert(k);
        c.rejections[j] = w;
        verify_fp128(&k.seed, &c)
    };
    let enc = |j: usize| cert(k).rejections[j].1;
    // the identity and G's order-2 class, (±i, 0), prove nothing
    let order2 = Affine {
        u: Fp::SQRT_M1,
        v: Fp::ZERO,
    }
    .encode();
    for e in [[0; BYTES], order2] {
        assert!(with(eight, (8, e)).is_err());
        assert!(with(three, (3, e)).is_err());
    }
    // a composite odd l is fine while P's order divides it, up to R_LO
    assert_eq!(R_LO % 6, 0);
    assert!(with(three, (9, enc(three))).is_ok());
    assert!(with(three, (R_LO - 3, enc(three))).is_ok());
    assert!(with(three, (R_LO + 3, enc(three))).is_err());
    assert!(with(three, (6, enc(three))).is_err());
}

#[derive(Clone, Debug)]
enum Tamper {
    R(u128),
    Index(u32),
    DropRejection(prop::sample::Index),
    Ell(prop::sample::Index),
    Point(prop::sample::Index, u8, u8),
    Seed(u8, u8),
}

fn tamper() -> impl Strategy<Value = Tamper> {
    prop_oneof![
        (1u128..1 << 20).prop_map(Tamper::R),
        (1u32..4).prop_map(Tamper::Index),
        any::<prop::sample::Index>().prop_map(Tamper::DropRejection),
        any::<prop::sample::Index>().prop_map(Tamper::Ell),
        (any::<prop::sample::Index>(), 0u8..BYTES as u8, 1u8..=255)
            .prop_map(|(i, b, x)| Tamper::Point(i, b, x)),
        (0u8..32, 1u8..=255).prop_map(|(b, x)| Tamper::Seed(b, x)),
    ]
}

/// As in tests/kat.rs: every tampering must be rejected. Witnesses aren't
/// unique (-P, or any other point of order l, would do), so a changed
/// witness byte is rejected with overwhelming probability rather than by
/// construction; `negated_witnesses_verify` pins a valid replacement.
fn apply(k: &CertKat, t: &Tamper) -> ([u8; 32], Certificate) {
    let (mut seed, mut c) = (k.seed, cert(k));
    match *t {
        Tamper::R(d) => c.r = c.r.wrapping_add(d),
        Tamper::Index(d) => c.index += d,
        Tamper::DropRejection(i) => {
            if !c.rejections.is_empty() {
                c.rejections.remove(i.index(c.rejections.len()));
            } else {
                c.index += 1;
            }
        }
        Tamper::Ell(i) => {
            if !c.rejections.is_empty() {
                let j = i.index(c.rejections.len());
                c.rejections[j].0 += 1;
            } else {
                c.r ^= 2;
            }
        }
        Tamper::Point(i, b, x) => {
            if !c.rejections.is_empty() {
                let j = i.index(c.rejections.len());
                c.rejections[j].1[b as usize] ^= x;
            } else {
                c.r = c.r.wrapping_sub(2);
            }
        }
        Tamper::Seed(b, x) => seed[b as usize] ^= x,
    }
    (seed, c)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn tampered_twisted128_certificates_fail(i in 0..kats128::FP128_CERTS.len(), t in tamper()) {
        let (seed, c) = apply(&kats128::FP128_CERTS[i], &t);
        prop_assert!(verify_fp128(&seed, &c).is_err());
    }
}
