//! Rust vs Sage (sage/kat107.sage) for the 14-byte families.

mod common;

use common::kats107::{self, CertKat};
use ephemeral_ecmh::curve::encoding::Signed;
use ephemeral_ecmh::curvegen::select107::{Certificate, verify_fp107, verify_weier107};
use ephemeral_ecmh::field::fp107::{BYTES, Fp};
use ephemeral_ecmh::group::{Accumulate, Encode, HashToCurve, SumBatch};
use ephemeral_ecmh::hash::Salted;
use proptest::prelude::*;

fn cert(k: &CertKat) -> Certificate {
    Certificate {
        index: k.index,
        r: k.r,
        rejections: k
            .rejections
            .iter()
            .map(|&(l, p)| (l, Fp::to_bytes(p)))
            .collect(),
    }
}

fn check_vectors<G: HashToCurve + Accumulate + SumBatch + Encode<Encoding = [u8; BYTES]>>(
    g: G,
    k: &CertKat,
) {
    let h = Salted::new(b"kat", &k.seed);
    let pts: Vec<G::Affine> = (0..k.hashes.len() as u8)
        .map(|i| g.hash(&h, &[i]))
        .collect();
    for (p, want) in pts.iter().zip(k.hashes) {
        assert_eq!(g.encode(p), Fp::to_bytes(want));
    }
    assert_eq!(g.encode(&g.sum_batch(&pts)), Fp::to_bytes(k.sum));
    let acc = pts
        .iter()
        .fold(g.identity(), |acc, p| g.add_affine(&acc, p));
    assert_eq!(g.encode(&g.to_affine(&acc)), Fp::to_bytes(k.sum));
}

#[test]
fn edwards107_certificates_and_vectors() {
    for k in kats107::FP107_CERTS {
        check_vectors(verify_fp107(&k.seed, &cert(k)).unwrap(), k);
    }
}

#[test]
fn weier107_certificates_and_vectors() {
    for k in kats107::WEIER107_CERTS {
        let c = verify_weier107(&k.seed, &cert(k)).unwrap();
        check_vectors(
            ephemeral_ecmh::curve::weier107::OddCurve::new(c).expect("fixture curve of even order"),
            k,
        );
    }
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

/// As in tests/kat.rs: every tampering must be rejected.
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
                // -P has P's order, so a lone sign flip (bit 107) leaves a valid witness
                let x = if b as usize == BYTES - 1 && x == 0x08 {
                    x | 1
                } else {
                    x
                };
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
    fn tampered_edwards107_certificates_fail(i in 0..kats107::FP107_CERTS.len(), t in tamper()) {
        let (seed, c) = apply(&kats107::FP107_CERTS[i], &t);
        prop_assert!(verify_fp107(&seed, &c).is_err());
    }

    #[test]
    fn tampered_weier107_certificates_fail(i in 0..kats107::WEIER107_CERTS.len(), t in tamper()) {
        let (seed, c) = apply(&kats107::WEIER107_CERTS[i], &t);
        prop_assert!(verify_weier107(&seed, &c).is_err());
    }
}
