//! Rust vs Sage (sage/kat109.sage) for the 14-byte binary variant.

#[path = "common/kats109.rs"]
mod kats109;

use ephemeral_ecmh::curve::binary109::{Affine, Curve, ENCODED_LEN, sum_batch};
use ephemeral_ecmh::curvegen::select109::{Certificate, verify};
use ephemeral_ecmh::hash::Salted;
use kats109::{CertKat, GF2_109_CERTS};
use proptest::prelude::*;

fn bytes(v: u128) -> [u8; ENCODED_LEN] {
    assert!(v >> (8 * ENCODED_LEN) == 0);
    v.to_le_bytes()[..ENCODED_LEN].try_into().unwrap()
}

fn cert(k: &CertKat) -> Certificate {
    Certificate {
        index: k.index,
        r: k.r,
        rejections: k.rejections.iter().map(|&(l, p)| (l, bytes(p))).collect(),
    }
}

fn check_vectors(c: Curve, k: &CertKat) {
    let h = Salted::new(b"kat", &k.seed);
    let msgs: Vec<[u8; 1]> = (0..k.hashes.len() as u8).map(|i| [i]).collect();
    let refs: Vec<&[u8]> = msgs.iter().map(|m| m.as_slice()).collect();
    let pts: Vec<Affine> = refs.iter().map(|m| c.hash_to_curve(&h, m)).collect();
    for (p, want) in pts.iter().zip(k.hashes) {
        assert_eq!(p.encode(), bytes(want));
        assert_eq!(c.decode(bytes(want)), Some(*p));
    }
    assert_eq!(c.hash_to_curve_batch(&h, &refs), pts);
    assert_eq!(sum_batch(&pts).encode(), bytes(k.sum));
    let acc = pts.iter().fold(c.neutral(), |acc, p| c.add_affine(&acc, p));
    assert_eq!(c.to_affine(&acc).encode(), bytes(k.sum));
}

#[test]
fn certificates_and_vectors() {
    for k in GF2_109_CERTS {
        check_vectors(verify(&k.seed, &cert(k)).unwrap(), k);
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
        (
            any::<prop::sample::Index>(),
            0u8..ENCODED_LEN as u8,
            1u8..=255
        )
            .prop_map(|(i, b, x)| Tamper::Point(i, b, x)),
        (0u8..32, 1u8..=255).prop_map(|(b, x)| Tamper::Seed(b, x)),
    ]
}

/// Every tampering must be rejected (or, for the index, land on a
/// different curve and be rejected there).
fn apply(k: &CertKat, t: &Tamper) -> ([u8; 32], Certificate) {
    let (mut seed, mut c) = (k.seed, cert(k));
    match *t {
        Tamper::R(d) => c.r = c.r.wrapping_add(d),
        Tamper::Index(d) => c.index += d,
        Tamper::DropRejection(i) => {
            c.rejections.remove(i.index(c.rejections.len()));
        }
        Tamper::Ell(i) => {
            // l is prime = ord(P), and l + 1 is never a multiple of it
            let j = i.index(c.rejections.len());
            c.rejections[j].0 += 1;
        }
        Tamper::Point(i, b, x) => {
            // Bit 109 alone would give -P, of the same order: move x too.
            let x = if b == 13 && x == 0x20 { 0x21 } else { x };
            let j = i.index(c.rejections.len());
            c.rejections[j].1[b as usize] ^= x;
        }
        Tamper::Seed(b, x) => seed[b as usize] ^= x,
    }
    (seed, c)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn tampered_certificates_fail(i in 0..GF2_109_CERTS.len(), t in tamper()) {
        let (seed, c) = apply(&GF2_109_CERTS[i], &t);
        prop_assert!(verify(&seed, &c).is_err());
    }
}
