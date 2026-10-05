//! Rust vs Sage (sage/kat.sage): certificates, hash-to-curve, sums.

mod common;

use common::kats::{self, CertKat};
use ephemeral_ecmh::curvegen::select::{
    Certificate, verify_fp127, verify_gf2_127, verify_weier127,
};
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
            .map(|&(l, p)| (l, p.to_le_bytes()))
            .collect(),
    }
}

fn check_vectors<G: HashToCurve + Accumulate + SumBatch + Encode<Encoding = [u8; 16]>>(
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
fn gf2_127_certificates_and_vectors() {
    for k in kats::GF2_127_CERTS {
        check_vectors(verify_gf2_127(&k.seed, &cert(k)).unwrap(), k);
    }
}

#[test]
fn edwards127_certificates_and_vectors() {
    for k in kats::FP127_CERTS {
        check_vectors(verify_fp127(&k.seed, &cert(k)).unwrap(), k);
    }
}

#[test]
fn weier127_certificates_and_vectors() {
    for k in kats::WEIER127_CERTS {
        let c = verify_weier127(&k.seed, &cert(k)).unwrap();
        check_vectors(
            ephemeral_ecmh::curve::weier127::OddCurve::new(c).expect("fixture curve of even order"),
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
        (any::<prop::sample::Index>(), 0u8..16, 1u8..=255)
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
            if !c.rejections.is_empty() {
                c.rejections.remove(i.index(c.rejections.len()));
            } else {
                c.index += 1;
            }
        }
        Tamper::Ell(i) => {
            if !c.rejections.is_empty() {
                let j = i.index(c.rejections.len());
                // l is prime or 8 = ord(P), and l + 1 is never a multiple of
                // it; larger shifts could land on a multiple and stay valid
                c.rejections[j].0 += 1;
            } else {
                c.r ^= 2;
            }
        }
        Tamper::Point(i, b, x) => {
            if !c.rejections.is_empty() {
                let j = i.index(c.rejections.len());
                // -P has P's order, so a lone sign flip leaves a valid witness
                let x = if b == 15 && x == 0x80 { x | 1 } else { x };
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
    fn tampered_gf2_127_certificates_fail(i in 0..kats::GF2_127_CERTS.len(), t in tamper()) {
        let (seed, c) = apply(&kats::GF2_127_CERTS[i], &t);
        prop_assert!(verify_gf2_127(&seed, &c).is_err());
    }

    #[test]
    fn tampered_edwards127_certificates_fail(i in 0..kats::FP127_CERTS.len(), t in tamper()) {
        let (seed, c) = apply(&kats::FP127_CERTS[i], &t);
        prop_assert!(verify_fp127(&seed, &c).is_err());
    }

    #[test]
    fn tampered_weier127_certificates_fail(i in 0..kats::WEIER127_CERTS.len(), t in tamper()) {
        let (seed, c) = apply(&kats::WEIER127_CERTS[i], &t);
        prop_assert!(verify_weier127(&seed, &c).is_err());
    }
}

/// Inputs for the map parity vectors: both values of the sign bit, and
/// both values of the low bits the map overwrites.
fn map_inputs() -> [u128; 16] {
    core::array::from_fn(|i| {
        (i as u128 + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15_f39c_c060_5ced_c835)
    })
}

/// Pornin's map on the first binary127 KAT curve: the recorded encodings of
/// `map_inputs()`, singly and batched.
const MAP: [u128; 16] = [
    0xb82b8938ff82e9605445881d7918d9e1,
    0x3920078f93f8ecc6e0cfb9053a425b8b,
    0x136f318dc37b26a360e14f6bbef45283,
    0x491e020e13483ce9db361212488742b3,
    0x32ed1738635b25749e856a8b2f8573cf,
    0x17408defbba53946af575cdce5c5aa01,
    0x2100b40cd4faca1b1d119d41b4be2a6d,
    0x1bd0e5de4448dba0093f44a28f0173f3,
    0x23c177eea5213f9c4b086650b7db7e63,
    0x59ff26d34a98012d02d65a22b4a592f3,
    0x2ae5ef3377c5dee82daa14782d0ee143,
    0xb546c2e8a69d0629c7c696a18c950a65,
    0x3c6c8fe8237f6ff367384dd5c38b2321,
    0x26c5554afc2054f3d3ca67a2b3906ce5,
    0xdc43b2c0ee81ebb30c711a45a3a4449f,
    0x48dcfd6ef7ffb6c82fe6166e7cb89d8f,
];

/// `hash_to_curve_map2` of the messages [0], ..., [15] on the same curve.
const MAP2: [u128; 16] = [
    0xfc2f4c06f0e49d98583a2852df25f2cd,
    0x12a55cef4911f801c283d0199b908a67,
    0x741fdf51c629876306cd37fe6b89638f,
    0x4ba85b74534313d470d69991e3629309,
    0x6f55c07fa753ae8389fa7bb74432fb93,
    0x51c4dd823feb69ea90c0fc8515576501,
    0x830cb9f4270e7ca82b5c51c6859a4321,
    0x67ab31bf6b33592d8d2cb12b2b96dd43,
    0x52e3be78ea3a6375cc59860be997b34f,
    0xc126976843a37d617931e346bb5d390b,
    0x6a3836025a94d2d0fefaccf341edcb33,
    0xca1be500b3c3e518be4716f936209513,
    0x621bdd2007c479af8a253da4a531d3ed,
    0x29fa73ae6c67b3a7784e18372f289463,
    0xa7def7ebc40b36d3f977ef86dfa04b43,
    0x85733f951b74e823d95b3e741380956f,
];

#[test]
fn gf2_127_pornin_map_is_unchanged() {
    use ephemeral_ecmh::curve::binary127::Point;
    use ephemeral_ecmh::curve::h2c::Map;
    let k = &kats::GF2_127_CERTS[0];
    let c = verify_gf2_127(&k.seed, &cert(k)).unwrap();
    let enc = |p: Point| u128::from_le_bytes(c.to_affine(&p).encode());
    let cs = map_inputs();
    let single: Vec<u128> = cs.iter().map(|&v| enc(c.map(v))).collect();
    let batch: Vec<u128> = c.map_to_curve_batch(&cs).into_iter().map(enc).collect();
    let h = Salted::new(b"kat", &k.seed);
    let map2: Vec<u128> = (0..16u8)
        .map(|i| enc(c.hash_to_curve_map2(&h, &[i])))
        .collect();
    assert_eq!(single, MAP);
    assert_eq!(batch, MAP);
    assert_eq!(map2, MAP2);
}

/// The map's λ-affine and unscaled addends are its points: taken back to
/// (x, y) and encoded, they give `MAP`.
#[test]
fn gf2_127_pornin_map_addends_are_its_points() {
    use ephemeral_ecmh::curve::binary::Model;
    use ephemeral_ecmh::curve::binary::lambda;
    use ephemeral_ecmh::curve::binary::unscaled;
    use ephemeral_ecmh::curve::binary127::{Affine, M127};
    use ephemeral_ecmh::field::batch::Invert;
    use ephemeral_ecmh::field::gf2_127::Gf;
    let k = &kats::GF2_127_CERTS[0];
    let c = verify_gf2_127(&k.seed, &cert(k)).unwrap();
    let u = unscaled::Curve::new(c);
    let xl =
        |a: lambda::Affine<M127>| u128::from_le_bytes(lambda::Point::from(a).to_affine().encode());
    // u = 1/x and v = (x + y/x + 1 + a)/x
    let xy = |a: unscaled::Affine<M127>| {
        let x = a.u.inv();
        let y = (a.v + Gf::ONE) * x.square() + (Gf::ONE + M127::A) * x;
        u128::from_le_bytes(Affine { x, y }.encode())
    };
    let cs = map_inputs();
    let single: Vec<u128> = cs.iter().map(|&v| xy(u.map_to_addend(v))).collect();
    let batch: Vec<u128> = u.map_to_addend_batch(&cs).into_iter().map(xy).collect();
    let single_l: Vec<u128> = cs.iter().map(|&v| xl(c.map_to_lambda(v))).collect();
    let batch_l: Vec<u128> = c.map_to_lambda_batch(&cs).into_iter().map(xl).collect();
    assert_eq!(single, MAP);
    assert_eq!(batch, MAP);
    assert_eq!(single_l, MAP);
    assert_eq!(batch_l, MAP);
}
