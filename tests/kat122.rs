//! Rust vs Sage (sage/kat122.sage) for the GF(2^122) families: the
//! certificates verify, and on each certified curve the four families
//! (Pornin and λ accumulators over the dense and the GLS constant) hash,
//! add and encode as Sage does.

#[path = "common/kats122.rs"]
mod kats122;

use ephemeral_ecmh::curve::binary::{Affine, Curve, Model, lambda, sum_batch};
use ephemeral_ecmh::curvegen::criteria::OrderError;
use ephemeral_ecmh::curvegen::select::{Certificate, Error};
use ephemeral_ecmh::curvegen::select122::{verify_dense, verify_gls};
use ephemeral_ecmh::ecmh::Ecmh;
use ephemeral_ecmh::group::{Accumulate, Encode, Group};
use ephemeral_ecmh::hash::Salted;
use kats122::{CertKat, GF2_122_CERTS, GF2_122_GLS_CERTS};
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

fn check_vectors<M: Model<Bytes = [u8; 16]>>(c: Curve<M>, k: &CertKat) {
    let h = Salted::new(b"kat", &k.seed);
    let msgs: Vec<[u8; 1]> = (0..k.hashes.len() as u8).map(|i| [i]).collect();
    let refs: Vec<&[u8]> = msgs.iter().map(|m| m.as_slice()).collect();
    let pts: Vec<Affine<M>> = refs.iter().map(|m| c.hash_to_curve(&h, m)).collect();
    for (p, want) in pts.iter().zip(k.hashes) {
        assert_eq!(p.encode(), want.to_le_bytes());
        assert_eq!(c.decode(want.to_le_bytes()), Some(*p));
    }
    assert_eq!(c.hash_to_curve_batch(&h, &refs), pts);
    let lams: Vec<lambda::Affine<M>> = pts.iter().map(lambda::Affine::lift).collect();
    assert_eq!(c.hash_to_lambda_batch(&h, &refs), lams);
    let sum = k.sum.to_le_bytes();
    assert_eq!(sum_batch(&pts).encode(), sum);
    let acc = pts.iter().fold(c.neutral(), |acc, p| c.add_affine(&acc, p));
    assert_eq!(c.to_affine(&acc).encode(), sum);
    let acc = pts
        .iter()
        .fold(c.neutral(), |acc, p| c.add(&acc, &c.from_affine(p)));
    assert_eq!(c.to_affine(&acc).encode(), sum);
    let acc = lams
        .iter()
        .fold(lambda::Point::<M>::IDENTITY, |acc, p| acc.add_affine(p));
    assert_eq!(acc.to_affine().encode(), sum);
    let (p0, p1) = (c.from_affine(&pts[0]), c.from_affine(&pts[1]));
    assert_eq!(
        c.to_affine(&c.add(&p0, &p0)).encode(),
        k.double.to_le_bytes()
    );
    assert_eq!(
        c.to_affine(&c.add(&p0, &p1.neg())).encode(),
        k.diff.to_le_bytes()
    );
    let l0 = lambda::Point::from(lams[0]);
    assert_eq!(
        l0.add_affine(&lams[0]).to_affine().encode(),
        k.double.to_le_bytes()
    );
    assert_eq!(
        l0.add_affine(&lams[1].neg()).to_affine().encode(),
        k.diff.to_le_bytes()
    );
    assert_eq!(pts[0].add(&pts[0]).encode(), k.double.to_le_bytes());
    assert_eq!(pts[0].add(&pts[1].neg()).encode(), k.diff.to_le_bytes());
    // the families' digests of the same multiset agree
    let salt = k.seed;
    let digest = |g: &dyn Fn(&[&[u8]]) -> [u8; 16]| g(&refs);
    let pornin = digest(&|xs| {
        let mut e = Ecmh::new(c, &salt);
        xs.iter().for_each(|x| e.insert(x));
        e.digest()
    });
    let lam = digest(&|xs| {
        let mut e = Ecmh::new(lambda::Curve(c), &salt);
        xs.iter().for_each(|x| e.insert(x));
        e.digest()
    });
    assert_eq!(pornin, lam);
    let g = lambda::Curve(c);
    let adds = g.prepare_batch(&pts);
    let acc = adds.iter().fold(g.identity(), |acc, a| g.add(&acc, a));
    assert_eq!(g.encode(&g.to_affine(&acc)), sum);
}

#[test]
fn dense_certificates_and_vectors() {
    for k in GF2_122_CERTS {
        check_vectors(verify_dense(&k.seed, &cert(k)).unwrap(), k);
    }
}

#[test]
fn gls_certificates_and_vectors() {
    for k in GF2_122_GLS_CERTS {
        check_vectors(verify_gls(&k.seed, &cert(k)).unwrap(), k);
    }
}

#[test]
fn the_families_do_not_cross() {
    // each family's certificate names another curve under the other's tag
    for (d, g) in GF2_122_CERTS.iter().zip(GF2_122_GLS_CERTS) {
        assert!(verify_gls(&d.seed, &cert(d)).is_err());
        assert!(verify_dense(&g.seed, &cert(g)).is_err());
    }
}

#[test]
fn a_small_embedding_degree_is_refused() {
    // r | q^k - 1 for k = 1: 2^122 = 1 mod r for any prime r dividing 2^122 - 1
    let r = 3;
    assert!(!ephemeral_ecmh::curvegen::select::embedding_degree_ok(
        (1u128 << 122) % r,
        r
    ));
    assert_eq!(
        ephemeral_ecmh::curvegen::select122::verify_dense(
            &GF2_122_CERTS[0].seed,
            &Certificate {
                r: 3,
                ..cert(&GF2_122_CERTS[0])
            }
        )
        .err(),
        Some(Error::Order(OrderError::Hasse))
    );
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
            c.rejections.remove(i.index(c.rejections.len()));
        }
        Tamper::Ell(i) => {
            // an odd l is prime = ord(P), and l + 1 is even, below the
            // Hasse interval; an even l is #E, and #E + 1 is odd, above r_min
            let j = i.index(c.rejections.len());
            c.rejections[j].0 += 1;
        }
        Tamper::Point(i, b, x) => {
            // Bit 127 alone would give -P, of the same order: move x too.
            // Only odd l: an order witness (l = #E) takes any point.
            let x = if b == 15 && x == 0x80 { 0x81 } else { x };
            let odd: Vec<usize> = (0..c.rejections.len())
                .filter(|&j| c.rejections[j].0 & 1 == 1)
                .collect();
            let j = odd[i.index(odd.len())];
            c.rejections[j].1[b as usize] ^= x;
        }
        Tamper::Seed(b, x) => seed[b as usize] ^= x,
    }
    (seed, c)
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(64))]

    #[test]
    fn tampered_dense_certificates_fail(i in 0..GF2_122_CERTS.len(), t in tamper()) {
        let (seed, c) = apply(&GF2_122_CERTS[i], &t);
        prop_assert!(verify_dense(&seed, &c).is_err());
    }

    #[test]
    fn tampered_gls_certificates_fail(i in 0..GF2_122_GLS_CERTS.len(), t in tamper()) {
        let (seed, c) = apply(&GF2_122_GLS_CERTS[i], &t);
        prop_assert!(verify_gls(&seed, &c).is_err());
    }
}
