//! The torsion sieve against Sage: the KAT certificates (sage/kat.sage),
//! and l | #E on a few hundred candidates (sage/sieve.sage).

mod common;

use common::kats::{self, CertKat};
use common::sieve_vectors as sv;
use ephemeral_ecmh::curve::edwards127;
use ephemeral_ecmh::curvegen::select::{self, Certificate};
use ephemeral_ecmh::curvegen::sieve::{self, Rejection};
use ephemeral_ecmh::field::gf2_127::is_zero;

const L_MAX: u32 = 11;

type Next = fn(&[u8; 32], u32, u32, &mut Vec<Rejection>) -> u32;

/// Rebuilds the KAT certificate the way a prover with the sieve would:
/// sieve rejections, and the KAT's own where the sieve passed a candidate
/// that point counting then rejected. The prover took the smallest l
/// too, so the sieve must reject exactly the candidates whose KAT l is at
/// most L_MAX, with that l, and pass the accepted one; and the result must
/// verify.
fn rebuild(k: &CertKat, has_entry: impl Fn(u32) -> bool, next: Next) -> Certificate {
    let mut entries = k.rejections.iter();
    let kat: Vec<Option<(u128, u128)>> = (0..k.index)
        .map(|j| has_entry(j).then(|| *entries.next().unwrap()))
        .collect();
    let mut rej = Vec::new();
    let mut j = 0;
    loop {
        let before = rej.len();
        let s = next(&k.seed, j, L_MAX, &mut rej);
        assert!(s <= k.index, "sieve rejected the accepted candidate");
        let want: Vec<u128> = kat[j as usize..s as usize]
            .iter()
            .flatten()
            .map(|e| e.0)
            .collect();
        let got: Vec<u128> = rej[before..].iter().map(|e| e.0).collect();
        assert_eq!(got, want, "candidates {j}..{s}");
        if s == k.index {
            break;
        }
        let (l, p) = kat[s as usize].unwrap();
        assert!(l > L_MAX as u128, "sieve passed candidate {s}, l = {l}");
        rej.push((l, p.to_le_bytes()));
        j = s + 1;
    }
    Certificate {
        index: k.index,
        r: k.r,
        rejections: rej,
    }
}

#[test]
fn gf2_127_certificates_from_the_sieve() {
    for k in kats::GF2_127_CERTS {
        let c = rebuild(k, |_| true, sieve::next_gf2_127);
        select::verify_gf2_127(&k.seed, &c).unwrap();
    }
}

#[test]
fn edwards127_certificates_from_the_sieve() {
    for k in kats::FP127_CERTS {
        let c = rebuild(
            k,
            |j| select::fp127_candidate(&k.seed, j).is_some(),
            sieve::next_fp127,
        );
        select::verify_fp127(&k.seed, &c).unwrap();
    }
}

#[test]
fn weier127_certificates_from_the_sieve() {
    for k in kats::WEIER127_CERTS {
        let c = rebuild(
            k,
            |j| select::weier127_candidate(&k.seed, j).is_some(),
            sieve::next_weier127,
        );
        select::verify_weier127(&k.seed, &c).unwrap();
    }
}

/// The odd l the vector tests check on every candidate.
const LS: [u32; 5] = [3, 5, 7, 11, 13];

/// Sage's #E says whether l divides; the sieve must find a point of order
/// l exactly then, and the point must pass select.rs's rejection checks.
fn agree<C>(
    vectors: &[(u32, u128)],
    curve: impl Fn(u32) -> C,
    ls: &[u32],
    torsion: impl Fn(&C, u32) -> Option<[u8; 16]>,
    order_l: impl Fn(&C, [u8; 16], u32) -> bool,
) {
    for &(j, n) in vectors {
        let c = curve(j);
        for &l in ls {
            let got = torsion(&c, l);
            assert_eq!(got.is_some(), n % l as u128 == 0, "j = {j}, l = {l}");
            if let Some(p) = got {
                assert!(order_l(&c, p, l), "j = {j}, l = {l}");
            }
        }
    }
}

#[test]
fn gf2_127_torsion_agrees_with_sage() {
    agree(
        sv::GF2_127,
        |j| select::gf2_127_candidate(&sv::SEED, j).unwrap(),
        &LS,
        sieve::gf2_127_torsion,
        |c, p, l| {
            let p = c.decode(p).unwrap();
            !p.is_identity() && is_zero(c.mul(&c.from_affine(&p), l as u128).x)
        },
    );
}

#[test]
fn edwards127_torsion_agrees_with_sage() {
    let o = edwards127::Point::IDENTITY;
    agree(
        sv::FP127,
        |j| select::fp127_candidate(&sv::SEED, j).unwrap(),
        &[&[8][..], &LS].concat(),
        sieve::edwards_torsion,
        |c, p, l| {
            let p = c.from_affine(&c.decode(p).unwrap());
            // l = 8: of order exactly 8; odd l: P != O
            let below = if l == 8 { 4 } else { 1 };
            c.mul(&p, l as u128).equals(&o) && !c.mul(&p, below).equals(&o)
        },
    );
}

#[test]
fn weier127_torsion_agrees_with_sage() {
    agree(
        sv::WEIER127,
        |j| select::weier127_candidate(&sv::SEED, j).unwrap(),
        &[&[2][..], &LS].concat(),
        sieve::weier_torsion,
        |c, p, l| {
            let p = c.decode(p).unwrap();
            !p.is_identity()
                && if l == 2 {
                    p.y.is_zero()
                } else {
                    c.mul(&c.from_affine(&p), l as u128).is_identity()
                }
        },
    );
}
