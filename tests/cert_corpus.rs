//! Acceptance equivalence for the certificate verifiers: for every family
//! and every KAT certificate, a fixed set of malformed certificates, and the
//! failure reason each verifier returns for it. The outcomes are a recorded
//! snapshot: a change to any of them changes which certificates are accepted
//! or how they fail. Tokens name failure reasons, independent of enum
//! wrappers.

#[path = "common/kats.rs"]
mod kats;

use ephemeral_ecmh::curvegen::select;
use ephemeral_ecmh::curvegen::select::{Certificate, Error};

/// A malformation applied to a valid (seed, certificate).
#[derive(Clone, Copy, Debug)]
enum Tamper {
    /// r + 1: even, or outside the range
    RPlusOne,
    /// r = the field's characteristic (Weierstrass: anomalous)
    RIsQ,
    /// r below the admissible range: 3
    RSmall,
    /// r prime above the range (2^127 + 1 is composite; a known prime above 2^127)
    RHuge,
    /// index + 1: one rejection short
    IndexPlusOne,
    /// index - 1 (if any): one rejection too many
    IndexMinusOne,
    /// drop the last rejection (if any)
    DropLast,
    /// label of the first odd-label rejection + 2 (stays odd)
    OddLabelPlusTwo,
    /// flip the low bit of the first odd-label rejection's point
    OddPointBit0,
    /// first rejection's label set to 1
    LabelOne,
    /// first rejection's point set to the all-zero encoding
    PointZero,
    /// seed byte 0 flipped
    SeedBit,
}

const TAMPERS: [Tamper; 12] = [
    Tamper::RPlusOne,
    Tamper::RIsQ,
    Tamper::RSmall,
    Tamper::RHuge,
    Tamper::IndexPlusOne,
    Tamper::IndexMinusOne,
    Tamper::DropLast,
    Tamper::OddLabelPlusTwo,
    Tamper::OddPointBit0,
    Tamper::LabelOne,
    Tamper::PointZero,
    Tamper::SeedBit,
];

/// A prime above 2^127, so that it is prime and out of every range.
const R_HUGE: u128 = (1u128 << 127) + 29;

/// The recorded file opens with how it is produced; the test compares the
/// whole file, so the header is part of what must not change.
const HEADER: &str = "\
# Verdicts of the certificate verifiers on tampered known-answer certificates,
# one line per (family, certificate, tamper), recorded by tests/cert_corpus.rs:
#   nix develop .#rust -c env CERT_CORPUS_RECORD=1 cargo test --test cert_corpus
# The test fails when the current verifiers disagree with this file, so the
# verdicts are those of any revision whose tests pass.
";

fn apply<const N: usize>(
    seed: &[u8; 32],
    c: &Certificate<N>,
    q: u128,
    t: Tamper,
) -> Option<([u8; 32], Certificate<N>)> {
    let (mut seed, mut c) = (*seed, c.clone());
    let first_odd = c.rejections.iter().position(|w| w.0 & 1 == 1);
    match t {
        Tamper::RPlusOne => c.r += 1,
        Tamper::RIsQ => c.r = q,
        Tamper::RSmall => c.r = 3,
        Tamper::RHuge => c.r = R_HUGE,
        Tamper::IndexPlusOne => c.index += 1,
        Tamper::IndexMinusOne => c.index = c.index.checked_sub(1)?,
        Tamper::DropLast => {
            c.rejections.pop()?;
        }
        Tamper::OddLabelPlusTwo => c.rejections[first_odd?].0 += 2,
        Tamper::OddPointBit0 => c.rejections[first_odd?].1[0] ^= 1,
        Tamper::LabelOne => c.rejections.first_mut()?.0 = 1,
        Tamper::PointZero => c.rejections.first_mut()?.1 = [0; N],
        Tamper::SeedBit => seed[0] ^= 1,
    }
    Some((seed, c))
}

fn reason(e: Error) -> String {
    match e {
        Error::Order(e) => format!("{e:?}"),
        e => format!("{e:?}"),
    }
}

/// One line per (family, certificate, tamper): the verifier's outcome.
fn outcomes<const N: usize, C>(
    family: &str,
    q: u128,
    certs: impl Iterator<Item = ([u8; 32], Certificate<N>)>,
    verify: impl Fn(&[u8; 32], &Certificate<N>) -> Result<C, Error>,
    out: &mut Vec<String>,
) {
    for (i, (seed, c)) in certs.enumerate() {
        assert!(verify(&seed, &c).is_ok(), "{family} KAT {i} must verify");
        for t in TAMPERS {
            let line = match apply(&seed, &c, q, t) {
                None => format!("{family} {i} {t:?}: n/a"),
                Some((seed, c)) => match verify(&seed, &c) {
                    Ok(_) => format!("{family} {i} {t:?}: Ok"),
                    Err(e) => format!("{family} {i} {t:?}: {}", reason(e)),
                },
            };
            out.push(line);
        }
    }
}

fn cert16(index: u32, r: u128, rejections: &[(u128, u128)]) -> Certificate<16> {
    Certificate {
        index,
        r,
        rejections: rejections
            .iter()
            .map(|&(l, p)| (l, p.to_le_bytes()))
            .collect(),
    }
}

const P127: u128 = (1 << 127) - 1;
const Q127: u128 = 1 << 127;

#[test]
fn malformed_certificates_keep_their_outcomes() {
    let mut out = Vec::new();
    macro_rules! family {
        ($name:literal, $q:expr, $certs:expr, $cert:ident, $verify:expr) => {
            outcomes(
                $name,
                $q,
                $certs
                    .iter()
                    .map(|k| (k.seed, $cert(k.index, k.r, k.rejections))),
                $verify,
                &mut out,
            );
        };
    }
    family!(
        "gf2_127",
        Q127,
        kats::GF2_127_CERTS,
        cert16,
        select::verify_gf2_127
    );
    family!(
        "edwards127",
        P127,
        kats::FP127_CERTS,
        cert16,
        select::verify_fp127
    );
    family!(
        "weier127",
        P127,
        kats::WEIER127_CERTS,
        cert16,
        select::verify_weier127
    );

    let got = format!("{HEADER}{}\n", out.join("\n"));
    let want = include_str!("fixtures/cert_corpus.txt");
    if std::env::var_os("CERT_CORPUS_RECORD").is_some() {
        std::fs::write(
            concat!(
                env!("CARGO_MANIFEST_DIR"),
                "/tests/fixtures/cert_corpus.txt"
            ),
            &got,
        )
        .unwrap();
    }
    assert_eq!(
        got, want,
        "verifier outcomes changed; see tests/fixtures/cert_corpus.txt"
    );
}
