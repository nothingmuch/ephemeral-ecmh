//! ECMH components and reference implementations: field arithmetic,
//! hash maps, point representations, and complete multiset digests.
//!
//!   nix develop -c cargo bench --bench compare
//!
//! Most iterations process N synthetic 36-byte items. Criterion records the
//! iteration time and element throughput; the report normalizes per element.
//! Ristretto255, libsecp256k1, and XOR-SHA256 provide fixed-group and hash
//! references with the specific APIs measured below.

mod common;
use common::each;

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use curve25519_dalek::ristretto::{CompressedRistretto, RistrettoPoint};
use ephemeral_ecmh::hash::Salted;
use secp256k1::PublicKey;
use secp256k1::ellswift::ElligatorSwift;
use sha2::{Digest, Sha512};

const N: usize = 1024;
/// Independent accumulators for the throughput-bound (RIBLT-like) adds.
const ACCS: usize = 8;

struct Setup {
    salt: Salted,
    items: Vec<[u8; 36]>,
}

fn setup() -> Setup {
    let seed = [0u8; 32];
    let salt = Salted::new(ephemeral_ecmh::ecmh::TAG_ITEM, &seed);
    let items = common::items(&[], N);
    Setup { salt, items }
}

impl Setup {
    fn refs(&self) -> Vec<&[u8]> {
        self.items.iter().map(|x| x.as_slice()).collect()
    }
}
/// Try-and-increment through libsecp256k1's compressed-point parser.
/// The fixed 0x02 prefix selects one y-coordinate sign, so the output
/// covers one representative of each sign pair.
fn secp_hash(h: &Salted, m: &[u8]) -> PublicKey {
    let mut buf = [2u8; 33];
    for ctr in 0.. {
        buf[1..].copy_from_slice(&h.digest(m, ctr));
        if let Ok(p) = PublicKey::from_byte_array_compressed(buf) {
            return p;
        }
    }
    unreachable!()
}

/// Decode two digests with libsecp256k1's BIP324 ElligatorSwift API.
fn secp_ellswift(h: &Salted, m: &[u8]) -> PublicKey {
    let mut b = [0u8; 64];
    b[..32].copy_from_slice(&h.digest(m, 0));
    b[32..].copy_from_slice(&h.digest(m, 1));
    PublicKey::from_ellswift(ElligatorSwift::from_byte_array(b))
}

/// Candidate x-coordinate checks for the reference groups. The study's
/// curve families use square tests, or trace tests with field inversion
/// for binary curves.
fn on_curve(c: &mut Criterion) {
    let s = setup();
    let mut g = c.benchmark_group("on_curve");
    let compressed: Vec<[u8; 33]> = s
        .items
        .iter()
        .map(|m| {
            let mut b = [2u8; 33];
            b[1..].copy_from_slice(&s.salt.digest(m, 0));
            b
        })
        .collect();
    each(
        &mut g,
        "secp256k1/x: parse compressed (sqrt)",
        &compressed,
        |b| PublicKey::from_byte_array_compressed(*b).is_ok(),
    );
    let rist: Vec<CompressedRistretto> = s
        .items
        .iter()
        .map(|m| RistrettoPoint::hash_from_bytes::<Sha512>(m).compress())
        .collect();
    each(
        &mut g,
        "ristretto255/decompress (valid, inv sqrt)",
        &rist,
        |b| b.decompress(),
    );
    g.finish();
}

/// The hash-to-curve alternatives; try-and-increment, the Group path, is
/// h2c/* in benches/group.rs.
fn hash_to_curve(c: &mut Criterion) {
    let s = setup();
    let refs = s.refs();
    let mut g = c.benchmark_group("hash_to_curve");
    each(&mut g, "sha256 (XOR baseline)", &refs, |m| {
        s.salt.digest(m, 0)
    });
    each(&mut g, "ristretto255/hash_from_bytes<Sha512>", &refs, |m| {
        RistrettoPoint::hash_from_bytes::<Sha512>(m)
    });
    each(&mut g, "secp256k1/try-and-increment (parse)", &refs, |m| {
        secp_hash(&s.salt, m)
    });
    each(&mut g, "secp256k1/ellswift decode", &refs, |m| {
        secp_ellswift(&s.salt, m)
    });
    g.finish();
}

/// Ristretto255's hash_from_bytes in its two steps: the SHA-512 digest,
/// then two Elligator maps and an add. The study's curves have no step
/// rows: their hashes are field operations, timed in `field`, around the
/// salted digest (`hash_to_curve`'s sha256).
fn h2c_parts(c: &mut Criterion) {
    let s = setup();
    let refs = s.refs();
    let mut g = c.benchmark_group("h2c_parts");
    each(&mut g, "digest/sha512 (ristretto input)", &refs, |m| {
        Sha512::digest(m)
    });
    let wide: Vec<[u8; 64]> = refs.iter().map(|m| Sha512::digest(m).into()).collect();
    each(
        &mut g,
        "ristretto255/from_uniform_bytes (2 maps + add)",
        &wide,
        RistrettoPoint::from_uniform_bytes,
    );
    g.finish();
}

/// Negations other than the addends' (group.neg in benches/group.rs).
fn negate(c: &mut Criterion) {
    let s = setup();
    let refs = s.refs();
    let hr: Vec<_> = refs
        .iter()
        .map(|m| RistrettoPoint::hash_from_bytes::<Sha512>(m))
        .collect();
    let hs: Vec<_> = refs.iter().map(|m| secp_hash(&s.salt, m)).collect();
    let mut g = c.benchmark_group("negate");
    each(&mut g, "ristretto255/-P", &hr, |p| -p);
    each(&mut g, "secp256k1/PublicKey::negate", &hs, |p| p.negate());
    g.finish();
}

fn add(c: &mut Criterion) {
    let s = setup();
    let refs = s.refs();
    let hr: Vec<_> = refs
        .iter()
        .map(|m| RistrettoPoint::hash_from_bytes::<Sha512>(m))
        .collect();
    let hx: Vec<[u8; 32]> = refs.iter().map(|m| s.salt.digest(m, 0)).collect();

    let mut g = c.benchmark_group("add");
    g.throughput(Throughput::Elements(N as u64));
    macro_rules! throughput {
        ($name:expr, $zero:expr, $xs:expr, $add:expr) => {
            g.bench_function(concat!($name, ", 8 accumulators"), |bn| {
                bn.iter(|| {
                    let mut acc = [$zero; ACCS];
                    for ch in $xs.chunks(ACCS) {
                        for (a, p) in acc.iter_mut().zip(ch) {
                            *a = $add(a, p);
                        }
                    }
                    acc
                })
            });
        };
    }
    g.bench_function("xor-sha256/xor 32B", |bn| {
        bn.iter(|| {
            hx.iter().fold([0u64; 4], |mut acc, d| {
                for (a, w) in acc.iter_mut().zip(d.as_chunks::<8>().0) {
                    *a ^= u64::from_le_bytes(*w);
                }
                acc
            })
        })
    });
    throughput!(
        "ristretto255/-=",
        RistrettoPoint::default(),
        hr,
        |a: &RistrettoPoint, p| a - p
    );
    g.finish();
}

/// Complete multiset digests, with batched or streaming accumulation as
/// named, compared with separately measured hash and sum components.
/// `riblt.encode` additionally maps each item to multiple cells and updates
/// their keys, counts, and point sums.
fn digest(c: &mut Criterion) {
    let s = setup();
    let refs = s.refs();
    let mut g = c.benchmark_group("digest");
    g.throughput(Throughput::Elements(N as u64));
    g.bench_function("xor-sha256", |bn| {
        bn.iter(|| {
            refs.iter().fold([0u8; 32], |mut acc, m| {
                for (a, d) in acc.iter_mut().zip(s.salt.digest(m, 0)) {
                    *a ^= d;
                }
                acc
            })
        })
    });
    g.bench_function("ristretto255/streaming", |bn| {
        bn.iter(|| {
            let acc = refs.iter().fold(RistrettoPoint::default(), |acc, m| {
                acc + RistrettoPoint::hash_from_bytes::<Sha512>(m)
            });
            acc.compress()
        })
    });
    g.finish();
}

criterion_group! {
    name = benches;
    config = common::config();
    targets = on_curve, hash_to_curve, h2c_parts, negate, add, digest
}
criterion_main!(benches);
