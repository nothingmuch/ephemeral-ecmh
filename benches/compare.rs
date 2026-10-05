//! ECMH components and reference implementations: field arithmetic,
//! hash maps, point representations, and complete multiset digests.
//!
//!   nix develop -c cargo bench --bench compare
//!
//! Curve fixtures have seed-derived dense parameters from
//! `tests/common/kats.rs`.
//! Most iterations process N synthetic 36-byte items. Criterion records the
//! iteration time and element throughput; the report normalizes per element.
//! Ristretto255, libsecp256k1, and XOR-SHA256 provide fixed-group and hash
//! references with the specific APIs measured below.

mod common;
use common::{each, packing, whole};

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use curve25519_dalek::ristretto::{CompressedRistretto, RistrettoPoint};
use ephemeral_ecmh::curve::h2c::{Elligator2, map1};
use ephemeral_ecmh::curve::{binary, binary127, edwards127};
use ephemeral_ecmh::ecmh::digest_batch;
use ephemeral_ecmh::field::fp127::Fp;
use ephemeral_ecmh::field::gf2_127::{self, Gf, MASK127, from_u128};
use ephemeral_ecmh::hash::{Salted, halves};
use secp256k1::PublicKey;
use secp256k1::ellswift::ElligatorSwift;
use sha2::{Digest, Sha512};
use std::hint::black_box;

const N: usize = 1024;
/// Independent accumulators for the throughput-bound (RIBLT-like) adds.
const ACCS: usize = 8;

struct Setup {
    seed: [u8; 32],
    bin: binary127::Curve,
    ed: edwards127::Curve,
    salt: Salted,
    items: Vec<[u8; 36]>,
    /// the first digest half of each item: a hash-to-curve candidate
    cands: Vec<u128>,
}

fn setup() -> Setup {
    let bin = common::families::binary127();
    let ed = common::families::edwards127();
    let seed = bin.seed;
    assert_eq!(seed, ed.seed);
    let salt = Salted::new(ephemeral_ecmh::ecmh::TAG_ITEM, &seed);
    let items = common::items(&[], N);
    let cands = items
        .iter()
        .map(|m| halves(&salt.digest(m, 0))[0])
        .collect();
    Setup {
        seed,
        bin: bin.group,
        ed: ed.group,
        salt,
        items,
        cands,
    }
}

impl Setup {
    fn refs(&self) -> Vec<&[u8]> {
        self.items.iter().map(|x| x.as_slice()).collect()
    }
    fn gf2_127_xs(&self) -> Vec<Gf> {
        self.cands.iter().map(|&c| from_u128(c | 1)).collect()
    }
    fn fp127_xs(&self) -> Vec<Fp> {
        self.cands.iter().map(|&c| Fp::new(c & MASK127)).collect()
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

// Fp has no MulAssign; `a = a * b` keeps the macro generic over both fields.
#[allow(clippy::assign_op_pattern)]
fn field(c: &mut Criterion) {
    let s = setup();
    let mut g = c.benchmark_group("field");
    let (gx, fx) = (s.gf2_127_xs(), s.fp127_xs());
    // Latency: one chain of dependent products. Throughput: ACCS chains.
    macro_rules! mul_chains {
        ($name:expr, $v:expr, $w:expr) => {
            g.throughput(Throughput::Elements(64));
            g.bench_function(concat!($name, "/mul latency (dependent chain)"), |bn| {
                bn.iter(|| (0..64).fold(black_box($v), |acc, _| acc * $w))
            });
            g.bench_function(concat!($name, "/mul throughput (8 chains)"), |bn| {
                bn.iter(|| {
                    // opaque as a whole: lanes known to be equal could share one chain
                    let mut acc = black_box([$v; ACCS]);
                    for _ in 0..64 / ACCS {
                        for a in acc.iter_mut() {
                            *a = *a * $w;
                        }
                    }
                    acc
                })
            });
        };
    }
    mul_chains!("gf2_127", gx[0], gx[1]);
    mul_chains!("fp127", fx[0], fx[1]);
    // squarings chain in a sqrt or an inversion's addition chain
    g.bench_function("gf2_127/square latency (dependent chain)", |bn| {
        bn.iter(|| (0..64).fold(black_box(gx[0]), |acc, _| acc.square()))
    });
    g.bench_function("gf2_127/square throughput (8 chains)", |bn| {
        bn.iter(|| {
            let mut acc = black_box([gx[0]; ACCS]);
            for _ in 0..64 / ACCS {
                for a in acc.iter_mut() {
                    *a = a.square();
                }
            }
            acc
        })
    });
    g.bench_function("fp127/square latency (dependent chain)", |bn| {
        bn.iter(|| (0..64).fold(black_box(fx[0]), |acc, _| acc.square()))
    });
    g.bench_function("fp127/square throughput (8 chains)", |bn| {
        bn.iter(|| {
            let mut acc = black_box([fx[0]; ACCS]);
            for _ in 0..64 / ACCS {
                for a in acc.iter_mut() {
                    *a = a.square();
                }
            }
            acc
        })
    });
    let fy: Vec<(Fp, Fp)> = fx
        .iter()
        .zip(fx.iter().rev())
        .map(|(&a, &b)| (a, b))
        .collect();
    each(&mut g, "fp127/add", &fy, |&(a, b)| a + b);
    each(&mut g, "fp127/neg", &fx, |&a| -a);
    each(&mut g, "gf2_127/invert", &gx, |v| v.invert());
    each(&mut g, "gf2_127/sqrt", &gx, |v| v.sqrt());
    each(&mut g, "gf2_127/halftrace", &gx, |v| v.halftrace());
    each(&mut g, "gf2_127/halftrace (byte tables)", &gx, |&v| {
        gf2_127::halftrace8(v)
    });
    each(&mut g, "gf2_127/halftrace (nibble tables)", &gx, |&v| {
        gf2_127::halftrace4(v)
    });
    each(&mut g, "gf2_127/trace", &gx, |v| v.trace());
    each(&mut g, "gf2_127/normalize (to_u128)", &gx, |&v| {
        gf2_127::to_u128(v)
    });
    whole(&mut g, "gf2_127/batch invert (product tree)", &gx, |v| {
        let mut w = v.to_vec();
        ephemeral_ecmh::field::batch::invert(&mut w);
        w
    });
    each(&mut g, "fp127/invert", &fx, |v| v.invert());
    each(&mut g, "fp127/sqrt (x^(2^125))", &fx, |v| v.sqrt());
    each(&mut g, "fp127/sqrt_ratio", &fy, |&(n, d)| {
        Fp::sqrt_ratio(n, d)
    });
    packing(&mut g, "fp127", &fx);
    each(&mut g, "fp127/pow_p34 (x^(2^125-1))", &fx, |v| v.pow_p34());
    whole(&mut g, "fp127/batch invert (product tree)", &fx, |v| {
        let mut w = v.to_vec();
        ephemeral_ecmh::field::batch::invert(&mut w);
        w
    });
    g.finish();
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
    each(&mut g, "gf2_127/pornin map x1", &refs, |m| {
        map1(&s.bin, &s.salt, m)
    });
    whole(&mut g, "gf2_127/pornin map x1, batched", &refs, |ms| {
        s.bin.hash_to_curve_map1_batch(&s.salt, ms)
    });

    // The map straight to a prepared addend: hash and prepare in one row, no
    // lift. One map, so the image is at most 2^(m-1) points of E[r] and not
    // uniform on it; see binary/map.rs.
    let u = binary::unscaled::Curve::new(s.bin);
    each(&mut g, "gf2_127-u/pornin map x1 to (u, v)", &refs, |m| {
        u.map_to_addend(halves(&s.salt.digest(m, 0))[0])
    });
    whole(
        &mut g,
        "gf2_127-u/pornin map x1 to (u, v), batched",
        &refs,
        |ms| {
            let cs: Vec<u128> = ms.iter().map(|m| halves(&s.salt.digest(m, 0))[0]).collect();
            u.map_to_addend_batch(&cs)
        },
    );
    each(
        &mut g,
        "gf2_127-lambda/pornin map x1 to (x, λ)",
        &refs,
        |m| s.bin.map_to_lambda(halves(&s.salt.digest(m, 0))[0]),
    );
    whole(
        &mut g,
        "gf2_127-lambda/pornin map x1 to (x, λ), batched",
        &refs,
        |ms| {
            let cs: Vec<u128> = ms.iter().map(|m| halves(&s.salt.digest(m, 0))[0]).collect();
            s.bin.map_to_lambda_batch(&cs)
        },
    );
    // the w codec's addend is the same (x, λ): its own rows, since the
    // report keys hash-to-addend recipes by family
    each(&mut g, "gf2_127-w/pornin map x1 to (x, λ)", &refs, |m| {
        s.bin.map_to_lambda(halves(&s.salt.digest(m, 0))[0])
    });
    whole(
        &mut g,
        "gf2_127-w/pornin map x1 to (x, λ), batched",
        &refs,
        |ms| {
            let cs: Vec<u128> = ms.iter().map(|m| halves(&s.salt.digest(m, 0))[0]).collect();
            s.bin.map_to_lambda_batch(&cs)
        },
    );
    // and try-and-increment to the same addend, as compare109.rs and
    // compare122.rs time it
    each(
        &mut g,
        "gf2_127-lambda/try-and-increment to (x, λ)",
        &refs,
        |m| s.bin.hash_to_lambda(&s.salt, m),
    );
    whole(
        &mut g,
        "gf2_127-lambda/try-and-increment to (x, λ), batched",
        &refs,
        |ms| s.bin.hash_to_lambda_batch(&s.salt, ms),
    );
    each(&mut g, "gf2_127/pornin map x2", &refs, |m| {
        s.bin.hash_to_curve_map2(&s.salt, m)
    });
    each(
        &mut g,
        "edwards127/edwards-native try-and-increment",
        &refs,
        |m| s.ed.hash_to_edwards(&s.salt, m),
    );
    let ell = Elligator2::new(s.ed).unwrap();
    each(&mut g, "edwards127/elligator2 x1", &refs, |m| {
        map1(&ell, &s.salt, m)
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
    let hb = s.bin.hash_to_curve_batch(&s.salt, &refs);
    let hm: Vec<_> = refs
        .iter()
        .map(|m| s.ed.hash_to_curve(&s.salt, m))
        .collect();
    let em: Vec<_> = hm.iter().map(|p| s.ed.from_affine(p)).collect();
    let hr: Vec<_> = refs
        .iter()
        .map(|m| RistrettoPoint::hash_from_bytes::<Sha512>(m))
        .collect();
    let hs: Vec<_> = refs.iter().map(|m| secp_hash(&s.salt, m)).collect();
    let mut g = c.benchmark_group("negate");
    each(&mut g, "gf2_127/affine (y += x)", &hb, |p| p.neg());
    each(&mut g, "edwards127/extended (-X, -T)", &em, |p| p.neg());
    each(&mut g, "ristretto255/-P", &hr, |p| -p);
    each(&mut g, "secp256k1/PublicKey::negate", &hs, |p| p.negate());
    g.finish();
}

fn add(c: &mut Criterion) {
    let s = setup();
    let refs = s.refs();
    let hb = s.bin.hash_to_curve_batch(&s.salt, &refs);
    let hm: Vec<_> = refs
        .iter()
        .map(|m| s.ed.hash_to_curve(&s.salt, m))
        .collect();
    let hr: Vec<_> = refs
        .iter()
        .map(|m| RistrettoPoint::hash_from_bytes::<Sha512>(m))
        .collect();
    let hx: Vec<[u64; 4]> = refs
        .iter()
        .map(|m| {
            let d = s.salt.digest(m, 0);
            core::array::from_fn(|i| u64::from_le_bytes(d.as_chunks::<8>().0[i]))
        })
        .collect();

    let mut g = c.benchmark_group("add");
    g.throughput(Throughput::Elements(N as u64));
    // The RIBLT forms (extended, Cached, RCB) are in benches/group.rs; these
    // are the other representations, throughput-bound: ACCS accumulators
    // round robin, like updating distinct RIBLT cells.
    macro_rules! streaming {
        ($name:expr, $zero:expr, $xs:expr, $add:expr) => {
            g.bench_function(concat!($name, ", 1 accumulator"), |bn| {
                bn.iter(|| $xs.iter().fold($zero, |acc, p| $add(&acc, p)))
            });
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
    streaming!(
        "xor-sha256/xor 32B",
        [0u64; 4],
        hx,
        |a: &[u64; 4], p: &[u64; 4]| { core::array::from_fn::<u64, 4, _>(|i| a[i] ^ p[i]) }
    );
    g.bench_function("gf2_127/batch affine (tree sum)", |bn| {
        bn.iter(|| binary127::sum_batch(black_box(&hb)))
    });
    throughput!("gf2_127/extended += affine", s.bin.neutral(), hb, |a, p| s
        .bin
        .add_affine(a, p));
    g.bench_function("edwards127/batch montgomery affine (tree sum)", |bn| {
        bn.iter(|| s.ed.sum_batch(black_box(&hm)))
    });
    throughput!(
        "edwards127/+= montgomery affine",
        edwards127::Point::IDENTITY,
        hm,
        |a, p| s.ed.add_ext(a, &s.ed.from_affine(p))
    );
    streaming!(
        "ristretto255/+=",
        RistrettoPoint::default(),
        hr,
        |a: &RistrettoPoint, p| a + p
    );
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
    g.bench_function("gf2_127/batch", |bn| {
        bn.iter(|| digest_batch(s.bin, &s.seed, &refs))
    });
    g.bench_function("edwards127/batch", |bn| {
        bn.iter(|| digest_batch(s.ed, &s.seed, &refs))
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
    targets = field, on_curve, hash_to_curve, h2c_parts, negate, add, digest
}
criterion_main!(benches);
