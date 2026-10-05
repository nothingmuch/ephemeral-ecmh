//! Field, hash, and representation costs for the 14-byte curve encodings
//! over GF(2^109), using the operation groups of `compare.rs`.
//!
//!   nix develop -c cargo bench --bench compare109
//!
//! Curve fixtures come from `tests/common/kats109.rs`. The direct hash to
//! lambda-affine (x, lambda) avoids the preparation inversion. This is
//! measured separately from `group.rs`'s `binary-lambda.109` hash rows, which
//! return ordinary affine points and share the `binary.109` hash.

mod common;
use common::{each, whole};

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use ephemeral_ecmh::curve::binary;
use ephemeral_ecmh::curve::binary109::{self, candidate};
use ephemeral_ecmh::curve::h2c::map1;
use ephemeral_ecmh::ecmh::digest_batch;
use ephemeral_ecmh::field::gf2_109::Gf;
use ephemeral_ecmh::hash::{Salted, halves};
use std::hint::black_box;

const N: usize = 1024;
/// Independent accumulators for the throughput-bound (RIBLT-like) adds.
const ACCS: usize = 8;

struct Setup {
    seed: [u8; 32],
    bin: binary109::Curve,
    salt: Salted,
    items: Vec<[u8; 36]>,
    /// the first digest half of each item: a hash-to-curve candidate
    cands: Vec<u128>,
}

fn setup() -> Setup {
    let k = common::families::binary109();
    let salt = Salted::new(ephemeral_ecmh::ecmh::TAG_ITEM, &k.seed);
    let items = common::items(&[], N);
    let cands = items
        .iter()
        .map(|m| halves(&salt.digest(m, 0))[0])
        .collect();
    Setup {
        seed: k.seed,
        bin: k.group,
        salt,
        items,
        cands,
    }
}

impl Setup {
    fn refs(&self) -> Vec<&[u8]> {
        self.items.iter().map(|x| x.as_slice()).collect()
    }
    fn xs(&self) -> Vec<Gf> {
        self.cands.iter().map(|&c| candidate(c).0).collect()
    }
}

fn field(c: &mut Criterion) {
    let s = setup();
    let mut g = c.benchmark_group("field");
    let gx = s.xs();
    // Latency: one chain of dependent products. Throughput: ACCS chains.
    let (v, w) = (gx[0], gx[1]);
    g.throughput(Throughput::Elements(64));
    g.bench_function("gf2_109/mul latency (dependent chain)", |bn| {
        bn.iter(|| (0..64).fold(black_box(v), |acc, _| acc * w))
    });
    g.bench_function("gf2_109/mul throughput (8 chains)", |bn| {
        bn.iter(|| {
            // opaque as a whole: lanes known to be equal could share one chain
            let mut acc = black_box([v; ACCS]);
            for _ in 0..64 / ACCS {
                for a in acc.iter_mut() {
                    *a *= w;
                }
            }
            acc
        })
    });
    // squarings chain in a sqrt or an inversion's addition chain
    g.bench_function("gf2_109/square latency (dependent chain)", |bn| {
        bn.iter(|| (0..64).fold(black_box(v), |acc, _| acc.square()))
    });
    g.bench_function("gf2_109/square throughput (8 chains)", |bn| {
        bn.iter(|| {
            let mut acc = black_box([v; ACCS]);
            for _ in 0..64 / ACCS {
                for a in acc.iter_mut() {
                    *a = a.square();
                }
            }
            acc
        })
    });
    each(&mut g, "gf2_109/invert", &gx, |v| v.invert());
    each(&mut g, "gf2_109/sqrt", &gx, |v| v.sqrt());
    each(&mut g, "gf2_109/halftrace", &gx, |v| v.halftrace());
    each(&mut g, "gf2_109/halftrace (byte tables)", &gx, |&v| {
        ephemeral_ecmh::field::gf2_109::halftrace8(v)
    });
    each(&mut g, "gf2_109/trace", &gx, |v| v.trace());
    each(&mut g, "gf2_109/normalize (to_u128)", &gx, |&v| {
        ephemeral_ecmh::field::gf2_109::to_u128(v)
    });
    whole(&mut g, "gf2_109/batch invert (product tree)", &gx, |v| {
        let mut w = v.to_vec();
        ephemeral_ecmh::field::batch::invert(&mut w);
        w
    });
    g.finish();
}

/// The λ family's own hash: (x, λ) out, no lift; and Pornin's map, as
/// `compare.rs` times it over GF(2^127).
fn h2c_parts(c: &mut Criterion) {
    let s = setup();
    let refs = s.refs();
    let mut g = c.benchmark_group("hash_to_curve");
    each(&mut g, "gf2_109/pornin map x1", &refs, |m| {
        map1(&s.bin, &s.salt, m)
    });
    whole(&mut g, "gf2_109/pornin map x1, batched", &refs, |ms| {
        s.bin.hash_to_curve_map1_batch(&s.salt, ms)
    });

    // The map straight to a prepared addend: hash and prepare in one row, no
    // lift. One map, so the image is at most 2^(m-1) points of E[r] and not
    // uniform on it; see binary/map.rs.
    let u = binary::unscaled::Curve::new(s.bin);
    each(&mut g, "gf2_109-u/pornin map x1 to (u, v)", &refs, |m| {
        u.map_to_addend(halves(&s.salt.digest(m, 0))[0])
    });
    whole(
        &mut g,
        "gf2_109-u/pornin map x1 to (u, v), batched",
        &refs,
        |ms| {
            let cs: Vec<u128> = ms.iter().map(|m| halves(&s.salt.digest(m, 0))[0]).collect();
            u.map_to_addend_batch(&cs)
        },
    );
    each(
        &mut g,
        "gf2_109-lambda/pornin map x1 to (x, λ)",
        &refs,
        |m| s.bin.map_to_lambda(halves(&s.salt.digest(m, 0))[0]),
    );
    whole(
        &mut g,
        "gf2_109-lambda/pornin map x1 to (x, λ), batched",
        &refs,
        |ms| {
            let cs: Vec<u128> = ms.iter().map(|m| halves(&s.salt.digest(m, 0))[0]).collect();
            s.bin.map_to_lambda_batch(&cs)
        },
    );
    // the w codec's addend is the same (x, λ): its own rows, since the
    // report keys hash-to-addend recipes by family
    each(&mut g, "gf2_109-w/pornin map x1 to (x, λ)", &refs, |m| {
        s.bin.map_to_lambda(halves(&s.salt.digest(m, 0))[0])
    });
    whole(
        &mut g,
        "gf2_109-w/pornin map x1 to (x, λ), batched",
        &refs,
        |ms| {
            let cs: Vec<u128> = ms.iter().map(|m| halves(&s.salt.digest(m, 0))[0]).collect();
            s.bin.map_to_lambda_batch(&cs)
        },
    );
    each(
        &mut g,
        "gf2_109-lambda/try-and-increment to (x, λ)",
        &refs,
        |m| s.bin.hash_to_lambda(&s.salt, m),
    );
    whole(
        &mut g,
        "gf2_109-lambda/try-and-increment to (x, λ), batched",
        &refs,
        |ms| s.bin.hash_to_lambda_batch(&s.salt, ms),
    );
    g.finish();
}

/// The affine negation; the addends' is group.neg in benches/group.rs.
fn negate(c: &mut Criterion) {
    let s = setup();
    let hb = s.bin.hash_to_curve_batch(&s.salt, &s.refs());
    let mut g = c.benchmark_group("negate");
    each(&mut g, "gf2_109/affine (y += x)", &hb, |p| p.neg());
    g.finish();
}

fn add(c: &mut Criterion) {
    let s = setup();
    let hb = s.bin.hash_to_curve_batch(&s.salt, &s.refs());
    let eb: Vec<_> = hb.iter().map(|p| s.bin.from_affine(p)).collect();

    let mut g = c.benchmark_group("add");
    g.throughput(Throughput::Elements(N as u64));
    // The RIBLT forms (extended, Cached, RCB) are in benches/group.rs; these
    // are the other representations, throughput-bound: ACCS accumulators
    // round robin, like updating distinct RIBLT cells.
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
    g.bench_function("gf2_109/batch affine (tree sum)", |bn| {
        bn.iter(|| binary109::sum_batch(black_box(&hb)))
    });
    throughput!("gf2_109/extended += affine", s.bin.neutral(), hb, |a, p| s
        .bin
        .add_affine(a, p));
    whole(&mut g, "gf2_109/normalize, batched", &eb, |ps| {
        s.bin.normalize_batch(ps)
    });
    g.finish();
}

fn digest(c: &mut Criterion) {
    let s = setup();
    let refs = s.refs();
    let mut g = c.benchmark_group("digest");
    g.throughput(Throughput::Elements(N as u64));
    g.bench_function("gf2_109/batch", |bn| {
        bn.iter(|| digest_batch(s.bin, &s.seed, &refs))
    });
    g.finish();
}

criterion_group! {
    name = benches;
    config = common::config();
    targets = field, h2c_parts, negate, add, digest
}
criterion_main!(benches);
