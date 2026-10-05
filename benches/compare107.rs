//! Field, hash, and representation costs for the 14-byte curve encodings
//! over F_p, p = 2^107 - 1, using the operation groups of `compare.rs`.
//!
//!   nix develop -c cargo bench --bench compare107
//!
//! The suite uses the first certificate fixture for each family and N
//! synthetic 36-byte items. Identifiers include `107` to distinguish these
//! measurements from the 127-bit instances. The report normalizes iteration
//! times per element.

mod common;
use common::{each, packing, whole};

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use ephemeral_ecmh::curve::h2c::{Elligator2, map1};
use ephemeral_ecmh::curve::weier::Sswu;
use ephemeral_ecmh::curve::{edwards107, weier107};
use ephemeral_ecmh::ecmh::digest_batch;
use ephemeral_ecmh::field::fp107::{Fp, MASK107};
use ephemeral_ecmh::hash::{Salted, halves};
use std::hint::black_box;

const N: usize = 1024;
/// Independent accumulators for the throughput-bound (RIBLT-like) adds.
const ACCS: usize = 8;

struct Setup {
    seed: [u8; 32],
    ed: edwards107::Curve,
    wei: weier107::OddCurve,
    salt: Salted,
    items: Vec<[u8; 36]>,
    /// the first digest half of each item: a hash-to-curve candidate
    cands: Vec<u128>,
}

fn setup() -> Setup {
    let ed = common::families::edwards107();
    let wei = common::families::weier107();
    let seed = ed.seed;
    assert_eq!(seed, wei.seed);
    let salt = Salted::new(ephemeral_ecmh::ecmh::TAG_ITEM, &seed);
    let items = common::items(&[], N);
    let cands = items
        .iter()
        .map(|m| halves(&salt.digest(m, 0))[0])
        .collect();
    Setup {
        seed,
        ed: ed.group,
        wei: wei.group,
        salt,
        items,
        cands,
    }
}

impl Setup {
    fn refs(&self) -> Vec<&[u8]> {
        self.items.iter().map(|x| x.as_slice()).collect()
    }
    fn fp107_xs(&self) -> Vec<Fp> {
        self.cands.iter().map(|&c| Fp::new(c & MASK107)).collect()
    }
}

fn field(c: &mut Criterion) {
    let s = setup();
    let mut g = c.benchmark_group("field");
    let fx = s.fp107_xs();
    // Latency: one chain of dependent products. Throughput: ACCS chains.
    macro_rules! mul_chains {
        ($name:expr, $v:expr, $w:expr, $mul:expr) => {
            g.throughput(Throughput::Elements(64));
            g.bench_function(concat!($name, "/mul latency (dependent chain)"), |bn| {
                bn.iter(|| (0..64).fold(black_box($v), |acc, _| $mul(acc, $w)))
            });
            g.bench_function(concat!($name, "/mul throughput (8 chains)"), |bn| {
                bn.iter(|| {
                    // opaque as a whole: lanes known to be equal could share one chain
                    let mut acc = black_box([$v; ACCS]);
                    for _ in 0..64 / ACCS {
                        for a in acc.iter_mut() {
                            *a = $mul(*a, $w);
                        }
                    }
                    acc
                })
            });
        };
    }
    mul_chains!("fp107", fx[0], fx[1], |a: Fp, b| a * b);
    mul_chains!("fp107 schoolbook", fx[0], fx[1], Fp::mul_schoolbook);
    // squarings chain in a sqrt or an inversion's addition chain
    g.bench_function("fp107/square latency (dependent chain)", |bn| {
        bn.iter(|| (0..64).fold(black_box(fx[0]), |acc, _| acc.square()))
    });
    g.bench_function("fp107/square throughput (8 chains)", |bn| {
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
    each(&mut g, "fp107/add", &fy, |&(a, b)| a + b);
    each(&mut g, "fp107/neg", &fx, |&a| -a);
    each(&mut g, "fp107/invert", &fx, |v| v.invert());
    each(&mut g, "fp107/sqrt (x^(2^105))", &fx, |v| v.sqrt());
    each(&mut g, "fp107/sqrt_ratio", &fy, |&(n, d)| {
        Fp::sqrt_ratio(n, d)
    });
    packing(&mut g, "fp107", &fx);
    each(&mut g, "fp107/pow_p34 (x^(2^105-1))", &fx, |v| v.pow_p34());
    whole(&mut g, "fp107/batch invert (product tree)", &fx, |v| {
        let mut w = v.to_vec();
        ephemeral_ecmh::field::batch::invert(&mut w);
        w
    });
    g.finish();
}

/// The alternatives to try-and-increment (h2c/* in benches/group.rs).
fn hash_to_curve(c: &mut Criterion) {
    let s = setup();
    let refs = s.refs();
    let mut g = c.benchmark_group("hash_to_curve");
    each(
        &mut g,
        "edwards107/edwards-native try-and-increment",
        &refs,
        |m| s.ed.hash_to_edwards(&s.salt, m),
    );
    let ell = Elligator2::new(s.ed).unwrap();
    each(&mut g, "edwards107/elligator2 x1", &refs, |m| {
        map1(&ell, &s.salt, m)
    });
    let sswu = Sswu::search(s.wei.curve(), 256).expect("fixture admits SSWU setup");
    each(&mut g, "weier107/sswu x1", &refs, |m| {
        map1(&sswu, &s.salt, m)
    });
    g.finish();
}

/// Negations other than the addends' (group.neg in benches/group.rs).
fn negate(c: &mut Criterion) {
    let s = setup();
    let refs = s.refs();
    let hm: Vec<_> = refs
        .iter()
        .map(|m| s.ed.hash_to_curve(&s.salt, m))
        .collect();
    let em: Vec<_> = hm.iter().map(|p| s.ed.from_affine(p)).collect();
    let hw: Vec<_> = refs
        .iter()
        .map(|m| s.wei.curve().hash_to_curve(&s.salt, m))
        .collect();
    let ew: Vec<_> = hw
        .iter()
        .map(|p| s.wei.curve().add_affine(&weier107::Point::IDENTITY, p))
        .collect();
    let mut g = c.benchmark_group("negate");
    each(&mut g, "edwards107/extended (-X, -T)", &em, |p| p.neg());
    each(&mut g, "weier107/projective (-Y)", &ew, |p| p.neg());
    g.finish();
}

fn add(c: &mut Criterion) {
    let s = setup();
    let refs = s.refs();
    let hm: Vec<_> = refs
        .iter()
        .map(|m| s.ed.hash_to_curve(&s.salt, m))
        .collect();
    let hw: Vec<_> = refs
        .iter()
        .map(|m| s.wei.curve().hash_to_curve(&s.salt, m))
        .collect();

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
    g.bench_function("edwards107/batch montgomery affine (tree sum)", |bn| {
        bn.iter(|| s.ed.sum_batch(black_box(&hm)))
    });
    throughput!(
        "edwards107/+= montgomery affine",
        edwards107::Point::IDENTITY,
        hm,
        |a, p| s.ed.add_ext(a, &s.ed.from_affine(p))
    );
    g.bench_function("weier107/batch affine (tree sum)", |bn| {
        bn.iter(|| s.wei.curve().sum_batch(black_box(&hw)))
    });
    g.finish();
}

fn digest(c: &mut Criterion) {
    let s = setup();
    let refs = s.refs();
    let mut g = c.benchmark_group("digest");
    g.throughput(Throughput::Elements(N as u64));
    g.bench_function("edwards107/batch", |bn| {
        bn.iter(|| digest_batch(s.ed, &s.seed, &refs))
    });
    g.bench_function("weier107/batch", |bn| {
        bn.iter(|| digest_batch(s.wei, &s.seed, &refs))
    });
    g.finish();
}

criterion_group! {
    name = benches;
    config = common::config();
    targets = field, hash_to_curve, negate, add, digest
}
criterion_main!(benches);
