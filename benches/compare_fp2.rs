//! Prime and quadratic-extension field arithmetic for ECMH candidate
//! selection: fp61x2 (p = 2^61 - 1), fp64x2 (p = 2^64 - 59), and goldilocks2
//! (p = 2^64 - 2^32 + 1), in the shared `field` benchmark group. The suite
//! also measures fp128 (p = 2^128 - 275), used by the a = -1 Edwards
//! curves, as a two-word prime field without an extension. The base fields
//! appear only inside their extensions' operations.
//!
//! Curves over these fields, their complete hashes and codecs, are measured
//! in `group.rs`. `hash_to_curve` here times one map per item: Elligator 2
//! on each fixture with a Montgomery model, the a = 1 Edwards curve over
//! fp61x2 and the a = -1 twisted curves over all four fields, and SSWU on
//! the prime-order Weierstrass curve over fp61x2.

mod common;
use common::{each, packing, whole};

use criterion::{Criterion, Throughput, criterion_group, criterion_main};
use ephemeral_ecmh::curve::h2c::{Elligator2, map1};
use ephemeral_ecmh::curve::twisted::MontgomeryModel;
use ephemeral_ecmh::curve::weier::Sswu;
use ephemeral_ecmh::field::{OddField, batch, fp61x2, fp64x2, fp128, goldilocks2};
use ephemeral_ecmh::hash::{Salted, halves};
use p3_field::PrimeCharacteristicRing;
use std::hint::black_box;

const N: usize = 1024;
/// Independent chains for the throughput-bound products.
const ACCS: usize = 8;

/// N pseudo-random digest halves, split into two words each.
fn words() -> Vec<(u64, u64)> {
    let h = Salted::new(b"bench-field", &[0; 32]);
    (0..N as u32)
        .map(|i| {
            let v = halves(&h.digest(&[], i))[0];
            (v as u64, (v >> 64) as u64)
        })
        .collect()
}

// the chains spell `acc * w` for every field; only Plonky3's have `*=`
#[allow(clippy::assign_op_pattern)]
fn field(c: &mut Criterion) {
    let w = words();
    let qx: Vec<fp61x2::Fq> = w
        .iter()
        .map(|&(a, b)| fp61x2::Fq::new(fp61x2::Fp::new(a), fp61x2::Fp::new(b)))
        .collect();
    let rx: Vec<fp64x2::Fq> = w
        .iter()
        .map(|&(a, b)| fp64x2::Fq::new(fp64x2::Fp::new(a), fp64x2::Fp::new(b)))
        .collect();
    let gx: Vec<goldilocks2::Fq> = w
        .iter()
        .map(|&(a, b)| goldilocks2::new(goldilocks2::Fp::new(a), goldilocks2::Fp::new(b)))
        .collect();
    let px: Vec<fp128::Fp> = w
        .iter()
        .map(|&(a, b)| fp128::Fp::new(a as u128 | (b as u128) << 64))
        .collect();
    let mut g = c.benchmark_group("field");
    // Latency: one chain of dependent products. Throughput: ACCS chains.
    macro_rules! chains {
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
            // squarings chain in a sqrt or an inversion's addition chain
            g.bench_function(concat!($name, "/square latency (dependent chain)"), |bn| {
                bn.iter(|| (0..64).fold(black_box($v), |acc, _| acc.square()))
            });
            g.bench_function(concat!($name, "/square throughput (8 chains)"), |bn| {
                bn.iter(|| {
                    let mut acc = black_box([$v; ACCS]);
                    for _ in 0..64 / ACCS {
                        for a in acc.iter_mut() {
                            *a = a.square();
                        }
                    }
                    acc
                })
            });
        };
    }
    chains!("fp61x2", qx[0], qx[1]);
    chains!("fp64x2", rx[0], rx[1]);
    chains!("goldilocks2", gx[0], gx[1]);
    chains!("fp128", px[0], px[1]);
    // goldilocks2's inversion and square root are free functions, Plonky3's
    // types being foreign
    macro_rules! ops {
        ($name:expr, $xs:expr) => {
            ops!($name, $xs, |v| v.invert(), |v| v.sqrt());
        };
        ($name:expr, $xs:expr, $invert:expr, $sqrt:expr) => {
            let ys: Vec<_> = $xs
                .iter()
                .zip($xs.iter().rev())
                .map(|(&a, &b)| (a, b))
                .collect();
            each(&mut g, concat!($name, "/add"), &ys, |&(a, b)| a + b);
            each(&mut g, concat!($name, "/sub"), &ys, |&(a, b)| a - b);
            each(&mut g, concat!($name, "/neg"), &$xs, |&a| -a);
            each(&mut g, concat!($name, "/invert"), &$xs, $invert);
            each(&mut g, concat!($name, "/sqrt"), &$xs, $sqrt);
            each(&mut g, concat!($name, "/sqrt_ratio"), &ys, |&(n, d)| {
                OddField::sqrt_ratio(n, d)
            });
            packing(&mut g, $name, &$xs);
        };
    }
    ops!("fp61x2", qx);
    ops!("fp64x2", rx);
    ops!("fp128", px);
    ops!("goldilocks2", gx, |&v| goldilocks2::invert(v), |&v| {
        goldilocks2::sqrt(v)
    });
    whole(&mut g, "fp61x2/batch invert (product tree)", &qx, |v| {
        let mut w = v.to_vec();
        batch::invert(&mut w);
        w
    });
    whole(&mut g, "fp64x2/batch invert (product tree)", &rx, |v| {
        let mut w = v.to_vec();
        batch::invert(&mut w);
        w
    });
    whole(
        &mut g,
        "goldilocks2/batch invert (product tree)",
        &gx,
        |v| {
            let mut w = v.to_vec();
            batch::invert(&mut w);
            w
        },
    );
    g.finish();
}

/// Elligator 2 on each fixture's Montgomery model, one map per item: not
/// uniform on its own, and not the hash `group.rs` times.
fn hash_to_curve(c: &mut Criterion) {
    use common::families;
    let ed = families::edwards61x2();
    let salt = Salted::new(ephemeral_ecmh::ecmh::TAG_ITEM, &ed.seed);
    let items = common::items(&[], N);
    let refs: Vec<&[u8]> = items.iter().map(|x| x.as_slice()).collect();
    let mut g = c.benchmark_group("hash_to_curve");
    let ell = Elligator2::new(ed.group).expect("fixture has a2 != 0");
    each(&mut g, "edwards61x2/elligator2 x1", &refs, |m| {
        map1(&ell, &salt, m)
    });
    macro_rules! twisted {
        ($name:literal, $fixture:expr) => {
            let ell = Elligator2::new(MontgomeryModel::new($fixture.group)).unwrap();
            each(&mut g, concat!($name, "/elligator2 x1"), &refs, |m| {
                map1(&ell, &salt, m)
            });
        };
    }
    twisted!("twisted128", families::twisted128());
    twisted!("twisted61x2", families::twisted61x2());
    twisted!("twisted64x2", families::twisted64x2());
    twisted!("twisted-goldilocks2", families::twisted_goldilocks2());
    let i = fp61x2::Fq::new(fp61x2::Fp::ZERO, fp61x2::Fp::ONE);
    let sswu = Sswu::search_from(families::weier61x2().group.curve(), i, 64)
        .expect("fixture admits SSWU setup");
    each(&mut g, "weier61x2/sswu x1", &refs, |m| {
        map1(&sswu, &salt, m)
    });
    g.finish();
}

criterion_group! {
    name = benches;
    config = common::config();
    targets = field, hash_to_curve
}
criterion_main!(benches);
