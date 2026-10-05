//! AGM point-counting and curve-selection costs, measured separately from
//! amortized ECMH operations. The ring-operation rows expose the components
//! of a count; the complete count measures their composition.
//!
//!   nix develop -c cargo bench --bench agm
//!
//! The nominal 127-bit schedule comprises 107 multiplications, 12 squarings,
//! 6 Frobenius applications and 63 inverse-Frobenius applications, followed
//! by an order consistency check using a 127-bit scalar multiplication.
//! `curvegen/find` measures the seed-dependent search; `curvegen/certify`
//! also constructs the rejection witnesses. Like `benches/curvegen.rs`,
//! both fix their sampling, which bench-run's profiles do not lengthen.

mod common;

use criterion::{Criterion, SamplingMode, criterion_group, criterion_main};
use ephemeral_ecmh::curvegen::agm::{Counter, Frobenius, Gf127, Zq, counter};
use ephemeral_ecmh::curvegen::prove::{certify_gf2_127, find_gf2_127};
use ephemeral_ecmh::curvegen::select::gf2_127_candidate;
use ephemeral_ecmh::field::gf2_127::to_u128;
use std::hint::black_box;
use std::time::Duration;

use common::families::kats;

fn count(c: &mut Criterion) {
    let k = &kats::GF2_127_CERTS[0];
    let curve = gf2_127_candidate(&k.seed, k.index).unwrap();
    let ctr = counter();
    let mut g = c.benchmark_group("agm");
    g.bench_function("order", |b| b.iter(|| ctr.order(black_box(&curve))));
    let b6 = to_u128(curve.big_b);
    g.bench_function("trace", |b| b.iter(|| ctr.trace(black_box(b6))));
    g.sample_size(10);
    g.bench_function("tables", |b| b.iter(Counter::<Gf127>::new));
    g.finish();
}

fn ring(c: &mut Criterion) {
    let f = Frobenius::<Gf127>::new();
    // dense 64-bit coefficients, like the lift's
    let x = f.sigma(&Zq::from_bits(0x1234_5678_9abc_def0_1357_9bdf_2468_ace0));
    let y = f.sigma(&x);
    let mut g = c.benchmark_group("zq");
    g.bench_function("mul", |b| b.iter(|| black_box(x) * black_box(y)));
    g.bench_function("square", |b| b.iter(|| black_box(x).square()));
    g.bench_function("sigma", |b| b.iter(|| f.sigma(black_box(&x))));
    g.bench_function("sigma_inv", |b| b.iter(|| f.sigma_inv(black_box(&x))));
    g.finish();
}

/// Seed-dependent curve search using the torsion sieve and AGM counts.
/// Each known-answer-test seed has its own identifier because the number
/// of candidates varies.
fn find(c: &mut Criterion) {
    let mut g = c.benchmark_group("curvegen/find");
    g.sampling_mode(SamplingMode::Flat)
        .sample_size(10)
        .warm_up_time(Duration::from_millis(500))
        .measurement_time(Duration::from_secs(1));
    for (i, k) in kats::GF2_127_CERTS.iter().enumerate() {
        g.bench_with_input(
            criterion::BenchmarkId::new(format!("gf2_127 agm+sieve/index {}", k.index), i),
            &k.seed,
            |b, seed| b.iter(|| find_gf2_127(black_box(seed))),
        );
    }
    g.finish();
    // the same counts, plus a witness per rejection: a hashed point, or a
    // point of order l multiplied down
    let mut g = c.benchmark_group("curvegen/certify");
    g.sampling_mode(SamplingMode::Flat)
        .sample_size(10)
        .warm_up_time(Duration::from_millis(500))
        .measurement_time(Duration::from_secs(1));
    for (i, k) in kats::GF2_127_CERTS.iter().enumerate() {
        g.bench_with_input(
            criterion::BenchmarkId::new(format!("gf2_127 agm+sieve/index {}", k.index), i),
            &k.seed,
            |b, seed| b.iter(|| certify_gf2_127(black_box(seed))),
        );
    }
    g.finish();
}

criterion_group! {
    name = benches;
    config = common::config();
    targets = count, ring, find
}
criterion_main!(benches);
