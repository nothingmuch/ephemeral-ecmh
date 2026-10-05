//! Shared criterion configuration, synthetic fixtures and timing wrappers
//! with element throughput.

// Each benchmark binary uses a different subset of these helpers.
#![allow(dead_code)]

use criterion::measurement::WallTime;
use criterion::{BenchmarkGroup, Criterion, Throughput};
use ephemeral_ecmh::hash::Salted;
use std::hint::black_box;
use std::time::Duration;

/// Every suite's configuration, a fifth of criterion's defaults' time (3 s
/// warm-up, 5 s measurement, 100 samples, 100k resamples). Groups may set
/// their own; `criterion_group!` applies the command line over this, so
/// flags override it.
pub fn config() -> Criterion {
    Criterion::default()
        .warm_up_time(Duration::from_millis(500))
        .measurement_time(Duration::from_secs(1))
        .sample_size(30)
        .nresamples(10_000)
        // criterion still plots with gnuplot when it finds one
        .without_plots()
}

/// Synthetic 36-byte items, the size of an outpoint. `message` distinguishes
/// the empty-message fixtures from the four-byte tags used for RIBLT sets.
pub fn items(message: &[u8], n: usize) -> Vec<[u8; 36]> {
    (0..n as u32)
        .map(|i| {
            let mut x = [0u8; 36];
            x[..32].copy_from_slice(&Salted::new(b"bench-item", &[0; 32]).digest(message, i));
            x[32..].copy_from_slice(&i.to_le_bytes());
            x
        })
        .collect()
}

/// Apply `f` to every element per iteration; declare xs.len() elements
/// of throughput so the report can normalize the measured time.
pub fn each<T, R>(g: &mut BenchmarkGroup<WallTime>, name: &str, xs: &[T], f: impl Fn(&T) -> R) {
    g.throughput(Throughput::Elements(xs.len() as u64));
    g.bench_function(name, |bn| {
        bn.iter(|| {
            for x in xs {
                black_box(f(black_box(x)));
            }
        })
    });
}

/// Apply `f` to the whole slice per iteration; declare xs.len() elements
/// of throughput so the report can normalize the measured batch time.
pub fn whole<T, R>(g: &mut BenchmarkGroup<WallTime>, name: &str, xs: &[T], f: impl Fn(&[T]) -> R) {
    g.throughput(Throughput::Elements(xs.len() as u64));
    g.bench_function(name, |bn| bn.iter(|| f(black_box(xs))));
}
