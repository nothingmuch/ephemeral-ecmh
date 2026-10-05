//! Torsion-sieve costs before curve point counting: a single torsion order
//! and the mean cost per candidate at the configured sieve bounds.
//!
//!   nix develop -c cargo bench --bench sieve
//!
//! Candidate fixtures and their curve orders come from `sage/sieve.sage`
//! and `tests/common/sieve_vectors.rs`. Operation-count estimates can be
//! compared with the measured costs.

mod common;

use criterion::measurement::WallTime;
use criterion::{
    BenchmarkGroup, BenchmarkId, Criterion, Throughput, criterion_group, criterion_main,
};
use ephemeral_ecmh::curvegen::select;
use ephemeral_ecmh::curvegen::sieve::{self, GF2_127_L_MAX};
use std::hint::black_box;

#[path = "../tests/common/sieve_vectors.rs"]
mod sv;

const LS: [u32; 5] = [3, 5, 7, 11, 13];

/// Measure a candidate whose curve order is not divisible by l. For odd
/// l, exclude divisibility of the twist order as well: twist points would
/// give roots of psi_l and require splitting the gcd. `two_q2` is
/// 2q + 2 modulo 2^128.
fn per_l<C>(
    g: &mut BenchmarkGroup<'_, WallTime>,
    family: &str,
    vectors: &[(u32, u128)],
    two_q2: u128,
    curve: impl Fn(u32) -> C,
    ls: &[u32],
    torsion: impl Fn(&C, u32) -> Option<[u8; 16]>,
) {
    for &l in ls {
        let l128 = l as u128;
        let &(j, _) = vectors
            .iter()
            .find(|&&(_, n)| n % l128 != 0 && (l % 2 == 0 || two_q2.wrapping_sub(n) % l128 != 0))
            .unwrap();
        let c = curve(j);
        g.bench_function(BenchmarkId::new(family, l), |b| {
            b.iter(|| torsion(black_box(&c), l))
        });
    }
}

fn torsion(c: &mut Criterion) {
    let mut g = c.benchmark_group("sieve/l");
    per_l(
        &mut g,
        "gf2_127",
        sv::GF2_127,
        2,
        |j| select::gf2_127_candidate(&sv::SEED, j).unwrap(),
        &LS,
        sieve::gf2_127_torsion,
    );
    g.finish();
}

/// Apply the full sieve to the fixture candidates at the configured
/// bounds; the report normalizes the iteration time per candidate.
fn per_candidate(c: &mut Criterion) {
    let mut g = c.benchmark_group("sieve/candidate");
    let gf2: Vec<_> = sv::GF2_127
        .iter()
        .map(|&(j, _)| select::gf2_127_candidate(&sv::SEED, j).unwrap())
        .collect();
    g.throughput(Throughput::Elements(gf2.len() as u64));
    g.bench_function(BenchmarkId::new("gf2_127", GF2_127_L_MAX), |b| {
        b.iter(|| {
            gf2.iter()
                .filter(|c| sieve::gf2_127(c, GF2_127_L_MAX).is_some())
                .count()
        })
    });
    g.finish();
}

criterion_group! {
    name = benches;
    config = common::config();
    targets = torsion, per_candidate
}
criterion_main!(benches);
