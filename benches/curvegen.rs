//! Seed-derived curve selection with PARI point counting. These setup
//! costs are measured separately from the amortized ECMH hot path.
//!
//!   nix develop -c cargo bench --features pari --bench curvegen
//!
//! Each group has an identifier per family and, where applicable, per
//! certificate fixture:
//!
//! - curvegen/count: a full PARI point count of one candidate, the accepted
//!   curve of the family's first certificate fixture. The 107-, 109-, 122-
//!   and 128-bit families have no Rust-side counter (their certificates
//!   come from Sage), so they have no count row.
//! - curvegen/verify_accept: acceptance checks on the selected candidate,
//!   including probable-prime screening, order witnesses, the Hasse bound,
//!   and the embedding-degree condition.
//! - curvegen/embedding: that condition alone on the accepted r, a
//!   baby-step giant-step search for q^k = 1 mod r with k <= 2^20, about
//!   2^11 products modulo r.
//! - curvegen/verify_full: selected-candidate acceptance checks and rejection
//!   checks for every preceding candidate in the deterministic sequence.
//!
//! Verification uses the first of each family's certificate fixtures from
//! `tests/common/kats*.rs`; verify_full names the number of rejections it
//! checks. Its cost over many seeds, single-shot, is in the curvegen_times
//! searches (`examples/curvegen_times.rs`), which also time whole searches
//! (find, prove), seconds to minutes depending on the seed.
//!
//! These are once-per-namespace costs, needed to a few percent, so every
//! group fixes its own sampling, criterion's minimum of ten samples over
//! a second, which bench-run's profiles do not lengthen.

mod common;

use criterion::{
    BenchmarkGroup, BenchmarkId, Criterion, SamplingMode, criterion_group, criterion_main,
    measurement::WallTime,
};
use ephemeral_ecmh::curvegen::criteria::Count;
use ephemeral_ecmh::curvegen::embedding_degree_ok;
use ephemeral_ecmh::curvegen::pari::{self, Pari};
use ephemeral_ecmh::curvegen::prove::{Binary127, Edwards127, Family, Weier127};
use ephemeral_ecmh::curvegen::select::{self, Certificate, Error};
use ephemeral_ecmh::curvegen::{select109, select122};
use std::hint::black_box;
use std::time::Duration;

use common::families::kats;
use common::families::kats109;
use common::families::kats122;

fn setup() {
    assert!(
        pari::seadata().is_some(),
        "no seadata: SEA would compute its own modular polynomials; set GP_DATA_DIR"
    );
}

/// A family's certificate fixtures: seed and certificate.
type Fixtures<const N: usize> = Vec<([u8; 32], Certificate<N>)>;

macro_rules! fixtures {
    ($kats:expr) => {
        $kats
            .iter()
            .map(|k| (k.seed, common::certificate(k.index, k.r, k.rejections)))
            .collect::<Fixtures<_>>()
    };
}

/// The accepted curve of a fixture, counted in full.
fn count_one<F: Family>(g: &mut BenchmarkGroup<WallTime>, name: &str, seed: &[u8; 32], j: u32)
where
    Pari: Count<F::Curve>,
{
    let curve = F::candidate(seed, j).unwrap();
    g.bench_function(name, |b| b.iter(|| Pari.order(black_box(&curve))));
}

/// A group sampled as the module states, whatever the profile.
fn group<'a>(c: &'a mut Criterion, name: &str) -> BenchmarkGroup<'a, WallTime> {
    let mut g = c.benchmark_group(name);
    g.sample_size(10)
        .warm_up_time(Duration::from_millis(500))
        .measurement_time(Duration::from_secs(1));
    g
}

fn count(c: &mut Criterion) {
    setup();
    let mut g = group(c, "curvegen/count");
    g.sampling_mode(SamplingMode::Flat);
    let k = &kats::GF2_127_CERTS[0];
    count_one::<Binary127>(&mut g, "gf2_127", &k.seed, k.index);
    let k = &kats::FP127_CERTS[0];
    count_one::<Edwards127>(&mut g, "edwards127", &k.seed, k.index);
    let k = &kats::WEIER127_CERTS[0];
    count_one::<Weier127>(&mut g, "weier127", &k.seed, k.index);
    g.finish();
}

type Accept<C> = fn(&[u8; 32], u32, u128) -> Result<C, Error>;
type Verify<C, const N: usize> = fn(&[u8; 32], &Certificate<N>) -> Result<C, Error>;

/// One family's verify_accept, embedding and verify_full rows, on its
/// first fixture, over a field of q elements.
fn family<C, const N: usize>(
    c: &mut Criterion,
    name: &str,
    q: u128,
    fixtures: &Fixtures<N>,
    accept: Accept<C>,
    verify: Verify<C, N>,
) {
    let mut g = group(c, "curvegen/verify_accept");
    for (i, (s, ct)) in fixtures.iter().enumerate().take(1) {
        g.bench_with_input(BenchmarkId::new(name, i), &(s, ct), |b, (s, ct)| {
            b.iter(|| accept(s, ct.index, ct.r).unwrap())
        });
    }
    g.finish();
    let mut g = group(c, "curvegen/embedding");
    for (i, (_, ct)) in fixtures.iter().enumerate().take(1) {
        assert!(
            embedding_degree_ok(q % ct.r, ct.r),
            "{name}: embedding degree"
        );
        g.bench_with_input(BenchmarkId::new(name, i), &ct.r, |b, &r| {
            b.iter(|| embedding_degree_ok(black_box(q) % r, r))
        });
    }
    g.finish();
    let mut g = group(c, "curvegen/verify_full");
    for (i, (s, ct)) in fixtures.iter().enumerate().take(1) {
        let id = format!("{name}/{} rejections", ct.rejections.len());
        g.bench_with_input(BenchmarkId::new(id, i), &(s, ct), |b, (s, ct)| {
            b.iter(|| verify(s, ct).unwrap())
        });
    }
    g.finish();
}

fn verify(c: &mut Criterion) {
    family(
        c,
        "gf2_127",
        1 << 127,
        &fixtures!(kats::GF2_127_CERTS),
        select::accept_gf2_127,
        select::verify_gf2_127,
    );
    family(
        c,
        "edwards127",
        (1 << 127) - 1,
        &fixtures!(kats::FP127_CERTS),
        select::accept_fp127,
        select::verify_fp127,
    );
    family(
        c,
        "weier127",
        (1 << 127) - 1,
        &fixtures!(kats::WEIER127_CERTS),
        select::accept_weier127,
        select::verify_weier127,
    );
    family(
        c,
        "gf2_109",
        1 << 109,
        &fixtures!(kats109::GF2_109_CERTS),
        select109::accept,
        select109::verify,
    );
    family(
        c,
        "gf2_122",
        1 << 122,
        &fixtures!(kats122::GF2_122_CERTS),
        select122::accept_dense,
        select122::verify_dense,
    );
    family(
        c,
        "gf2_122-gls",
        1 << 122,
        &fixtures!(kats122::GF2_122_GLS_CERTS),
        select122::accept_gls,
        select122::verify_gls,
    );
}

criterion_group! {
    name = benches;
    config = common::config();
    targets = count, verify
}
criterion_main!(benches);
