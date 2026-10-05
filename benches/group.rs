//! Group capabilities relevant to ECMH and RIBLT costs, measured through
//! one generic harness for each included curve and representation.
//!
//!   nix develop -c cargo bench --bench group
//!
//! Identifiers are `<layer>.<op>/<family>.<bits>/<parameters>`, with `h2c`
//! as the hash layer:
//!
//! - h2c: hash one item (mode=indep), or chunks of n (mode=batch,n=...).
//!   Preparation, encoding, and decoding use this sweep.
//! - h2c.id: hash N 32-byte IDs to addends as one batch, under each
//!   construction benches/common/h2c.rs lists (`h2c=<token>`) and each
//!   per-salt hash of an ID (`proj=`): salted SHA-256 (`none`), or its
//!   projection in F_(2^130 - 5), F_(2^127 - 1) or GF(2^127) (`fp130`,
//!   `fp127`, `gf2_127`).
//! - group.prepare: convert a hash output to the addend reused across cells:
//!   extended on binary curves, lambda-affine on binary-lambda, Cached on
//!   Edwards, and affine on Weierstrass. The binary-w and binary-u hashes
//!   already return their addend representation, so preparation is a copy.
//! - group.add: update one dependent accumulator (mode=latency), or ACCS
//!   independent accumulators in round-robin order (mode=throughput).
//! - group.sub: the same accumulation with addend negation on every update,
//!   mode=throughput only: subtraction is negation and addition, both
//!   measured, so the row is decoding's check on their composition.
//! - group.neg: negate a prepared addend.
//! - group.is_identity: the empty-cell test, applied to running-sum
//!   accumulator fixtures.
//! - group.equals: the purity check, comparing a cell's sum with an addend:
//!   mode=mismatch on running sums against the next addend (impure cells,
//!   peeling's common case, where a comparison may stop at the first
//!   coordinate), mode=match on sums equal to the addend at a non-unit
//!   scale (pure cells).
//! - group.encode: produce accumulators' canonical point bytes for
//!   transmission of cells, through the affine form unless the group
//!   encodes accumulators directly. Batch measurements share normalization
//!   work.
//! - group.decode: recover the group's affine representation from point bytes,
//!   independently or in batches; the w codec returns lambda-affine points.
//!
//! Coverage is planned from a first pass, as for the RIBLT suite.
//! `GROUP_PASS=core` times every family's core rows: batches at N only and
//! the mismatch equality mode. From those and the comparison maps,
//! `bench-report --group-plan` names the families whose modeled insertion
//! is within the RIBLT buffer ratio of the cheapest (`lead`), and, if
//! none of them is over a field of odd characteristic, the cheapest that
//! is (`contrast`). `GROUP_PASS=wide`, with `GROUP_PLAN` naming that TSV of
//! `<family>\t<role>` lines, times the planned families' remaining rows:
//! the batch sweep at n = 8 and 64, subtraction, the match equality mode
//! and h2c.id. Unset, `GROUP_PASS` takes both passes and `GROUP_PLAN`
//! every family, for exploration and for `--list`.
//!
//! Curve fixtures come from known-answer-test certificates; inputs are N
//! synthetic 36-byte items, the size of outpoints. binary.122 and binary.122-gls
//! use dense and GF(2^61)-valued constants. binary-lambda variants change accumulator
//! coordinates; binary-w variants also change the hash and encoding, and
//! binary-u variants keep binary-w's curves, hashes and encodings under
//! unscaled accumulators. weier-jacobian variants use Jacobian accumulators
//! on the same curves. *.61x2, *.64x2 and *.goldilocks2 are over GF(p^2),
//! p = 2^61 - 1, 2^64 - 59 and 2^64 - 2^32 + 1, with curves from the
//! `select_fp2` certificates; twisted.* are a = -1 curves with twisted.128's
//! model and codec.

mod common;
use common::each;

use common::h2c::{Constructions, PROJECTIONS};
use criterion::measurement::WallTime;
use criterion::{BenchmarkGroup, Criterion, Throughput, criterion_group, criterion_main};
use ephemeral_ecmh::ecmh::TAG_ITEM;
use ephemeral_ecmh::group::{Decode, Group, HashToCurve, Negate};
use ephemeral_ecmh::hash::{self, Salted};
use std::collections::HashSet;
use std::hint::black_box;

const N: usize = 1024;
/// Batch sizes for mode=batch.
const BATCHES: [usize; 3] = [8, 64, N];
/// Independent accumulators for mode=throughput.
const ACCS: usize = 8;

/// The rows `GROUP_PASS` and `GROUP_PLAN` select.
struct Pass {
    core: bool,
    wide: bool,
    /// The families whose wide rows run; None for every family.
    plan: Option<HashSet<String>>,
}

impl Pass {
    fn new() -> Self {
        let pass = std::env::var("GROUP_PASS").ok();
        let (core, wide) = match pass.as_deref() {
            None => (true, true),
            Some("core") => (true, false),
            Some("wide") => (false, true),
            Some(p) => panic!("GROUP_PASS: {p:?}, not core or wide"),
        };
        let plan = std::env::var_os("GROUP_PLAN").map(|path| {
            let text = std::fs::read_to_string(&path)
                .unwrap_or_else(|e| panic!("GROUP_PLAN={}: {e}", path.to_string_lossy()));
            text.lines()
                .filter(|l| !l.is_empty() && !l.starts_with('#'))
                .map(|line| match line.split('\t').collect::<Vec<_>>()[..] {
                    [family, "lead" | "contrast"] => family.to_string(),
                    _ => panic!("GROUP_PLAN: expected family and role in {line:?}"),
                })
                .collect()
        });
        Self { core, wide, plan }
    }

    fn wide(&self, family: &str) -> bool {
        self.wide && self.plan.as_ref().is_none_or(|p| p.contains(family))
    }
}

/// Process all xs in chunks of each size in `sizes`. Criterion records
/// the iteration time with xs.len() elements of throughput; the report
/// normalizes per element to compare the effect of batch size.
/// An encoder can batch a whole set; a decoder has the candidates in one
/// peeling wave, often far fewer.
fn batched<T, R>(
    g: &mut BenchmarkGroup<WallTime>,
    name: &str,
    sizes: &[usize],
    xs: &[T],
    f: impl Fn(&[T]) -> R,
) {
    g.throughput(Throughput::Elements(xs.len() as u64));
    for &n in sizes {
        g.bench_function(format!("{name}/mode=batch,n={n}"), |bn| {
            bn.iter(|| {
                for ch in xs.chunks(n) {
                    black_box(f(black_box(ch)));
                }
            })
        });
    }
}

/// Folds `step` over xs, one accumulator (latency) and ACCS (throughput).
fn accumulate<G: Group, T>(
    g: &mut BenchmarkGroup<WallTime>,
    name: &str,
    group: &G,
    xs: &[T],
    latency: bool,
    step: impl Fn(&G::Point, &T) -> G::Point,
) {
    g.throughput(Throughput::Elements(xs.len() as u64));
    if latency {
        g.bench_function(format!("{name}/mode=latency"), |bn| {
            bn.iter(|| xs.iter().fold(group.identity(), |acc, a| step(&acc, a)))
        });
    }
    g.bench_function(format!("{name}/mode=throughput"), |bn| {
        bn.iter(|| {
            let mut acc = [group.identity(); ACCS];
            for ch in xs.chunks(ACCS) {
                for (p, a) in acc.iter_mut().zip(ch) {
                    *p = step(p, a);
                }
            }
            acc
        })
    });
}

fn family<G: HashToCurve + Negate + Decode + Constructions>(
    c: &mut Criterion,
    pass: &Pass,
    fixture: common::families::Fixture<G>,
) {
    let name = fixture.family.id();
    let name = name.as_str();
    let group = fixture.group;
    let seed = &fixture.seed;
    let r = fixture.r;
    let info = fixture.family.curve;
    let (core, wide) = (pass.core, pass.wide(name));
    if !(core || wide) {
        return;
    }
    let sizes: Vec<usize> = BATCHES
        .into_iter()
        .filter(|&n| if n == N { core } else { wide })
        .collect();
    let sizes = sizes.as_slice();
    // Parameters belong to the verified fixture, not to its display name.
    if core {
        eprintln!(
            "group-fixture\t{name}\tr={r}\tcofactor={}\tautomorphisms={}",
            info.cofactor, info.automorphisms
        );
    }
    let salt = Salted::new(TAG_ITEM, seed);
    let items = common::items(&[], N);
    let refs: Vec<&[u8]> = items.iter().map(|x| x.as_slice()).collect();
    let hs = group.hash_batch(&salt, &refs);
    let adds = group.prepare_batch(&hs);
    // Running sums provide accumulator inputs for identity and codec timings.
    let pts: Vec<G::Point> = adds
        .iter()
        .scan(group.identity(), |p, a| {
            *p = group.add(p, a);
            Some(*p)
        })
        .collect();

    let mut g = c.benchmark_group("h2c");
    if core {
        each(&mut g, &format!("{name}/mode=indep"), &refs, |m| {
            group.hash(&salt, m)
        });
    }
    batched(&mut g, name, sizes, &refs, |ms| group.hash_batch(&salt, ms));
    g.finish();

    if wide {
        let ids: Vec<[u8; 32]> = items.iter().map(|x| hash::id(x)).collect();
        let ids: Vec<&[u8]> = ids.iter().map(|x| x.as_slice()).collect();
        let hashes = [("none", Salted::new(TAG_ITEM, seed))]
            .into_iter()
            .chain(PROJECTIONS.map(|(t, f)| (t, Salted::projected(TAG_ITEM, seed, f))));
        let mut g = c.benchmark_group("h2c.id");
        g.throughput(Throughput::Elements(N as u64));
        for (proj, h) in hashes {
            for k in group.constructions() {
                let id = format!("{name}/proj={proj},h2c={},mode=batch,n={N}", k.token);
                g.bench_function(id, |bn| bn.iter(|| (k.batch)(&group, &h, black_box(&ids))));
            }
        }
        g.finish();
    }

    let mut g = c.benchmark_group("group.prepare");
    if core {
        each(&mut g, &format!("{name}/mode=indep"), &hs, |a| {
            group.prepare(a)
        });
    }
    batched(&mut g, name, sizes, &hs, |a| group.prepare_batch(a));
    g.finish();

    if core {
        let mut g = c.benchmark_group("group.add");
        accumulate(&mut g, name, &group, &adds, true, |p, a| group.add(p, a));
        g.finish();
    }

    if wide {
        let mut g = c.benchmark_group("group.sub");
        accumulate(&mut g, name, &group, &adds, false, |p, a| {
            group.add(p, &group.neg_addend(a))
        });
        g.finish();
    }

    if core {
        let mut g = c.benchmark_group("group.neg");
        each(&mut g, name, &adds, |a| group.neg_addend(a));
        g.finish();

        let mut g = c.benchmark_group("group.is_identity");
        each(&mut g, name, &pts, |p| group.is_identity(p));
        g.finish();
    }

    // A pure cell's sum equals the candidate's addend, but as an accumulator
    // two updates away from it, so with a scale that is not one; an impure
    // cell's sum differs, and comparisons can stop at the first coordinate.
    let pure: Vec<(G::Point, G::Addend)> = adds
        .iter()
        .zip(adds.iter().cycle().skip(1))
        .map(|(a, b)| {
            let p = group.add(&group.add(&group.identity(), a), b);
            (group.add(&p, &group.neg_addend(b)), *a)
        })
        .collect();
    let impure: Vec<(G::Point, G::Addend)> = pts
        .iter()
        .copied()
        .zip(adds.iter().copied().cycle().skip(1))
        .collect();
    assert!(
        pure.iter().all(|(p, a)| group.equals_addend(p, a)),
        "{name}: impure match row"
    );
    assert!(
        !impure.iter().any(|(p, a)| group.equals_addend(p, a)),
        "{name}: pure mismatch row"
    );
    let mut g = c.benchmark_group("group.equals");
    if wide {
        each(&mut g, &format!("{name}/mode=match"), &pure, |(p, a)| {
            group.equals_addend(p, a)
        });
    }
    if core {
        each(
            &mut g,
            &format!("{name}/mode=mismatch"),
            &impure,
            |(p, a)| group.equals_addend(p, a),
        );
    }
    g.finish();

    let mut g = c.benchmark_group("group.encode");
    if core {
        each(&mut g, &format!("{name}/mode=indep"), &pts, |p| {
            group.encode_point(p)
        });
    }
    batched(&mut g, name, sizes, &pts, |ps| group.encode_batch(ps));
    g.finish();

    let es = group.encode_batch(&pts);
    let mut g = c.benchmark_group("group.decode");
    if core {
        each(&mut g, &format!("{name}/mode=indep"), &es, |e| {
            group.decode(e)
        });
    }
    batched(&mut g, name, sizes, &es, |es| group.decode_batch(es));
    g.finish();
}

fn groups(c: &mut Criterion) {
    struct Bench<'a>(&'a mut Criterion, Pass);
    impl common::families::Visitor for Bench<'_> {
        fn visit<G: HashToCurve + Negate + Decode + Constructions>(
            &mut self,
            fixture: common::families::Fixture<G>,
        ) {
            family(self.0, &self.1, fixture);
        }
    }
    common::families::Fixtures::new().visit(&mut Bench(c, Pass::new()));
}

criterion_group! {
    name = benches;
    config = common::config();
    targets = groups
}
criterion_main!(benches);
