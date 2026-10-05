//! Finite and rateless RIBLT workloads, used to compare the combined costs
//! of hashing, reusable addends, and signed accumulation.
//!
//!   nix develop -c cargo bench --bench riblt
//!
//! In the repeated reconciliation of docs/workload.md, riblt.encode per
//! item stands for a new item's cost to a party that keeps its cells, though
//! it also pays for allocating m cells, amortized over the n items, which
//! kept cells do not. riblt.cells per item stands for a rebuilt round's cost
//! per item of the set, including that allocation, the mapping and the
//! bookkeeping. riblt.peel is a round's decoding.
//!
//! Hashing is a parameter, not a property of the curve. The group suite
//! measures each curve's own hash, the comparison suites the alternatives,
//! and the rows that hash (encode, peel, stream) carry `h2c=<token>` for
//! the construction they ran under, as `benches/common/h2c.rs` lists them
//! per family: try-and-increment (`ti`, the groups' own hash; `ti-addend`
//! straight to the addend), Pornin's map (`pornin`, `pornin-addend`),
//! Elligator 2 and SSWU. riblt.cells prepares its addends outside the
//! timed routine and runs once per family without the parameter.
//!
//! Coverage is planned from the elementary measurements, not fixed here:
//! the TSV `RIBLT_PLAN` names gives each curve family its constructions, a
//! scope each and, for some, projection fields for the checksum and the
//! mapping. `full` sweeps every
//! dimension: encode and cells at every m, peel at every d under all four
//! `Peel` settings, stream at every difference. `buffer` keeps each swept
//! dimension's endpoints and the spot point, peel with the prefilter and
//! batched hashing only. `spot` is one point per dimension: encode and
//! cells at m = 1350 with n = 3500, the summary's insertion workload; peel
//! at d = 1350 with the prefilter and batched hashing; stream at d = 1000.
//! A family the plan omits is spot-checked under every construction it has;
//! with no plan, every family runs every construction in full. A line with
//! projection fields adds ID-keyed rows (`proj=<field>,mapproj=<field>`):
//! encode and peel at the buffer points under `Riblt::projected`, over the
//! same items' 32-byte IDs, which are computed outside the timed routine
//! since an ID outlives every salt. The mapping's field is chosen from
//! riblt.mapping, which bench-run times first (`RIBLT_PASS=mapping`) and
//! the rest after it (`RIBLT_PASS=workload`); unset, `RIBLT_PASS` takes
//! both. The reference groups xor-sha256.64, xor-siphash.64
//! and ristretto255 are not planned: they have one hash each, no parameter,
//! and keep the sweeps Go's benchmarks have, encode at every m and stream
//! at every difference, with one point of the rest.
//!
//! - riblt.encode: allocate m cells and insert n items, including hashing,
//!   preparation, mapping, key XOR, counts, and point updates. The expected
//!   mapping degree is k(m) = 2(H_{m+1} - 1), approximately 2.90, 5.29,
//!   9.20, 13.57, and 17.96 for m = 5, 20, 150, 1350, and 12150.
//!   The largest table exceeds a 128 KiB L1 data cache from its item keys
//!   alone; point sums, counts, and padding increase its working set.
//!   Go's BenchmarkSketchAddSymbol corresponds to these rows per item: each
//!   of its operations hashes a 64-byte symbol to a ristretto255 point and
//!   adds it to the k(m) cells of a preallocated sketch, for m from 10^3 to
//!   10^7 (k(m) from 12.97 to 31.39). Here items have 36 bytes, the m cells
//!   are allocated within the timed routine, and m = 1350 (k = 13.57) is the
//!   row nearest Go's m = 1000. Go's sketches from m = 10^5 exceed the
//!   largest table here by one to three orders of magnitude, so they include
//!   memory effects that no row here measures.
//! - riblt.cells: start with prepared addends; retain allocation, mapping,
//!   key XOR, counts, and point updates.
//! - riblt.peel: decode d differences under the `Peel` settings (prefilter
//!   and batch, each on or off), using the first successful tested cell
//!   count in 2d, 2d + d/2, ... . Fixture construction and cloning are
//!   outside the timed decoding operation. Small differences can require
//!   proportionally more cells.
//! - riblt.stream: Go's BenchmarkEncodeAndDecode over `riblt::rateless`. The
//!   decoder holds d/2 items of its own and d common ones, the encoder d/2 of
//!   its own and the d common ones, all 64-byte symbols with a little-endian
//!   counter, fresh for every iteration; the counter starts at 0 in every row,
//!   so the families reconcile the same sequence of sets. Insertion, which
//!   hashes and prepares every item, is the untimed setup; the timed routine
//!   produces, receives and tries to decode one coded symbol at a time until
//!   the decoder has decoded every received symbol. Rows are normalized per
//!   difference, where Go reports time per reconciliation. Go's symbols/diff
//!   is a property of the RIBLT, not of the checksum, and is not recorded.
//!   d ranges over Go's 10, 20, 40, 100, 1000, 10000, 50000 and 100000.
//!   The group takes 10 samples in flat sampling mode, so that a row of
//!   d = 100000 can run one reconciliation per sample.
//! - riblt.mapping: the mapping's index generators alone, per index of one
//!   mapping advanced 1024 times (`next`, as Go's BenchmarkMapping), and per
//!   item from its map digest through the indices below m (`item`), at
//!   every m; and the salted map digest of an item, which every insertion
//!   computes besides its hash to the curve (`salted-sha256/digest`). Per
//!   item, the unsalted SHA-256 ID that outlives salts (`sha256/id`), and
//!   the per-salt map digest of that 32-byte ID: salted SHA-256, or its
//!   projection in F_(2^130 - 5), F_(2^127 - 1) or GF(2^127)
//!   (`<hash>/digest of id`), which is a projected mapping's seed; and per
//!   salt, what each hash derives before its first item (`<hash>/keys`).
//!   xoshiro256pp is the default, chacha8 riblt-ecmh's generator, mcg64
//!   upstream's 64-bit generator, and sha256-ctr a SHA-256 counter stream
//!   over the digest, a third wide-seed construction.
//!
//! Every workload row uses the xoshiro256++ mapping. riblt.encode and
//! riblt.stream also run xor-siphash.64, ristretto255 and binary.127 under
//! upstream's Mcg64 mapping, with `map=mcg64` in the parameter, which
//! isolates the mapping's share where hashing costs least, on the
//! reference's group and on a curve family. riblt.stream runs the same
//! three under ChaCha8, `map=chacha8`, at the two largest differences,
//! where the window holds one generator state per item and the state's
//! size shows. ristretto255 under ChaCha8 is riblt-ecmh's own combination,
//! the matched comparison with the Go numbers.
//!
//! The curve families are those of `benches/group.rs`. The binary and
//! binary-lambda variants share curves and hashes while using different
//! accumulator coordinates. The binary-w variants use their own hash and
//! codec and need no addend conversion; the binary-u variants keep binary-w's
//! curves, hashes and cells under unscaled accumulators. The weier-jacobian
//! variants share weier's curves and hashes under Jacobian accumulators. The
//! twisted families accumulate in the quotient by the rational 2-torsion;
//! twisted.128 is that construction over F_p, p = 2^128 - 275.
//!
//! ristretto255 is curve25519-dalek's group through
//! benches/common/ristretto.rs, the checksum group of Buchanan's riblt-ecmh
//! <https://github.com/DavidBuchanan314/riblt-ecmh>; it hashes salted
//! items, and has no prepared addend form.
//!
//! xor-sha256.64 and xor-siphash.64 use 64-bit XOR checksums through the same
//! encoder and decoder. The first uses salted SHA-256, the second SipHash-2-4
//! keyed per seed. The XOR of 64-bit item hashes is the RIBLT checksum that
//! <https://github.com/yangl1996/riblt/issues/3> proposes replacing.
//! They are computational reference points for the point-sum workload,
//! not lower bounds. Curve fixtures come from known-answer-test
//! certificates. Encoder rows are normalized per item; peeling rows per input
//! difference.

mod common;

use std::collections::HashMap;

use common::h2c::{self, Hashed, Scope};
use criterion::measurement::WallTime;
use criterion::{
    BatchSize, BenchmarkGroup, BenchmarkId, Criterion, SamplingMode, Throughput, criterion_group,
    criterion_main,
};
use ephemeral_ecmh::group::{Accumulate, Decode, Encode, Group, HashToCurve, Negate, SumBatch};
use ephemeral_ecmh::hash::{self, Field, Salted};
use ephemeral_ecmh::riblt::rateless::{Decoder, Encoder};
use ephemeral_ecmh::riblt::{ChaCha8, Mapping, Mcg64, Peel, Prng, Riblt, TAG_MAP, Xoshiro256pp};
use sha2::{Digest, Sha256 as Sha256Hasher};
use siphasher::sip::SipHasher24;
use std::cell::Cell;
use std::hint::black_box;

/// Synthetic 36-byte items per encoder input, a representative block's
/// transaction count.
const N: usize = 3500;
const MS: [usize; 5] = [5, 20, 150, 1350, 12150];
const DS: [usize; 5] = [4, 20, 150, 1350, 10000];
/// Go's BenchmarkEncodeAndDecode differences.
const STREAM_DS: [usize; 8] = [10, 20, 40, 100, 1000, 10000, 50000, 100000];
/// The differences the `map=chacha8` comparison runs.
const CHACHA8_DS: [usize; 2] = [10000, 100000];

/// Which points of each dimension a family measures.
struct Sweep {
    /// encode's cell counts.
    ms: &'static [usize],
    /// cells' cell counts.
    cell_ms: &'static [usize],
    /// peel's differences.
    ds: &'static [usize],
    /// peel under all four `Peel` settings, or the prefiltered batched one.
    ablation: bool,
    /// stream's differences.
    stream: &'static [usize],
}

/// Every point of every dimension.
const FULL: Sweep = Sweep {
    ms: &MS,
    cell_ms: &MS,
    ds: &DS,
    ablation: true,
    stream: &STREAM_DS,
};
/// Each dimension's endpoints and its spot point.
const BUFFER: Sweep = Sweep {
    ms: &[5, 1350, 12150],
    cell_ms: &[5, 1350, 12150],
    ds: &[4, 1350, 10000],
    ablation: false,
    stream: &[10, 1000, 100000],
};
/// One point per dimension.
const SPOT: Sweep = Sweep {
    ms: &[1350],
    cell_ms: &[1350],
    ds: &[1350],
    ablation: false,
    stream: &[1000],
};
/// The reference groups: Go's sweeps, one point of the rest.
const PARITY: Sweep = Sweep {
    ms: &MS,
    cell_ms: &[1350],
    ds: &[1350],
    ablation: false,
    stream: &STREAM_DS,
};

impl Scope {
    fn sweep(self) -> &'static Sweep {
        match self {
            Scope::Full => &FULL,
            Scope::Buffer => &BUFFER,
            Scope::Spot => &SPOT,
        }
    }
}

/// XOR of 64-bit item hashes produced by `H`.
#[derive(Clone, Copy)]
struct Xor64<H>(H);

/// A 64-bit item hash for `Xor64`.
trait Hash64: Copy {
    fn hash64(&self, h: &Salted, msg: &[u8]) -> u64;
}

/// The first 8 bytes of the salted SHA-256 the curves hash with.
#[derive(Clone, Copy)]
struct Sha256;

impl Hash64 for Sha256 {
    fn hash64(&self, h: &Salted, msg: &[u8]) -> u64 {
        u64::from_le_bytes(h.digest(msg, 0)[..8].try_into().unwrap())
    }
}

/// SipHash-2-4 with a key derived once from the fixture seed.
#[derive(Clone, Copy)]
struct Sip(SipHasher24);

impl Sip {
    fn new(seed: &[u8; 32]) -> Self {
        let key = Salted::new(b"bench-siphash-key", seed).digest(&[], 0);
        Sip(SipHasher24::new_with_key(key[..16].try_into().unwrap()))
    }
}

impl Hash64 for Sip {
    fn hash64(&self, _: &Salted, msg: &[u8]) -> u64 {
        self.0.hash(msg)
    }
}

impl<H: Hash64> Group for Xor64<H> {
    type Affine = u64;
    type Point = u64;
    fn identity(&self) -> u64 {
        0
    }
    fn is_identity(&self, p: &u64) -> bool {
        *p == 0
    }
    fn to_affine(&self, p: &u64) -> u64 {
        *p
    }
}

impl<H: Hash64> Accumulate for Xor64<H> {
    type Addend = u64;
    fn prepare(&self, a: &u64) -> u64 {
        *a
    }
    fn add(&self, p: &u64, a: &u64) -> u64 {
        p ^ a
    }
    fn add_affine(&self, p: &u64, a: &u64) -> u64 {
        p ^ a
    }
}

impl<H: Hash64> Negate for Xor64<H> {
    fn neg(&self, a: &u64) -> u64 {
        *a
    }
    fn neg_addend(&self, a: &u64) -> u64 {
        *a
    }
}

impl<H: Hash64> HashToCurve for Xor64<H> {
    fn hash(&self, h: &Salted, msg: &[u8]) -> u64 {
        self.0.hash64(h, msg)
    }
}

impl<H: Hash64> SumBatch for Xor64<H> {
    fn sum_batch(&self, a: &[u64]) -> u64 {
        a.iter().fold(0, |s, x| s ^ x)
    }
}

impl<H: Hash64> Encode for Xor64<H> {
    type Encoding = [u8; 8];
    fn encode(&self, a: &u64) -> [u8; 8] {
        a.to_le_bytes()
    }
}

/// SHA-256 in counter mode keyed by the map digest, four samples per
/// compression: the third wide-seed generator in riblt.mapping.
#[derive(Clone)]
struct Sha256Ctr {
    key: [u8; 32],
    ctr: u32,
    buf: [u64; 4],
    pos: usize,
}

impl Prng for Sha256Ctr {
    fn seed(digest: &[u8; 32]) -> Self {
        Self {
            key: *digest,
            ctr: 0,
            buf: [0; 4],
            pos: 4,
        }
    }

    fn next_u64(&mut self) -> u64 {
        if self.pos == 4 {
            let h: [u8; 32] = Sha256Hasher::new()
                .chain_update(self.key)
                .chain_update(self.ctr.to_le_bytes())
                .finalize()
                .into();
            for (b, w) in self.buf.iter_mut().zip(h.as_chunks::<8>().0) {
                *b = u64::from_le_bytes(*w);
            }
            self.ctr += 1;
            self.pos = 0;
        }
        self.pos += 1;
        self.buf[self.pos - 1]
    }
}

/// The cell keys of a workload, distinct per tag: outpoint-sized items,
/// or for the ID-keyed rows those items' IDs.
trait Keys: Sized {
    fn keys(tag: u32, n: usize) -> Vec<Self>;
}

impl Keys for [u8; 36] {
    fn keys(tag: u32, n: usize) -> Vec<Self> {
        common::items(&tag.to_le_bytes(), n)
    }
}

impl Keys for [u8; 32] {
    fn keys(tag: u32, n: usize) -> Vec<Self> {
        <[u8; 36]>::keys(tag, n)
            .iter()
            .map(|x| hash::id(x))
            .collect()
    }
}

/// `suffix` is appended to the parameter: `,h2c=<token>` for the hash
/// construction, `,map=<generator>` for a mapping other than the default
/// and `,proj=<field>,mapproj=<field>` for ID-keyed cells under
/// projections.
fn encode<G: HashToCurve + Negate, P: Prng, const L: usize>(
    g: &mut BenchmarkGroup<WallTime>,
    name: &str,
    suffix: &str,
    r: &Riblt<G, P>,
    xs: &[[u8; L]],
    ms: &[usize],
) {
    g.throughput(Throughput::Elements(xs.len() as u64));
    for &m in ms {
        let id = BenchmarkId::new(name, format!("m={m},n={}{suffix}", xs.len()));
        g.bench_with_input(id, &m, |b, &m| {
            b.iter(|| {
                let mut cells = r.cells(m);
                r.encode(&mut cells, black_box(xs), 1);
                cells
            })
        });
    }
}

fn cells<G: HashToCurve + Negate, P: Prng>(
    g: &mut BenchmarkGroup<WallTime>,
    name: &str,
    r: &Riblt<G, P>,
    xs: &[[u8; 36]],
    ms: &[usize],
) {
    let adds = r.addends(xs);
    g.throughput(Throughput::Elements(xs.len() as u64));
    for &m in ms {
        let id = BenchmarkId::new(name, format!("m={m},n={}", xs.len()));
        g.bench_with_input(id, &m, |b, &m| {
            b.iter(|| {
                let mut cells = r.cells(m);
                for (x, a) in xs.iter().zip(&adds) {
                    r.apply(&mut cells, x, a, 1);
                }
                cells
            })
        });
    }
}

/// d differences (half each side) over 1000 common items.
fn peel<G: HashToCurve + Negate, P: Prng, const L: usize>(
    g: &mut BenchmarkGroup<WallTime>,
    name: &str,
    suffix: &str,
    r: &Riblt<G, P>,
    ds: &[usize],
    ablation: bool,
) where
    [u8; L]: Keys + Ord,
{
    let common = <[u8; L]>::keys(0, 1000);
    let settings: &[Peel] = if ablation {
        &[
            Peel {
                prefilter: false,
                batch: false,
            },
            Peel {
                prefilter: false,
                batch: true,
            },
            Peel {
                prefilter: true,
                batch: false,
            },
            Peel {
                prefilter: true,
                batch: true,
            },
        ]
    } else {
        &[Peel {
            prefilter: true,
            batch: true,
        }]
    };
    for &d in ds {
        let (a, b) = (<[u8; L]>::keys(1, d / 2), <[u8; L]>::keys(2, d - d / 2));
        let (a, b) = ([common.clone(), a].concat(), [common.clone(), b].concat());
        let diff = (2 * d..)
            .step_by(d / 2)
            .map(|m| {
                let mut diff = r.cells(m);
                r.encode(&mut diff, &a, 1);
                r.encode(&mut diff, &b, -1);
                diff
            })
            .find(|diff| {
                r.peel(
                    &mut diff.clone(),
                    Peel {
                        prefilter: false,
                        batch: true,
                    },
                )
                .is_some()
            })
            .unwrap();
        let m = diff.len();
        g.throughput(Throughput::Elements(d as u64));
        for &how in settings {
            let Peel { prefilter, batch } = how;
            assert!(
                r.peel(&mut diff.clone(), how).is_some(),
                "{name} d={d}: no decode"
            );
            let id = BenchmarkId::new(
                name,
                format!("d={d},m={m},prefilter={prefilter},batch={batch}{suffix}"),
            );
            g.bench_with_input(id, &how, |bn, &how| {
                bn.iter_batched_ref(|| diff.clone(), |c| r.peel(c, how), BatchSize::SmallInput)
            });
        }
    }
}

/// Go's test symbol: a little-endian counter in 64 bytes.
fn symbol(i: u64) -> [u8; 64] {
    let mut x = [0; 64];
    x[..8].copy_from_slice(&i.to_le_bytes());
    x
}

/// Go's BenchmarkEncodeAndDecode: `suffix` is appended to the parameter.
fn stream<G: HashToCurve + Negate, P: Prng>(
    g: &mut BenchmarkGroup<WallTime>,
    name: &str,
    suffix: &str,
    r: &Riblt<G, P>,
    ds: &[usize],
) {
    for &d in ds {
        let half = d as u64 / 2;
        let next = Cell::new(0u64);
        let setup = || {
            let id = next.get();
            next.set(id + 4 * half);
            let own = |from| (from..from + half).map(symbol);
            let common: Vec<_> = (id + 2 * half..id + 4 * half).map(symbol).collect();
            let mut enc = Encoder::new(r);
            let mut dec = Decoder::new(r);
            enc.insert(&own(id).chain(common.iter().copied()).collect::<Vec<_>>());
            dec.insert(&own(id + half).chain(common).collect::<Vec<_>>());
            (enc, dec)
        };
        g.throughput(Throughput::Elements(d as u64));
        let id = BenchmarkId::new(name, format!("d={d}{suffix}"));
        g.bench_function(id, |b| {
            b.iter_batched(
                setup,
                |(mut enc, mut dec)| {
                    loop {
                        dec.receive(enc.produce());
                        dec.try_decode();
                        if dec.decoded() {
                            break;
                        }
                    }
                    (enc, dec)
                },
                BatchSize::LargeInput,
            )
        });
    }
}

fn riblt(c: &mut Criterion) {
    use common::families::{Fixture, Fixtures, Visitor};
    // Curve certificates use their own seeds; every RIBLT row uses
    // binary127's seed as the workload salt.
    let fixtures = Fixtures::new();
    let seed = fixtures.binary127.seed;
    let xs = <[u8; 36]>::keys(0, N);
    let ids = <[u8; 32]>::keys(0, N);
    let plan = h2c::plan();
    let (mapping_rows, workload) = match std::env::var("RIBLT_PASS").ok().as_deref() {
        None => (true, true),
        Some("mapping") => (true, false),
        Some("workload") => (false, true),
        Some(p) => panic!("RIBLT_PASS: {p:?}, not mapping or workload"),
    };

    enum Stage {
        Encode,
        Cells,
        Peel,
        Stream,
    }
    struct Bench<'a, 'b> {
        group: &'a mut BenchmarkGroup<'b, WallTime>,
        stage: Stage,
        seed: &'a [u8; 32],
        items: &'a [[u8; 36]],
        ids: &'a [[u8; 32]],
        plan: Option<&'a HashMap<String, Vec<h2c::Line>>>,
    }
    impl Bench<'_, '_> {
        fn measure<G: HashToCurve + Negate, P: Prng>(
            &mut self,
            name: &str,
            suffix: &str,
            r: &Riblt<G, P>,
            sweep: &Sweep,
        ) {
            match self.stage {
                Stage::Encode => encode(self.group, name, suffix, r, self.items, sweep.ms),
                Stage::Cells => cells(self.group, name, r, self.items, sweep.cell_ms),
                Stage::Peel => {
                    peel::<_, _, 36>(self.group, name, suffix, r, sweep.ds, sweep.ablation)
                }
                Stage::Stream => stream(self.group, name, suffix, r, sweep.stream),
            }
        }

        /// The mapping-generator rows, which only encode and stream measure.
        fn generators<G: HashToCurve + Negate>(&mut self, name: &str, suffix: &str, group: G) {
            let mcg = Riblt::<_, Mcg64>::with_prng(group, self.seed);
            let mcg_map = format!("{suffix},map=mcg64");
            match self.stage {
                Stage::Encode => encode(self.group, name, &mcg_map, &mcg, self.items, &MS),
                Stage::Stream => {
                    stream(self.group, name, &mcg_map, &mcg, &STREAM_DS);
                    let chacha = Riblt::<_, ChaCha8>::with_prng(group, self.seed);
                    let chacha_map = format!("{suffix},map=chacha8");
                    stream(self.group, name, &chacha_map, &chacha, &CHACHA8_DS);
                }
                Stage::Cells | Stage::Peel => {}
            }
        }

        /// ID-keyed cells, the checksum projecting in `item` and the mapping
        /// in `map`, at the buffer points of encode and peel, the stages
        /// that hash and seed mappings.
        fn projected<G: HashToCurve + Negate>(
            &mut self,
            name: &str,
            suffix: &str,
            group: G,
            (item, map): (Field, Field),
        ) {
            let r = Riblt::projected(group, self.seed, item, map);
            let token = |field| {
                h2c::PROJECTIONS
                    .iter()
                    .find(|(_, f)| *f == field)
                    .unwrap()
                    .0
            };
            let suffix = format!("{suffix},proj={},mapproj={}", token(item), token(map));
            match self.stage {
                Stage::Encode => encode(self.group, name, &suffix, &r, self.ids, BUFFER.ms),
                Stage::Peel => peel::<_, _, 32>(self.group, name, &suffix, &r, BUFFER.ds, false),
                Stage::Cells | Stage::Stream => {}
            }
        }

        /// A reference group: Go's sweeps and the generator rows.
        fn reference<G: HashToCurve + Negate>(&mut self, name: &str, group: G, generators: bool) {
            self.measure(name, "", &Riblt::new(group, self.seed), &PARITY);
            if generators {
                self.generators(name, "", group);
            }
        }
    }
    impl Visitor for Bench<'_, '_> {
        fn visit<G: HashToCurve + Negate + Decode + h2c::Constructions>(
            &mut self,
            fixture: Fixture<G>,
        ) {
            let name = fixture.family.id();
            let planned = h2c::planned(&fixture.group, &name, self.plan);
            for (construction, scope, fields) in &planned {
                let hashed = Hashed {
                    group: fixture.group,
                    h2c: construction,
                };
                let sweep = scope.sweep();
                // cells prepares its addends outside the timed routine: one
                // row per family, under no construction, at the widest
                // planned scope
                if let Stage::Cells = self.stage {
                    let widest = planned.iter().map(|(_, s, _)| *s).min().unwrap();
                    cells(
                        self.group,
                        &name,
                        &Riblt::new(hashed, self.seed),
                        self.items,
                        widest.sweep().cell_ms,
                    );
                    break;
                }
                let suffix = format!(",h2c={}", construction.token);
                self.measure(&name, &suffix, &Riblt::new(hashed, self.seed), sweep);
                for &fields in fields {
                    self.projected(&name, &suffix, hashed, fields);
                }
                if name == "binary.127" {
                    self.generators(&name, &suffix, hashed);
                }
            }
        }
    }
    if workload {
        for (name, stage) in [
            ("riblt.encode", Stage::Encode),
            ("riblt.cells", Stage::Cells),
            ("riblt.peel", Stage::Peel),
            ("riblt.stream", Stage::Stream),
        ] {
            let mut group = c.benchmark_group(name);
            if let Stage::Stream = stage {
                group.sample_size(10).sampling_mode(SamplingMode::Flat);
            }
            let mut bench = Bench {
                group: &mut group,
                stage,
                seed: &seed,
                items: &xs,
                ids: &ids,
                plan: plan.as_ref(),
            };
            bench.reference("xor-sha256.64", Xor64(Sha256), false);
            bench.reference("xor-siphash.64", Xor64(Sip::new(&seed)), true);
            bench.reference("ristretto255", common::ristretto::Ristretto255, true);
            fixtures.visit(&mut bench);
            group.finish();
        }
    }
    if !mapping_rows {
        return;
    }

    let mut group = c.benchmark_group("riblt.mapping");
    let rx = Riblt::new(Xor64(Sha256), &seed);
    let digests: Vec<[u8; 32]> = xs.iter().map(|x| rx.map.digest(x, 0)).collect();
    group.throughput(Throughput::Elements(xs.len() as u64));
    group.bench_function(BenchmarkId::new("salted-sha256", "digest"), |b| {
        b.iter(|| {
            black_box(&xs)
                .iter()
                .fold(0, |s, x| s ^ rx.map.digest(x, 0)[0])
        })
    });
    // the per-salt seed of an ID, the item's unsalted SHA-256, against the
    // ID's own cost, which every salt shares
    group.bench_function(BenchmarkId::new("sha256", "id"), |b| {
        b.iter(|| black_box(&xs).iter().fold(0, |s, x| s ^ hash::id(x)[0]))
    });
    let new = |s: &[u8; 32], field: Option<Field>| match field {
        None => Salted::new(TAG_MAP, s),
        Some(f) => Salted::projected(TAG_MAP, s, f),
    };
    let hashes = [("salted-sha256".to_string(), None)]
        .into_iter()
        .chain(h2c::PROJECTIONS.map(|(t, f)| (format!("projection-{t}"), Some(f))));
    for (name, field) in hashes {
        // what a salt costs before its first item, once per tag
        group.throughput(Throughput::Elements(1));
        group.bench_function(BenchmarkId::new(&name, "keys"), |b| {
            b.iter(|| new(black_box(&seed), field))
        });
        let map = new(&seed, field);
        group.throughput(Throughput::Elements(ids.len() as u64));
        group.bench_function(BenchmarkId::new(&name, "digest of id"), |b| {
            b.iter(|| {
                black_box(&ids)
                    .iter()
                    .fold(0, |s, x| s ^ map.digest(x, 0)[0])
            })
        });
    }
    mapping::<Xoshiro256pp>(&mut group, "xoshiro256pp", &digests);
    mapping::<ChaCha8>(&mut group, "chacha8", &digests);
    mapping::<Mcg64>(&mut group, "mcg64", &digests);
    mapping::<Sha256Ctr>(&mut group, "sha256-ctr", &digests);
    group.finish();
}

/// One mapping advanced 1024 times per iteration, and each digest's
/// mapping walked below m.
fn mapping<P: Prng>(g: &mut BenchmarkGroup<WallTime>, name: &str, digests: &[[u8; 32]]) {
    const STEPS: u64 = 1024;
    g.throughput(Throughput::Elements(STEPS));
    g.bench_function(BenchmarkId::new(name, "next"), |b| {
        b.iter_batched_ref(
            || Mapping::new(P::seed(&digests[0])),
            |m| (0..STEPS).fold(0, |s, _| s ^ m.advance()),
            BatchSize::SmallInput,
        )
    });
    g.throughput(Throughput::Elements(digests.len() as u64));
    for m in MS {
        g.bench_function(BenchmarkId::new(name, format!("item,m={m}")), |b| {
            b.iter(|| {
                black_box(digests)
                    .iter()
                    .map(|d| Mapping::new(P::seed(d)).below(m).count())
                    .sum::<usize>()
            })
        });
    }
}

criterion_group! {
    name = benches;
    config = common::config();
    targets = riblt
}
criterion_main!(benches);
