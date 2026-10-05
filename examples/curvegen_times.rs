//! Wall-clock cost of selecting a curve, per seed and per family: one
//! `find`, one `prove` and one verification of the proved certificate,
//! each timed once. A search takes seconds to minutes and its length
//! depends on the seed, so it is measured single-shot rather than with
//! criterion's repeated samples; the count columns say how much work each
//! seed took, so a time can be read against them.
//!
//!   nix develop -c cargo run --release --features pari --example curvegen_times -- OUT.csv [SEEDS]
//!
//! SEEDS defaults to 4; they are the first seeds of `curvegen_stats`. One
//! CSV row per family, method and seed:
//!
//! - method: `pari+sieve` counts with PARI the candidates that the torsion
//!   sieve (`sieve`) leaves, at the family's bounds for `find` and for
//!   `prove` (`sieve::Bounds`), in the families it covers; `pari` counts every
//!   candidate with PARI and no sieve, as `prove` does by default;
//!   `agm+sieve` is the binary family in Rust alone (`find_gf2_127`,
//!   `certify_gf2_127`); `agm` counts every candidate in Rust and factors
//!   nothing: Rust alone, labelling a composite by its order; `agm+factor`
//!   (gf2_109) counts in Rust and labels by PARI's smallest prime factor,
//!   as sage/kat109.sage does, so its certificates are Sage's. The 109- and
//!   122-bit families have no sieve bounds; their `pari` rows count with
//!   PARI, factoring as `agm+factor` does for gf2_109 and not at all for
//!   the 122-bit families.
//! - candidates: indices walked, the accepted one included.
//! - find_counts, find_aborted: point counts `find` started, and how many
//!   of them SEA's early abort cut short.
//! - prove_counts: full point counts `prove` made.
//! - rejections: entries in the certificate.
//! - find_s, prove_s, verify_s: wall-clock seconds.
//!
//! Then, for `find` and for `prove` (prefixed find_, prove_), the seconds
//! spent in each phase of the search, and the work some of them did:
//!
//! - candidate_s: deriving candidates, structural exclusions included.
//! - quick_s, quick: `Family::quick_reject`, and the rejections it made.
//! - sieve_s, sieved: the sieve, and the rejections it made.
//! - count_s: point counting, aborted counts included.
//! - factor_s, factored (prove only): factoring composite orders.
//! - witness_s (prove only): rejection witnesses after a count.
//! - accept_s: verdicts on the counted orders, and the acceptance check.
//!
//! The 107- and 128-bit families have no Rust prover (their certificates
//! come from Sage), so they have no rows.

use ephemeral_ecmh::curve::binary109;
use ephemeral_ecmh::curvegen::agm::Agm;
use ephemeral_ecmh::curvegen::criteria::{self, AdmissibleR, Count, Criteria};
use ephemeral_ecmh::curvegen::pari::{self, Pari};
use ephemeral_ecmh::curvegen::prove::{
    self, Binary109, Binary127, Dense122, Edwards127, Factor, Family, Gls122, NoFactor, NoSieve,
    Rejection, Sieve, SmallL, Verdict, Weier127,
};
use ephemeral_ecmh::curvegen::select::{self, Certificate, Error, Policy};
use ephemeral_ecmh::curvegen::{select109, select122, sieve};
use ephemeral_ecmh::hash::Salted;
use std::cell::Cell;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::marker::PhantomData;
use std::time::{Duration, Instant};

/// Time and work per phase of one search.
#[derive(Clone, Copy, Debug, Default)]
struct Phases {
    candidate: Duration,
    quick: Duration,
    sieve: Duration,
    count: Duration,
    factor: Duration,
    witness: Duration,
    accept: Duration,
    quick_n: usize,
    sieved: usize,
    full: usize,
    aborted: usize,
    factored: usize,
}

thread_local! {
    static PHASES: Cell<Phases> = Cell::default();
}

fn add(f: impl FnOnce(&mut Phases)) {
    PHASES.with(|c| {
        let mut p = c.get();
        f(&mut p);
        c.set(p);
    });
}

/// Runs f, adding its time to one phase.
fn timed<R>(phase: fn(&mut Phases) -> &mut Duration, f: impl FnOnce() -> R) -> R {
    let start = Instant::now();
    let r = f();
    let t = start.elapsed();
    add(|p| *phase(p) += t);
    r
}

/// The phases since the last call.
fn take() -> Phases {
    PHASES.with(|c| c.take())
}

const FIND_COLUMNS: &str = "candidate_s,quick_s,sieve_s,count_s,accept_s,quick,sieved";
const PROVE_COLUMNS: &str =
    "candidate_s,quick_s,sieve_s,count_s,factor_s,witness_s,accept_s,quick,sieved,factored";

impl Phases {
    fn find_cells(&self) -> String {
        let s = |d: Duration| format!("{:.6}", d.as_secs_f64());
        [
            s(self.candidate),
            s(self.quick),
            s(self.sieve),
            s(self.count),
            s(self.accept),
            self.quick_n.to_string(),
            self.sieved.to_string(),
        ]
        .join(",")
    }

    fn prove_cells(&self) -> String {
        let s = |d: Duration| format!("{:.6}", d.as_secs_f64());
        [
            s(self.candidate),
            s(self.quick),
            s(self.sieve),
            s(self.count),
            s(self.factor),
            s(self.witness),
            s(self.accept),
            self.quick_n.to_string(),
            self.sieved.to_string(),
            self.factored.to_string(),
        ]
        .join(",")
    }
}

/// F, with the time of each of its hooks added to its phase.
struct Timed<F>(PhantomData<F>);

impl<F: Family> Criteria for Timed<F> {
    const TORS: u32 = F::TORS;
    fn check_order(n: u128) -> Result<u128, criteria::OrderError> {
        timed(|p| &mut p.accept, || F::check_order(n))
    }
    fn reject_without_count(c: &Self::Curve) -> bool {
        let rejected = timed(|p| &mut p.quick, || F::reject_without_count(c));
        add(|p| p.quick_n += usize::from(rejected));
        rejected
    }
    type Curve = F::Curve;
    const TAG: &'static [u8] = F::TAG;
    const R: AdmissibleR = F::R;

    fn candidate(seed: &[u8; 32], j: u32) -> Option<Self::Curve> {
        timed(|p| &mut p.candidate, || F::candidate(seed, j))
    }
}

impl<F: Family> Policy<16> for Timed<F> {
    const CLEAR: u128 = F::CLEAR;
    const ENTRY_PER_INDEX: bool = F::ENTRY_PER_INDEX;

    fn marker(c: &Self::Curve, p: &select::Affine<16, F>, l: u128) -> Option<bool> {
        F::marker(c, p, l)
    }
}

impl<F: Family> Family for Timed<F> {
    fn verdict(n: u128) -> Verdict {
        timed(|p| &mut p.accept, || F::verdict(n))
    }

    fn accept(seed: &[u8; 32], j: u32, r: u128) -> Result<Self::Curve, select::Error> {
        timed(|p| &mut p.accept, || F::accept(seed, j, r))
    }

    fn witness(c: &Self::Curve, h: &Salted, j: u32, n: u128, l: u128) -> [u8; 16] {
        timed(|p| &mut p.witness, || F::witness(c, h, j, n, l))
    }

    fn quick_reject(c: &Self::Curve) -> Option<Rejection> {
        let r = timed(|p| &mut p.quick, || F::quick_reject(c));
        add(|p| p.quick_n += usize::from(r.is_some()));
        r
    }

    fn order_witness(c: &Self::Curve, h: &Salted, j: u32, n: u128) -> Option<[u8; 16]> {
        timed(|p| &mut p.witness, || F::order_witness(c, h, j, n))
    }
}

/// Wraps a counter, timing and tallying full counts and aborts.
struct Tally<K>(K);

impl<C, K: Count<C>> Count<C> for Tally<K> {
    fn order(&mut self, c: &C) -> u128 {
        add(|p| p.full += 1);
        timed(|p| &mut p.count, || self.0.order(c))
    }

    fn order_early_abort(&mut self, c: &C, tors: u32) -> Option<u128> {
        let n = timed(|p| &mut p.count, || self.0.order_early_abort(c, tors));
        match n {
            Some(_) => add(|p| p.full += 1),
            None => add(|p| p.aborted += 1),
        }
        n
    }
}

impl<K: Factor> Factor for Tally<K> {
    fn smallest_prime_factor(&mut self, m: u128) -> Option<u128> {
        add(|p| p.factored += 1);
        timed(|p| &mut p.factor, || self.0.smallest_prime_factor(m))
    }
}

impl<C, S: Sieve<C>> Sieve<C> for Tally<S> {
    fn reject(&mut self, c: &C) -> Option<Rejection> {
        let r = timed(|p| &mut p.sieve, || self.0.reject(c));
        add(|p| p.sieved += usize::from(r.is_some()));
        r
    }
}

type Verify<C> = fn(&[u8; 32], &Certificate) -> Result<C, Error>;

/// `select109::verify` on the prover's zero-padded witnesses.
fn verify_gf2_109(seed: &[u8; 32], cert: &Certificate) -> Result<binary109::Curve, Error> {
    select109::verify(seed, &prove::certificate109(cert.clone()))
}

/// One row per seed for family F, searched with `count` and `factor`, and
/// sieved by `sieves`, find's and prove's.
#[allow(clippy::too_many_arguments)]
fn family<F: Family, S: Sieve<F::Curve>>(
    out: &mut impl Write,
    name: &str,
    method: &str,
    seeds: &[[u8; 32]],
    count: impl Count<F::Curve>,
    factor: impl Factor,
    sieves: (S, S),
    verify: Verify<F::Curve>,
) {
    let (mut count, mut factor) = (Tally(count), Tally(factor));
    let (mut find_sieve, mut prove_sieve) = (Tally(sieves.0), Tally(sieves.1));
    for (i, seed) in seeds.iter().enumerate() {
        take();
        let start = Instant::now();
        let found = criteria::find::<Timed<F>>(seed, &mut count, &mut find_sieve);
        let find_s = start.elapsed().as_secs_f64();
        let find = take();

        let start = Instant::now();
        let (cert, _) = prove::prove::<Timed<F>>(seed, &mut count, &mut factor, &mut prove_sieve);
        let prove_s = start.elapsed().as_secs_f64();
        let prove = take();
        assert_eq!((cert.index, cert.r), (found.index, found.r));

        let start = Instant::now();
        verify(seed, &cert).expect("the proved certificate verifies");
        let verify_s = start.elapsed().as_secs_f64();

        let row = format!(
            "{name},{method},{i},{},{},{},{},{},{find_s:.6},{prove_s:.6},{verify_s:.6},{},{}",
            found.index + 1,
            find.full + find.aborted,
            find.aborted,
            prove.full,
            cert.rejections.len(),
            find.find_cells(),
            prove.prove_cells(),
        );
        eprintln!("{row}");
        writeln!(out, "{row}").and_then(|()| out.flush()).unwrap();
    }
}

fn main() {
    let mut args = std::env::args().skip(1);
    let path = args.next().expect("usage: curvegen_times OUT.csv [SEEDS]");
    let seeds: u32 = args.next().map_or(4, |s| s.parse().expect("SEEDS"));
    assert!(
        pari::seadata().is_some(),
        "no seadata: SEA would compute its own modular polynomials; set GP_DATA_DIR"
    );
    let h = Salted::new(b"curvegen-stats", &[0; 32]);
    let seeds: Vec<[u8; 32]> = (0..seeds).map(|i| h.digest(&[], i)).collect();

    let mut out = BufWriter::new(File::create(&path).expect("create the CSV"));
    let prefixed = |p: &str, cols: &str| {
        cols.split(',')
            .map(|c| format!("{p}_{c}"))
            .collect::<Vec<_>>()
            .join(",")
    };
    writeln!(
        out,
        "family,method,seed,candidates,find_counts,find_aborted,prove_counts,rejections,find_s,prove_s,verify_s,{},{}",
        prefixed("find", FIND_COLUMNS),
        prefixed("prove", PROVE_COLUMNS),
    )
    .unwrap();
    let o = &mut out;
    let s = &seeds;
    let gf2 = SmallL(sieve::GF2_127_L_MAX);
    family::<Binary127, _>(
        o,
        "gf2_127",
        "agm+sieve",
        s,
        Agm,
        NoFactor,
        (gf2, gf2),
        select::verify_gf2_127,
    );
    let n = (NoSieve, NoSieve);
    let l = |b: sieve::Bounds| (SmallL(b.find), SmallL(b.prove));
    let v = select::verify_gf2_127;
    family::<Binary127, _>(o, "gf2_127", "pari", s, Pari, Pari, n, v);
    let v = verify_gf2_109;
    family::<Binary109, _>(o, "gf2_109", "agm", s, Agm, NoFactor, n, v);
    family::<Binary109, _>(o, "gf2_109", "agm+factor", s, Agm, Pari, n, v);
    family::<Binary109, _>(o, "gf2_109", "pari", s, Pari, Pari, n, v);
    let v = select122::verify_dense;
    family::<Dense122, _>(o, "gf2_122", "agm", s, Agm, NoFactor, n, v);
    family::<Dense122, _>(o, "gf2_122", "pari", s, Pari, NoFactor, n, v);
    let v = select122::verify_gls;
    family::<Gls122, _>(o, "gf2_122-gls", "agm", s, Agm, NoFactor, n, v);
    family::<Gls122, _>(o, "gf2_122-gls", "pari", s, Pari, NoFactor, n, v);
    let (v, b) = (select::verify_fp127, l(sieve::EDWARDS127));
    family::<Edwards127, _>(o, "edwards127", "pari+sieve", s, Pari, Pari, b, v);
    let (v, b) = (select::verify_weier127, l(sieve::WEIER127));
    family::<Weier127, _>(o, "weier127", "pari+sieve", s, Pari, Pari, b, v);
}
