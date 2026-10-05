//! How much work the curve search is, in counts rather than time: for 64
//! seeds per family, the candidates walked, how many `quick_reject`
//! settled uncounted (8 | #E on the Edwards side), how many were
//! point-counted, and how many of those SEA's early abort threw out
//! before finishing.
//!
//!   nix develop -c cargo run --release --features pari --example curvegen_stats
//!
//! With PARI 2.17.3 and seadata-small:
//!
//! ```text
//!  gf2_127: 64 seeds, 4599 candidates, 0 quick rejections, 4599 counted (71.9/seed, max 306), 0 aborted early, 4599 counted in full (4535 of them rejected)
//!    fp127: 64 seeds, 13159 candidates, 3324 quick rejections, 3295 counted (51.5/seed, max 327), 2454 aborted early, 841 counted in full (777 of them rejected)
//! weier127: 64 seeds, 9550 candidates, 0 quick rejections, 9550 counted (149.2/seed, max 562), 8759 aborted early, 791 counted in full (727 of them rejected)
//! ```
//!
//! `prove` counts every one of them in full. For weier127 that includes the
//! even orders, which the sieve's 2-torsion rule (a root of x^3 - 3x + b)
//! rejects without a count.

use ephemeral_ecmh::curvegen::criteria::{self, Count};
use ephemeral_ecmh::curvegen::pari::Pari;
use ephemeral_ecmh::curvegen::prove::{Binary127, Edwards127, Family, NoSieve, Weier127};
use ephemeral_ecmh::hash::Salted;

const SEEDS: u32 = 64;

/// Wraps a counter, tallying full counts and aborts.
#[derive(Default)]
struct Tally {
    full: usize,
    aborted: usize,
}

impl<C> Count<C> for Tally
where
    Pari: Count<C>,
{
    fn order(&mut self, c: &C) -> u128 {
        self.full += 1;
        Pari.order(c)
    }

    fn order_early_abort(&mut self, c: &C, tors: u32) -> Option<u128> {
        let n = Pari.order_early_abort(c, tors);
        match n {
            Some(_) => self.full += 1,
            None => self.aborted += 1,
        }
        n
    }
}

fn stats<F: Family>(name: &str)
where
    Pari: Count<F::Curve>,
{
    let h = Salted::new(b"curvegen-stats", &[0; 32]);
    let (mut walked, mut quick, mut t) = (0, 0, Tally::default());
    let mut max = 0;
    for i in 0..SEEDS {
        let seed = h.digest(&[], i);
        let before = t.full + t.aborted;
        let found = criteria::find::<F>(&seed, &mut t, &mut NoSieve);
        walked += found.index as usize + 1;
        quick += (0..found.index)
            .filter_map(|j| F::candidate(&seed, j))
            .filter(|c| F::quick_reject(c).is_some())
            .count();
        max = max.max(t.full + t.aborted - before);
    }
    let counted = t.full + t.aborted;
    println!(
        "{name:>8}: {SEEDS} seeds, {walked} candidates, {quick} quick rejections, \
         {counted} counted ({:.1}/seed, max {max}), {} aborted early, \
         {} counted in full ({} of them rejected)",
        counted as f64 / f64::from(SEEDS),
        t.aborted,
        t.full,
        t.full - SEEDS as usize,
    );
}

fn main() {
    stats::<Binary127>("gf2_127");
    stats::<Edwards127>("fp127");
    stats::<Weier127>("weier127");
}
