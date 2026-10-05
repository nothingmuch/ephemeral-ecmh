//! The prover against Sage, without PARI: point counts come from a table
//! (sage/prove_orders.sage), and must give kat.sage's certificates to the
//! byte: index, r, labels and witnesses.
//! Only candidates that `Family::quick_reject` leaves may be counted. The
//! 109- and 122-bit families count with the AGM instead.

#[path = "common/kats.rs"]
mod kats;
#[path = "common/orders.rs"]
mod orders;

use ephemeral_ecmh::curvegen::criteria::{self, Count};
use std::collections::HashMap;

use ephemeral_ecmh::curve::{binary127, edwards127};
use ephemeral_ecmh::curvegen::prove::{
    self, Binary127, Edwards127, Factor, Family, NoFactor, NoSieve, Rejection, Sieve, SmallL,
};
use ephemeral_ecmh::curvegen::select::{self, Certificate};
use ephemeral_ecmh::curvegen::sieve;
use ephemeral_ecmh::field::gf2_127::to_u128;
use kats::CertKat;
use orders::Orders;

/// The candidate's parameter as the tables key it: B, d or b.
trait Param {
    fn param(&self) -> u128;
}

impl Param for binary127::Curve {
    fn param(&self) -> u128 {
        to_u128(self.big_b)
    }
}

impl Param for edwards127::Curve {
    fn param(&self) -> u128 {
        self.d.value()
    }
}

/// Counts by lookup, and tallies the lookups. `quick` says from the order
/// which candidates `quick_reject` should have settled uncounted.
struct Table {
    orders: HashMap<u128, u128>,
    quick: fn(u128) -> bool,
    counted: usize,
    aborted: usize,
}

impl Table {
    fn new(o: &Orders, quick: fn(u128) -> bool) -> Self {
        Self {
            orders: o.orders.iter().copied().collect(),
            quick,
            counted: 0,
            aborted: 0,
        }
    }

    fn lookup(&self, p: u128) -> u128 {
        let n = self.orders[&p];
        assert!(!(self.quick)(n), "counted a quick rejection: #E = {n}");
        n
    }

    /// The candidates not settled by `quick_reject`: all that get counted.
    fn slow(&self) -> usize {
        self.orders.values().filter(|&&n| !(self.quick)(n)).count()
    }
}

struct Factors(HashMap<u128, u128>);

impl Factors {
    fn new(o: &Orders) -> Self {
        Self(o.factors.iter().copied().collect())
    }
}

impl<C: Param> Count<C> for Table {
    fn order(&mut self, c: &C) -> u128 {
        self.counted += 1;
        self.lookup(c.param())
    }

    /// ellsea's contract, with its Elkies primes replaced by all l < 50.
    fn order_early_abort(&mut self, c: &C, tors: u32) -> Option<u128> {
        let n = self.lookup(c.param());
        let small = [2u128, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47];
        if small
            .iter()
            .any(|&l| n.is_multiple_of(l) && !u128::from(tors).is_multiple_of(l))
        {
            self.aborted += 1;
            return None;
        }
        self.counted += 1;
        Some(n)
    }
}

impl Factor for Factors {
    fn smallest_prime_factor(&mut self, m: u128) -> Option<u128> {
        Some(self.0[&m])
    }
}

/// The KAT's own rejections, served as a sieve.
struct KatSieve(HashMap<u128, Rejection>);

impl KatSieve {
    fn new<F: Family>(k: &CertKat) -> Self
    where
        F::Curve: Param,
    {
        let live = (0..k.index).filter_map(|j| F::candidate(&k.seed, j));
        let entries = k.rejections.iter().map(|&(l, p)| (l, p.to_le_bytes()));
        if F::ENTRY_PER_INDEX {
            assert!((0..k.index).all(|j| F::candidate(&k.seed, j).is_some()));
        }
        Self(live.map(|c| c.param()).zip(entries).collect())
    }
}

impl<C: Param> Sieve<C> for KatSieve {
    fn reject(&mut self, c: &C) -> Option<Rejection> {
        self.0.get(&c.param()).copied()
    }
}

fn ls(c: &Certificate) -> Vec<u128> {
    c.rejections.iter().map(|e| e.0).collect()
}

fn kat_ls(k: &CertKat) -> Vec<u128> {
    k.rejections.iter().map(|e| e.0).collect()
}

fn kat_rejections(k: &CertKat) -> Vec<Rejection> {
    k.rejections
        .iter()
        .map(|&(l, p)| (l, p.to_le_bytes()))
        .collect()
}

fn check<F: Family>(
    o: &Orders,
    k: &CertKat,
    quick: fn(u128) -> bool,
    l_max: u32,
    verify: fn(&[u8; 32], &Certificate) -> Result<F::Curve, select::Error>,
) where
    F::Curve: Param,
    SmallL: Sieve<F::Curve>,
{
    assert_eq!(o.seed, k.seed);
    for c in (0..=k.index).filter_map(|j| F::candidate(&k.seed, j)) {
        assert_eq!(F::reject_without_count(&c), F::quick_reject(&c).is_some());
    }

    let mut t = Table::new(o, quick);
    let (cert, _) = prove::prove::<F>(&o.seed, &mut t, &mut Factors::new(o), &mut NoSieve);
    assert_eq!((cert.index, cert.r), (k.index, k.r));
    // Sage derives the witnesses as the prover does, with its own arithmetic
    assert_eq!(cert.rejections, kat_rejections(k));
    assert_eq!(t.counted, t.slow());
    verify(&o.seed, &cert).unwrap();

    // proving is deterministic
    let mut t2 = Table::new(o, quick);
    let (again, _) = prove::prove::<F>(&o.seed, &mut t2, &mut Factors::new(o), &mut NoSieve);
    assert_eq!(again.rejections, cert.rejections);

    // a sieve's rejections go into the certificate uncounted, after the
    // quick ones
    let mut t = Table::new(o, quick);
    let mut sieve = KatSieve::new::<F>(k);
    let (sieved, _) = prove::prove::<F>(&o.seed, &mut t, &mut Factors::new(o), &mut sieve);
    assert_eq!(t.counted, 1);
    assert_eq!((sieved.index, sieved.r), (k.index, k.r));
    assert_eq!(sieved.rejections, kat_rejections(k));
    verify(&o.seed, &sieved).unwrap();

    // so does the torsion sieve, whose witnesses the verifier takes
    let mut t = Table::new(o, quick);
    let mut sieve = SmallL(l_max);
    let (sieved, _) = prove::prove::<F>(&o.seed, &mut t, &mut Factors::new(o), &mut sieve);
    assert_eq!((sieved.index, sieved.r), (k.index, k.r));
    assert_eq!(ls(&sieved), kat_ls(k));
    assert!(ls(&sieved).iter().any(|&l| l % 2 == 1 && l <= l_max.into()));
    verify(&o.seed, &sieved).unwrap();
    let found = criteria::find::<F>(&o.seed, &mut Table::new(o, quick), &mut sieve);
    assert_eq!((found.index, found.r), (k.index, k.r));

    // early abort finds the same curve, counting less
    let mut t = Table::new(o, quick);
    let found = criteria::find::<F>(&o.seed, &mut t, &mut NoSieve);
    assert_eq!((found.index, found.r), (k.index, k.r));
    assert_eq!(t.counted + t.aborted, t.slow());
    assert!(t.aborted > 0);
}

fn none(_: u128) -> bool {
    false
}

#[test]
fn gf2_127_matches_sage() {
    check::<Binary127>(
        &orders::GF2_127,
        &kats::GF2_127_CERTS[0],
        none,
        sieve::GF2_127_L_MAX,
        select::verify_gf2_127,
    );
}

#[test]
fn fp127_matches_sage() {
    let eight = |n: u128| n.is_multiple_of(8);
    check::<Edwards127>(
        &orders::FP127,
        &kats::FP127_CERTS[0],
        eight,
        sieve::EDWARDS127.prove,
        select::verify_fp127,
    );
}

#[test]
fn rust_alone_finds_the_binary127_kat_curves() {
    for k in kats::GF2_127_CERTS {
        let f = prove::find_gf2_127(&k.seed);
        assert_eq!((f.index, f.r), (k.index, k.r));
    }
}

/// Without factoring, the composites Sage factored are rejected by their
/// order instead; the rest of the certificate is unchanged.
#[test]
fn gf2_127_rejects_composites_by_order() {
    let (o, k) = (&orders::GF2_127, &kats::GF2_127_CERTS[0]);
    let mut t = Table::new(o, none);
    let (cert, _) = prove::prove::<Binary127>(&o.seed, &mut t, &mut NoFactor, &mut NoSieve);
    assert_eq!((cert.index, cert.r), (k.index, k.r));
    select::verify_gf2_127(&o.seed, &cert).unwrap();
    // trial division settles factors below 2^10; the rest are orders
    let factored: Vec<u128> = o
        .factors
        .iter()
        .map(|f| f.1)
        .filter(|&l| l >= 1 << 10)
        .collect();
    let mut by_order = 0;
    for (got, kat) in ls(&cert).iter().zip(kat_ls(k)) {
        if got & 1 == 0 {
            // an order, where Sage gave the factor it found
            assert!(factored.contains(&kat), "{got} for {kat}");
            by_order += 1;
        } else {
            assert_eq!(*got, kat);
        }
    }
    assert_eq!(by_order, factored.len());
}

/// An order witness never rejects a candidate of order 2r, r prime: its
/// only admissible n is 2r, and r is prime. On a composite, the order is
/// the one n that works.
#[test]
fn gf2_127_order_witnesses() {
    let k = &kats::GF2_127_CERTS[0];
    let valid = select::gf2_127_candidate(&k.seed, k.index).unwrap();
    let h = ephemeral_ecmh::hash::Salted::new(b"test", &k.seed);
    let p = valid.hash_to_curve(&h, b"p".as_slice()).encode();
    // 2r passes all but primality; its neighbours fail n*P = O
    for n in [2 * k.r, 2 * k.r - 2, 2 * k.r + 2, 2 * k.r + 4] {
        assert!(!select::gf2_127_rejects(&valid, n, p), "{n}");
    }
    let orders: HashMap<u128, u128> = orders::GF2_127.orders.iter().copied().collect();
    let factored: Vec<u128> = orders::GF2_127
        .factors
        .iter()
        .filter(|f| f.1 >= 1 << 10)
        .map(|f| f.0)
        .collect();
    let (c, n) = (0..k.index)
        .filter_map(|j| select::gf2_127_candidate(&k.seed, j))
        .map(|c| (c, orders[&to_u128(c.big_b)]))
        .find(|(_, n)| factored.contains(&(n / 2)))
        .expect("a composite with no small factor");
    let p = c.hash_to_curve(&h, b"p".as_slice()).encode();
    assert!(select::gf2_127_rejects(&c, n, p));
    // a multiple of ord(P) other than #E is out of the Hasse interval
    for m in [n - 2, n + 2, n / 2, 2 * n] {
        assert!(!select::gf2_127_rejects(&c, m, p), "{m}");
    }
    assert!(!select::gf2_127_rejects(&c, n, [0; 16]));
}

#[test]
fn rust_alone_certifies_the_binary127_kat_curves() {
    for k in kats::GF2_127_CERTS {
        let (cert, _) = prove::certify_gf2_127(&k.seed);
        assert_eq!((cert.index, cert.r), (k.index, k.r));
        select::verify_gf2_127(&k.seed, &cert).unwrap();
    }
}
