//! A finite RIBLT simulator with an ECMH checksum, using the mapping from
//! Yang, Gilad, and Alizadeh, "Practical Rateless Set Reconciliation"
//! (SIGCOMM 2024; yangl1996/riblt).
//!
//! Each cell stores the XOR of item bytes, a signed count, and the signed sum
//! of their hashed curve points. Item bytes remain the recovery payload;
//! curve points provide a separate checksum. This simulator does not implement
//! the point-identifier protocol or a wire format.
//!
//! The caller supplies the salt used by separate item-hashing and mapping tags.
//! Each item's mapping is drawn from xoshiro256++ seeded by its salted map
//! digest. Changing the salt changes both the hash namespace and the mapping
//! inputs. [`Mcg64`], upstream's generator seeded by 64 bits of the digest,
//! and [`ChaCha8`], riblt-ecmh's, remain available through
//! [`Riblt::with_prng`] for comparison.
//!
//! One `encode` and one `peel` are one round of reconciliation. Over repeated
//! rounds, as docs/workload.md (Cost in repeated reconciliation) models
//! them, the parts are counted separately.
//! An item maps into k(m) = 2(H_{m+1} - 1) of the first m cells on average.
//! Its [`Riblt::addends`] are valid for as long as the curve, salt and mapping
//! are, so a round that rebuilds its cells from n items through
//! [`Riblt::apply`] costs n k(m) additions and hashes only the new items. A
//! party that keeps its cells between rounds instead performs one hash, one
//! preparation and k(m) additions per new item. `peel` is the decoding either
//! way.

use crate::ecmh::TAG_ITEM;
use crate::group::{Group, HashToCurve, Negate};
use crate::hash::Salted;
use core::marker::PhantomData;
use rand_chacha::ChaCha8Rng;
use rand_chacha::rand_core::{RngCore, SeedableRng};
use rand_xoshiro::Xoshiro256PlusPlus;

pub const TAG_MAP: &[u8] = b"ephemeral-ecmh/riblt-map";

/// The uniform 64-bit samples from which a [`Mapping`] draws its gaps,
/// seeded by an item's salted map digest.
pub trait Prng: Clone {
    fn seed(digest: &[u8; 32]) -> Self;
    fn next_u64(&mut self) -> u64;
}

/// The multiplicative congruential generator of yangl1996/riblt, seeded by
/// the digest's first 8 bytes. Two items whose 64-bit seeds are equal have
/// equal mappings, so under a known salt such a pair is found after about
/// 2^32 hashes.
#[derive(Clone, Copy, Debug)]
pub struct Mcg64(pub u64);

impl Prng for Mcg64 {
    fn seed(digest: &[u8; 32]) -> Self {
        Self(u64::from_le_bytes(digest[..8].try_into().unwrap()))
    }

    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(0xda942042e4dd58b5);
        self.0
    }
}

/// xoshiro256++ whose state is the 32-byte digest itself, so equal states
/// require equal digests, with one exception in 2^256: rand_xoshiro reseeds
/// the all-zero digest from SplitMix64, so it shares a state with the one
/// digest that encodes that state. The generator runs through its 2^256 - 1
/// nonzero states in one cycle, so distinct states give distinct sample
/// streams. From a state of few set bits its first samples are poorly mixed,
/// as for every xoshiro; the salted digest is modelled as uniform, under
/// which such a state is no likelier than any other. Cryptographic strength
/// does not enter: the digest is a salted hash and the mapping is public.
#[derive(Clone, Debug)]
pub struct Xoshiro256pp(Xoshiro256PlusPlus);

impl Prng for Xoshiro256pp {
    fn seed(digest: &[u8; 32]) -> Self {
        Self(Xoshiro256PlusPlus::from_seed(*digest))
    }

    fn next_u64(&mut self) -> u64 {
        self.0.next_u64()
    }
}

/// ChaCha8 keyed by the whole 32-byte digest: the generator of riblt-ecmh,
/// as Go's ChaCha8Rand, whose output differs from rand_chacha's.
///
/// A ChaCha8 block is 8 rounds of 32-bit additions, XORs and rotations and
/// yields eight samples; rand_chacha computes four blocks per refill, so one
/// refill covers an item's k(m) + 1 samples, 19 on average for m = 12150. A
/// SHA-256 counter stream yields four samples per 64-round compression, with
/// hardware support on aarch64 through the SHA-2 extension and on x86-64
/// through SHA-NI. benches/riblt.rs (riblt.mapping) times all of them.
#[derive(Clone, Debug)]
pub struct ChaCha8(ChaCha8Rng);

impl Prng for ChaCha8 {
    fn seed(digest: &[u8; 32]) -> Self {
        Self(ChaCha8Rng::from_seed(*digest))
    }

    fn next_u64(&mut self) -> u64 {
        self.0.next_u64()
    }
}

/// The RIBLT `randomMapping` recurrence, starting at cell 0, over the
/// samples of `P`.
/// Its idealized inclusion probability at index i is 1 / (1 + i/2), giving
/// an expected 2(H_{m+1} - 1) assigned cells among the first m cells.
#[derive(Clone, Copy, Debug)]
pub struct Mapping<P = Xoshiro256pp> {
    prng: P,
    index: u64,
}

impl<P: Prng> Mapping<P> {
    pub fn new(prng: P) -> Self {
        Self { prng, index: 0 }
    }

    pub fn index(&self) -> u64 {
        self.index
    }

    /// The next index. For the 2^10 values of r nearest u64::MAX, r + 1 is
    /// 2^64 in f64 and the gap is zero, so the item enters the same cell
    /// again: upstream's recurrence does the same, and every peer that
    /// computes it in IEEE 754 doubles does alike, so cell sums stay
    /// consistent. It happens with probability
    /// about 2^-54 per step.
    pub fn advance(&mut self) -> u64 {
        let r = self.prng.next_u64();
        let gap = ((self.index as f64 + 1.5)
            * ((1u64 << 32) as f64 / (r as f64 + 1.0).sqrt() - 1.0))
            .ceil();
        // Both the float-to-integer conversion and the index addition saturate.
        self.index = self.index.saturating_add(gap as u64);
        self.index
    }

    /// The indices below m, from the current one.
    pub fn below(mut self, m: usize) -> impl Iterator<Item = usize> {
        core::iter::from_fn(move || {
            let i = self.index;
            (i < m as u64).then(|| {
                self.advance();
                i as usize
            })
        })
    }
}

/// A simulator cell containing an L-byte XOR payload, signed count, and checksum.
#[derive(Clone, Copy, Debug)]
pub struct Cell<G: Group, const L: usize> {
    pub key: [u8; L],
    pub count: i64,
    pub sum: G::Point,
}

/// Options for testing candidate cells whose signed count is +1 or -1.
#[derive(Clone, Copy, Debug)]
pub struct Peel {
    /// Reject a candidate when its recovered key's mapping excludes the cell,
    /// before evaluating that candidate's hash-to-curve checksum.
    pub prefilter: bool,
    /// Hash each wave's candidate keys together (`HashToCurve::hash_batch`).
    pub batch: bool,
}

/// A decoded difference: the items only in A, and those only in B.
pub type Difference<const L: usize> = (Vec<[u8; L]>, Vec<[u8; L]>);

/// A finite RIBLT whose item mappings draw from `P`.
pub struct Riblt<G: Group, P = Xoshiro256pp> {
    pub group: G,
    pub item: Salted,
    pub map: Salted,
    prng: PhantomData<P>,
}

fn xor<const L: usize>(a: &mut [u8; L], b: &[u8; L]) {
    a.iter_mut().zip(b).for_each(|(a, b)| *a ^= b);
}

impl<G: HashToCurve + Negate> Riblt<G> {
    pub fn new(group: G, salt: &[u8; 32]) -> Self {
        Self::with_prng(group, salt)
    }
}

impl<G: HashToCurve + Negate, P: Prng> Riblt<G, P> {
    pub fn with_prng(group: G, salt: &[u8; 32]) -> Self {
        Self {
            group,
            item: Salted::new(TAG_ITEM, salt),
            map: Salted::new(TAG_MAP, salt),
            prng: PhantomData,
        }
    }

    pub fn mapping(&self, item: &[u8]) -> Mapping<P> {
        Mapping::new(P::seed(&self.map.digest(item, 0)))
    }

    pub fn cells<const L: usize>(&self, m: usize) -> Vec<Cell<G, L>> {
        let empty = Cell {
            key: [0; L],
            count: 0,
            sum: self.group.identity(),
        };
        vec![empty; m]
    }

    /// Prepared addends for items, hashed as a batch.
    pub fn addends<const L: usize>(&self, items: &[[u8; L]]) -> Vec<G::Addend> {
        let refs: Vec<&[u8]> = items.iter().map(|x| x.as_slice()).collect();
        let g = &self.group;
        g.prepare_batch(&g.hash_batch(&self.item, &refs))
    }

    /// Apply item `x` and its prepared checksum addend `a` to all assigned cells.
    /// The caller must supply the matching addend and a sign of +1 or -1.
    pub fn apply<const L: usize>(
        &self,
        cells: &mut [Cell<G, L>],
        x: &[u8; L],
        a: &G::Addend,
        sign: i64,
    ) {
        let a = if sign > 0 {
            *a
        } else {
            self.group.neg_addend(a)
        };
        for i in self.mapping(x).below(cells.len()) {
            let c = &mut cells[i];
            xor(&mut c.key, x);
            c.count += sign;
            c.sum = self.group.add(&c.sum, &a);
        }
    }

    /// Adds (sign 1) or removes (sign -1) items in every cell they map to.
    pub fn encode<const L: usize>(&self, cells: &mut [Cell<G, L>], items: &[[u8; L]], sign: i64) {
        for (x, a) in items.iter().zip(self.addends(items)) {
            self.apply(cells, x, &a, sign);
        }
    }

    /// Peel cells representing A minus B into the item lists A \ B and B \ A.
    /// Return `None` if peeling stops with a nonempty count or checksum.
    /// A candidate passes the purity test when its count is +1 or -1 and its
    /// checksum equals the hashed XOR payload with that sign. The final emptiness
    /// test requires count zero and the group identity; it does not inspect the
    /// XOR payload. Cells are updated in place, including when peeling stalls.
    pub fn peel<const L: usize>(
        &self,
        cells: &mut [Cell<G, L>],
        how: Peel,
    ) -> Option<Difference<L>> {
        let g = &self.group;
        let (mut plus, mut minus) = (vec![], vec![]);
        let mut wave: Vec<usize> = (0..cells.len()).collect();
        while !wave.is_empty() {
            let cand: Vec<(usize, [u8; L])> = wave
                .iter()
                .filter(|&&i| cells[i].count.abs() == 1)
                .filter(|&&i| {
                    !how.prefilter || self.mapping(&cells[i].key).below(i + 1).any(|j| j == i)
                })
                .map(|&i| (i, cells[i].key))
                .collect();
            // Share hashes, but retain each cell: equal XOR payloads can have
            // different checksums.
            let mut keys: Vec<[u8; L]> = cand.iter().map(|c| c.1).collect();
            keys.sort_unstable();
            keys.dedup();
            let adds = if how.batch {
                self.addends(&keys)
            } else {
                keys.iter()
                    .map(|k| g.prepare(&g.hash(&self.item, k)))
                    .collect()
            };
            let mut next = vec![];
            for (i, key) in cand {
                // Recheck the current count and payload: earlier peels in this
                // wave may have changed the cell. Affected cells are also queued
                // for the next wave.
                let c = &cells[i];
                if c.count.abs() != 1 || c.key != key {
                    continue;
                }
                let a = adds[keys.binary_search(&key).unwrap()];
                let sign = c.count;
                let item = if sign > 0 { a } else { g.neg_addend(&a) };
                if !g.equals_addend(&c.sum, &item) {
                    continue;
                }
                if sign > 0 { &mut plus } else { &mut minus }.push(key);
                self.apply(cells, &key, &a, -sign);
                next.extend(self.mapping(&key).below(cells.len()));
            }
            next.sort_unstable();
            next.dedup();
            wave = next;
        }
        cells
            .iter()
            .all(|c| c.count == 0 && g.is_identity(&c.sum))
            .then_some((plus, minus))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curve::binary::{lambda, unscaled, wcodec};
    use crate::curve::weier::jacobian;
    use crate::curve::{
        binary109, binary122, binary127, edwards, edwards61x2, edwards107, edwards127, twisted,
        twisted_goldilocks2, twisted61x2, twisted64x2, twisted128, weier, weier61x2, weier107,
        weier127,
    };
    use crate::field::{Packed, fp61x2, fp64x2, fp107, fp127, fp128, goldilocks2};
    use core::convert::identity;
    use proptest::prelude::*;

    /// Distinct map digests, as uniform as the salted hash makes them: from
    /// a state with 192 zero bits, xoshiro256++'s second sample repeats its
    /// first.
    fn digest(s: u64) -> [u8; 32] {
        Salted::new(b"riblt-test-map", &[0; 32]).digest(&s.to_le_bytes(), 0)
    }

    /// The mean number of indices below m, against k(m) = 2(H_{m+1} - 1).
    fn density<P: Prng>() {
        for m in [20usize, 150, 1350] {
            let k = 2.0 * ((1..=m + 1).map(|i| 1.0 / i as f64).sum::<f64>() - 1.0);
            let n = 20_000;
            let total: usize = (0..n)
                .map(|s| Mapping::new(P::seed(&digest(s))).below(m).count())
                .sum();
            let mean = total as f64 / n as f64;
            assert!((mean / k - 1.0).abs() < 0.02, "m {m}: mean {mean}, k {k}");
        }
    }

    #[test]
    fn mapping_density_matches_riblt() {
        density::<Mcg64>();
        density::<ChaCha8>();
        density::<Xoshiro256pp>();
    }

    /// How many of n mappings include each index below 101.
    fn inclusions<P: Prng>(n: u64) -> [u64; 101] {
        let mut hits = [0; 101];
        for s in 0..n {
            for i in Mapping::new(P::seed(&digest(s))).below(hits.len()) {
                hits[i] += 1;
            }
        }
        hits
    }

    /// Every generator includes each index with upstream's frequency, within
    /// five standard deviations of the difference of two binomial counts.
    /// The frequencies are those of the recurrence, not the idealized
    /// 1 / (1 + i/2): index 1 follows index 0 with probability 0.64, not 2/3.
    #[test]
    fn mapping_inclusion_agrees_between_generators() {
        let n = 20_000;
        let a = inclusions::<Mcg64>(n);
        for b in [inclusions::<ChaCha8>(n), inclusions::<Xoshiro256pp>(n)] {
            for i in [0, 1, 2, 3, 10, 30, 100] {
                let p = (a[i] + b[i]) as f64 / (2 * n) as f64;
                let sd = (2.0 * n as f64 * p * (1.0 - p)).sqrt();
                let diff = a[i].abs_diff(b[i]) as f64;
                assert!(diff <= 5.0 * sd.max(1.0), "i {i}: {} vs {}", a[i], b[i]);
            }
            assert!((b[1] as f64 / n as f64 - 0.64).abs() < 0.02, "{}", b[1]);
        }
    }

    fn indices<P: Prng>(d: &[u8; 32]) -> Vec<usize> {
        Mapping::new(P::seed(d)).below(1 << 20).collect()
    }

    #[test]
    fn mapping_starts_at_0_and_increases() {
        for v in [
            indices::<Mcg64>(&digest(12345)),
            indices::<ChaCha8>(&digest(12345)),
            indices::<Xoshiro256pp>(&digest(12345)),
        ] {
            assert_eq!(v[0], 0);
            assert!(v.windows(2).all(|w| w[0] < w[1]));
        }
    }

    #[test]
    fn wide_mappings_depend_on_the_whole_digest() {
        let (a, mut b) = (digest(7), digest(7));
        b[31] ^= 1;
        assert_eq!(indices::<Mcg64>(&a), indices::<Mcg64>(&b));
        assert_ne!(indices::<ChaCha8>(&a), indices::<ChaCha8>(&b));
        assert_ne!(indices::<Xoshiro256pp>(&a), indices::<Xoshiro256pp>(&b));
    }

    #[test]
    fn mapping_repeats_a_cell_where_upstream_does() {
        // the state whose next r is u64::MAX: r + 1 rounds to 2^64
        let prng: u64 = 0x747c_72fc_ab15_2a63;
        assert_eq!(prng.wrapping_mul(0xda942042e4dd58b5), u64::MAX);
        let mut m = Mapping {
            prng: Mcg64(prng),
            index: 7,
        };
        assert_eq!(m.advance(), 7);
        assert!(m.advance() > 7);
    }

    fn item(tag: u32, i: u32) -> [u8; 36] {
        let mut x = [0u8; 36];
        x[..32].copy_from_slice(&Salted::new(b"riblt-test", &[0; 32]).digest(&[], tag));
        x[32..].copy_from_slice(&i.to_le_bytes());
        x
    }

    /// A = common + only_a, B = common + only_b, in m cells.
    fn reconcile<G: HashToCurve + Negate, P: Prng>(
        g: G,
        salt: &[u8; 32],
        common: u32,
        only_a: u32,
        only_b: u32,
        m: usize,
        how: Peel,
    ) -> Option<Difference<36>> {
        let r = Riblt::<G, P>::with_prng(g, salt);
        let c: Vec<_> = (0..common).map(|i| item(0, i)).collect();
        let a: Vec<_> = (0..only_a).map(|i| item(1, i)).collect();
        let b: Vec<_> = (0..only_b).map(|i| item(2, i)).collect();
        let mut cells = r.cells(m);
        r.encode(&mut cells, &[c.clone(), a].concat(), 1);
        r.encode(&mut cells, &[c, b].concat(), -1);
        r.peel(&mut cells, how)
    }

    /// Check returned item lists against the generated difference; report success.
    fn check<G: HashToCurve + Negate, P: Prng>(
        g: G,
        salt: [u8; 32],
        common: u32,
        da: u32,
        db: u32,
        m: usize,
    ) -> Result<bool, TestCaseError> {
        let expect = |tag, n| {
            let mut v: Vec<_> = (0..n).map(|i| item(tag, i)).collect();
            v.sort();
            v
        };
        let mut results = vec![];
        for prefilter in [false, true] {
            for batch in [false, true] {
                let got = reconcile::<G, P>(g, &salt, common, da, db, m, Peel { prefilter, batch });
                if let Some((mut a, mut b)) = got.clone() {
                    a.sort();
                    b.sort();
                    prop_assert_eq!(&a, &expect(1, da));
                    prop_assert_eq!(&b, &expect(2, db));
                }
                results.push(got.is_some());
            }
        }
        prop_assert!(results.iter().all(|&r| r == results[0]));
        Ok(results[0])
    }

    /// The curve with the smallest admissible constant among those with a
    /// set bit in the packing's upper half: for the quadratic extensions
    /// that keeps the constant out of the base field, where every element
    /// is a square and so admits neither Edwards model.
    fn first<F: Packed, C>(new: impl Fn(F) -> Option<C>) -> C {
        (2u128..1 << 16)
            .find_map(|k| new(F::unpack(k | 1 << (F::BITS / 2))?))
            .expect("no admissible constant among the first 2^16 searched")
    }

    /// Fixed curves of every family, each under every accumulator and codec
    /// the benches use, passed one at a time to a generic `$f(group)`.
    macro_rules! every_group {
        ($f:ident) => {{
            let bin = binary127::Curve::new(crate::field::gf2_127::from_u128(
                0x1234_5678_9abc_def0_1234_5678_9abc_def1,
            ));
            let bin109 = binary109::Curve::new(crate::field::gf2_109::from_u128(
                0x1234_5678_9abc_def0_1234_5678_9abc_def1 & crate::field::gf2_109::MASK109,
            ));
            let (dense, gls) = binary122_curves();
            $f(bin);
            $f(lambda::Curve(bin));
            $f(wcodec::Curve::new(bin));
            $f(unscaled::Curve::new(bin));
            $f(dense);
            $f(gls);
            $f(lambda::Curve(dense));
            $f(lambda::Curve(gls));
            $f(wcodec::Curve::new(dense));
            $f(wcodec::Curve::new(gls));
            $f(unscaled::Curve::new(dense));
            $f(unscaled::Curve::new(gls));
            $f(bin109);
            $f(lambda::Curve(bin109));
            $f(wcodec::Curve::new(bin109));
            $f(unscaled::Curve::new(bin109));

            $f(first(edwards::Curve::<fp127::Fp>::new));
            $f(first(edwards::Curve::<fp107::Fp>::new));
            $f(first(edwards::Curve::<fp61x2::Fq>::new));
            $f(first(twisted::Curve::<fp128::Fp>::new));
            $f(first(twisted::Curve::<fp61x2::Fq>::new));
            $f(first(twisted::Curve::<fp64x2::Fq>::new));
            $f(first(twisted::Curve::<goldilocks2::Fq>::new));
            let we127 = first(|b: fp127::Fp| weier::Curve::new(b).and_then(weier::OddCurve::new));
            let we107 = first(|b: fp107::Fp| weier::Curve::new(b).and_then(weier::OddCurve::new));
            let we61 = first(|b: fp61x2::Fq| weier::Curve::new(b).and_then(weier::OddCurve::new));
            $f(we127);
            $f(we107);
            $f(we61);
            $f(jacobian::Curve(we127));
            $f(jacobian::Curve(we107));
            $f(jacobian::Curve(we61));
        }};
    }

    fn peels_50_in_100<G: HashToCurve + Negate>(g: G) {
        assert!(check::<G, Mcg64>(g, [3; 32], 100, 30, 20, 100).unwrap());
        assert!(check::<G, ChaCha8>(g, [3; 32], 100, 30, 20, 100).unwrap());
        assert!(check::<G, Xoshiro256pp>(g, [3; 32], 100, 30, 20, 100).unwrap());
    }

    /// Recover 50 differing items from 100 cells with fixed curves and salt.
    #[test]
    fn peels_50_differences_in_100_cells() {
        every_group!(peels_50_in_100);
    }

    fn binary122_curves() -> (binary122::Dense, binary122::Gls) {
        let dense = binary122::Dense::new(crate::field::gf2_122::from_u128(
            0x0234_5678_9abc_def0_1234_5678_9abc_def1,
        ));
        let beta = crate::field::gf2_122::gf2_61::from_u64(0x0123_4567_89ab_cdef);
        (dense, binary122::Gls::new(beta.xsquare(2)))
    }

    /// Check selected inconsistent count/checksum cases and valid singletons.
    fn rejects_forgeries<G: HashToCurve + Negate>(g: G) {
        let r = Riblt::new(g, &[5; 32]);
        let (x, y) = (item(1, 0), item(2, 0));
        let hx = g.prepare(&g.hash(&r.item, &x));
        let hy = g.prepare(&g.hash(&r.item, &y));
        let id = g.identity();
        let cell = |key, count, sum| Cell { key, count, sum };
        let peel = |cells: &[Cell<G, 36>]| {
            [false, true].map(|prefilter| {
                let how = Peel {
                    prefilter,
                    batch: prefilter,
                };
                r.peel(&mut cells.to_vec(), how)
            })
        };
        // the genuine cell, either sign
        for (count, sum, want) in [
            (1, g.add(&id, &hx), (vec![x], vec![])),
            (-1, g.add(&id, &g.neg_addend(&hx)), (vec![], vec![x])),
        ] {
            for got in peel(&[cell(x, count, sum)]) {
                assert_eq!(got, Some(want.clone()));
            }
        }
        for forged in [
            // count +1, no checksum
            [cell(x, 1, id), cell([0; 36], 0, id)],
            // count +1, checksum of -x
            [
                cell(x, 1, g.add(&id, &g.neg_addend(&hx))),
                cell([0; 36], 0, id),
            ],
            // checksum of x, count 2
            [cell(x, 2, g.add(&id, &hx)), cell([0; 36], 0, id)],
            // count 0 and key 0, checksum x - y: not empty
            [
                cell([0; 36], 0, id),
                cell([0; 36], 0, g.add(&g.add(&id, &hx), &g.neg_addend(&hy))),
            ],
        ] {
            for got in peel(&forged) {
                assert_eq!(got, None);
            }
        }
    }

    #[test]
    fn peel_rejects_forged_cells() {
        every_group!(rejects_forgeries);
    }

    /// The XOR of the first 8 bytes of each item's salted digest, the checksum
    /// yangl1996/riblt#3 is about, in the simulator's cells.
    #[derive(Clone, Copy)]
    struct Xor64;

    impl Group for Xor64 {
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

    impl crate::group::Accumulate for Xor64 {
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

    impl Negate for Xor64 {
        fn neg(&self, a: &u64) -> u64 {
            *a
        }
        fn neg_addend(&self, a: &u64) -> u64 {
            *a
        }
    }

    impl HashToCurve for Xor64 {
        fn hash(&self, h: &Salted, msg: &[u8]) -> u64 {
            u64::from_le_bytes(h.digest(msg, 0)[..8].try_into().unwrap())
        }
    }

    /// The indices of a subset of `rows` whose XOR is `target`, if one exists,
    /// by Gaussian elimination over F_2: each basis vector carries the mask of
    /// the rows it combines. `v ^ b < v` exactly when `b`'s leading bit is set
    /// in `v`, so with the basis in decreasing order one pass reduces fully.
    fn dependent_subset(rows: &[u128], target: u128) -> Option<Vec<usize>> {
        let mut basis: Vec<(u128, u128)> = vec![];
        let reduce = |basis: &[(u128, u128)], mut v: u128, mut m: u128| {
            for &(b, bm) in basis {
                if v ^ b < v {
                    v ^= b;
                    m ^= bm;
                }
            }
            (v, m)
        };
        for (i, &row) in rows.iter().enumerate() {
            let (v, m) = reduce(&basis, row, 1 << i);
            if v != 0 {
                basis.push((v, m));
                basis.sort_unstable_by_key(|a| core::cmp::Reverse(a.0));
            }
        }
        let (v, m) = reduce(&basis, target, 0);
        (v == 0).then(|| (0..rows.len()).filter(|i| m >> i & 1 == 1).collect())
    }

    /// Two sets whose XOR-checksum difference peels as one item in neither.
    /// The 120 candidates `item(3, i)` share 32 bytes, so a set's XOR payload
    /// is determined by its parity and the XOR of its indices; with the
    /// 64-bit hash that is 97 linear constraints over F_2, and the solution is
    /// a dependent subset of odd size whose payload and checksum are those of
    /// the absent `e = item(3, 0)`. Split so that A holds one more of them
    /// than B, one cell holding A - B reads as the pure cell of `e`.
    fn forged_sets() -> (Vec<[u8; 36]>, Vec<[u8; 36]>, [u8; 36]) {
        let r = Riblt::new(Xor64, &FORGERY_SALT);
        let row = |i: u32| Xor64.hash(&r.item, &item(3, i)) as u128 | (i as u128) << 64 | 1 << 96;
        // indices spread over 32 bits, none of them 0
        let index = |k: u32| k.wrapping_mul(0x9e37_79b1);
        let rows: Vec<u128> = (1..=120).map(|k| row(index(k))).collect();
        let s = dependent_subset(&rows, row(0)).expect("97 constraints over 120 rows");
        assert_eq!(s.len() % 2, 1);
        let side = |p| {
            s.iter()
                .skip(p)
                .step_by(2)
                .map(|&j| item(3, index(j as u32 + 1)))
                .collect()
        };
        (side(0), side(1), item(3, 0))
    }

    const FORGERY_SALT: [u8; 32] = [7; 32];

    /// A - B in one cell, with common items on both sides.
    fn forged_cell<G: HashToCurve + Negate>(g: G, m: usize) -> (Riblt<G>, Vec<Cell<G, 36>>) {
        let (a, b, _) = forged_sets();
        let common: Vec<_> = (0..10).map(|i| item(0, i)).collect();
        let r = Riblt::new(g, &FORGERY_SALT);
        let mut cells = r.cells(m);
        r.encode(&mut cells, &[common.clone(), a].concat(), 1);
        r.encode(&mut cells, &[common, b].concat(), -1);
        (r, cells)
    }

    const PEELS: [Peel; 2] = [
        Peel {
            prefilter: false,
            batch: false,
        },
        Peel {
            prefilter: true,
            batch: true,
        },
    ];

    /// yangl1996/riblt#3 through the simulator: the XOR checksum is linear, so
    /// the forged sets decode to a difference of one item that neither set
    /// holds, and the mapping prefilter passes it since index 0 begins every
    /// mapping.
    #[test]
    fn xor_checksum_decodes_a_forged_difference() {
        let (a, b, e) = forged_sets();
        assert!(!a.contains(&e) && !b.contains(&e));
        for how in PEELS {
            let (r, mut cells) = forged_cell(Xor64, 1);
            assert_eq!(r.peel(&mut cells, how), Some((vec![e], vec![])));
        }
    }

    /// The same sets under a curve checksum: one cell is not pure, so peeling
    /// reports failure rather than an item, and enough cells decode A - B.
    fn rejects_the_forgery<G: HashToCurve + Negate>(g: G) {
        let (mut a, mut b, _) = forged_sets();
        a.sort();
        b.sort();
        for how in PEELS {
            let (r, mut cells) = forged_cell(g, 1);
            assert_eq!(r.peel(&mut cells, how), None);
            let (r, mut cells) = forged_cell(g, 100);
            let (mut got_a, mut got_b) = r.peel(&mut cells, how).expect("100 cells decode");
            got_a.sort();
            got_b.sort();
            assert_eq!((got_a, got_b), (a.clone(), b.clone()));
        }
    }

    #[test]
    fn curve_checksums_reject_the_forged_difference() {
        every_group!(rejects_the_forgery);
    }

    /// Reconciliation over random curves, salts and set sizes, for a
    /// representative subset of the (family, representation) pairs that the
    /// fixed-curve tests cover: `name: curve strategy => wrapper`.
    macro_rules! peels {
        ($($name:ident: $curve:expr => $wrap:expr;)*) => {
            proptest! {
                #![proptest_config(ProptestConfig::with_cases(16))]
                $(
                #[test]
                fn $name(c in $curve, salt in any::<[u8; 32]>(), common in 0u32..40, da in 0u32..20, db in 0u32..20, extra in 0usize..40) {
                    check::<_, Xoshiro256pp>($wrap(c), salt, common, da, db, (da + db) as usize + extra)?;
                }
                )*
            }
        };
    }

    peels! {
        binary127_peels: binary127::tests::curve() => identity;
        binary127_lambda_peels: binary127::tests::curve() => lambda::Curve;
        binary127_w_peels: binary127::tests::curve() => wcodec::Curve::new;
        binary127_u_peels: binary127::tests::curve() => unscaled::Curve::new;
        binary122_peels: binary122::tests::dense::curve() => identity;
        binary122_gls_lambda_peels: binary122::tests::gls::curve() => lambda::Curve;
        binary122_u_peels: binary122::tests::dense::curve() => unscaled::Curve::new;
        binary109_peels: binary109::tests::curve() => identity;
        binary109_u_peels: binary109::tests::curve() => unscaled::Curve::new;
        edwards127_peels: edwards127::tests::curve() => identity;
        edwards107_peels: edwards107::tests::curve() => identity;
        edwards61x2_peels: edwards61x2::tests::curve() => identity;
        twisted128_peels: twisted128::tests::curve() => identity;
        twisted61x2_peels: twisted61x2::tests::curve() => identity;
        twisted64x2_peels: twisted64x2::tests::curve() => identity;
        twisted_goldilocks2_peels: twisted_goldilocks2::tests::curve() => identity;
        weier127_peels: weier127::tests::odd_curve() => identity;
        weier127_jacobian_peels: weier127::tests::odd_curve() => jacobian::Curve;
        weier107_peels: weier107::tests::odd_curve() => identity;
        weier61x2_peels: weier61x2::tests::odd_curve() => identity;
        weier61x2_jacobian_peels: weier61x2::tests::odd_curve() => jacobian::Curve;
    }
}
