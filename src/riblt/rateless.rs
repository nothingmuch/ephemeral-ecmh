//! A rateless RIBLT encoder and decoder, after yangl1996/riblt's encoder.go
//! and decoder.go.
//!
//! The encoder produces the unbounded sequence of coded symbols of its set A;
//! the decoder holds a set B, subtracts it from each received symbol, and
//! peels the difference out of the symbols received so far. Coded symbols are
//! the finite simulator's [`Cell`]s, and an item enters the symbols that its
//! [`Mapping`] lists, so symbol i equals cell i of a finite table of any size
//! m > i. Item hashes, mapping seeds and prepared addends are those of the
//! [`Riblt`] that the encoder and decoder borrow.
//!
//! Upstream's algorithm is retained: each set is a window of items with a
//! priority queue keyed by the next symbol each maps to; a received symbol
//! is decodable when its count is +1 or -1 and its checksum is the payload's
//! hashed point with that sign, or when its count is 0 and its checksum the
//! identity; a recovered item is removed from the symbols received so far and
//! joins a window that removes it from later ones. It differs in four ways:
//!
//! - Checksum and mapping hashes are salted and separately tagged, as in the
//!   finite simulator, where upstream derives both from one unsalted hash.
//! - The decodable queue keeps the addend that the purity test prepared, and
//!   recovery adds it; the ristretto255 version of riblt-ecmh hashes the item
//!   to the curve again when it recovers it.
//! - A window's addends carry its sign, so applying it is one addition per
//!   item, with no subtraction.
//! - Insertion hashes and prepares items as a batch.

use super::{Cell, Difference, Mapping, Prng, Riblt, xor};
use crate::group::{Accumulate, HashToCurve, Negate};
use core::cmp::Reverse;
use std::collections::BinaryHeap;

fn update<G: Accumulate, const L: usize>(
    g: &G,
    c: &mut Cell<G, L>,
    x: &[u8; L],
    a: &G::Addend,
    sign: i64,
) {
    xor(&mut c.key, x);
    c.count += sign;
    c.sum = g.add(&c.sum, a);
}

/// Items whose addends carry the window's sign, and a queue of the next coded
/// symbol each maps to: upstream's codingWindow.
struct Window<G: Accumulate, P, const L: usize> {
    sign: i64,
    items: Vec<([u8; L], G::Addend)>,
    maps: Vec<Mapping<P>>,
    /// (next symbol, item), least first.
    queue: BinaryHeap<Reverse<(u64, usize)>>,
    /// The index of the next symbol to apply the window to.
    next: u64,
}

impl<G: Accumulate, P: Prng, const L: usize> Window<G, P, L> {
    fn new(sign: i64) -> Self {
        Self {
            sign,
            items: vec![],
            maps: vec![],
            queue: BinaryHeap::new(),
            next: 0,
        }
    }

    /// Add item `x`, whose mapping `m` is at or after the next symbol.
    fn push(&mut self, x: [u8; L], a: G::Addend, m: Mapping<P>) {
        debug_assert!(m.index() >= self.next);
        self.queue.push(Reverse((m.index(), self.items.len())));
        self.items.push((x, a));
        self.maps.push(m);
    }

    /// Apply the items mapped to the next symbol to c and advance their mappings.
    /// A zero gap (see Mapping::advance) requeues an item at the same symbol,
    /// so it is applied again, as upstream does.
    fn apply(&mut self, g: &G, c: &mut Cell<G, L>) {
        while let Some(mut top) = self.queue.peek_mut() {
            let Reverse((i, s)) = *top;
            if i != self.next {
                break;
            }
            let (x, a) = &self.items[s];
            update(g, c, x, a, self.sign);
            *top = Reverse((self.maps[s].advance(), s));
        }
        self.next += 1;
    }
}

/// The coded symbols of a set A.
pub struct Encoder<'r, G: HashToCurve + Negate, P: Prng, const L: usize> {
    riblt: &'r Riblt<G, P>,
    window: Window<G, P, L>,
}

impl<'r, G: HashToCurve + Negate, P: Prng, const L: usize> Encoder<'r, G, P, L> {
    pub fn new(riblt: &'r Riblt<G, P>) -> Self {
        Self {
            riblt,
            window: Window::new(1),
        }
    }

    /// Add items to A. The coded symbols already produced would not include
    /// them, so this panics after the first symbol.
    pub fn insert(&mut self, items: &[[u8; L]]) {
        assert_eq!(self.window.next, 0, "insertion after a coded symbol");
        let r = self.riblt;
        for (x, a) in items.iter().zip(r.addends(items)) {
            self.window.push(*x, a, r.mapping(x));
        }
    }

    /// The next coded symbol.
    pub fn produce(&mut self) -> Cell<G, L> {
        let g = &self.riblt.group;
        let mut c = Cell {
            key: [0; L],
            count: 0,
            sum: g.identity(),
        };
        self.window.apply(g, &mut c);
        c
    }
}

/// The difference between a remote set A, received as coded symbols, and the
/// local set B.
pub struct Decoder<'r, G: HashToCurve + Negate, P: Prng, const L: usize> {
    riblt: &'r Riblt<G, P>,
    /// The received symbols of A, less B and the items recovered so far.
    cells: Vec<Cell<G, L>>,
    /// B, subtracted from every received symbol.
    local_set: Window<G, P, L>,
    /// Recovered items of A \ B, subtracted.
    remote: Window<G, P, L>,
    /// Recovered items of B \ A, added back.
    local: Window<G, P, L>,
    /// Symbols that tested decodable, with the payload's unsigned addend for
    /// counts of +1 and -1, in the order found.
    decodable: Vec<(usize, Option<G::Addend>)>,
    decoded: usize,
}

/// The prepared addend of `c`'s payload, if `c` has count +1 or -1 and a
/// checksum equal to that addend with the count's sign.
fn pure<G: HashToCurve + Negate, P: Prng, const L: usize>(
    r: &Riblt<G, P>,
    c: &Cell<G, L>,
) -> Option<G::Addend> {
    let g = &r.group;
    let a = g.prepare(&g.hash(&r.item, &c.key));
    let signed = if c.count > 0 { a } else { g.neg_addend(&a) };
    g.equals_addend(&c.sum, &signed).then_some(a)
}

impl<'r, G: HashToCurve + Negate, P: Prng, const L: usize> Decoder<'r, G, P, L> {
    pub fn new(riblt: &'r Riblt<G, P>) -> Self {
        Self {
            riblt,
            cells: vec![],
            local_set: Window::new(-1),
            remote: Window::new(-1),
            local: Window::new(1),
            decodable: vec![],
            decoded: 0,
        }
    }

    /// Add items to B. The received symbols would not have them removed, so
    /// this panics after the first symbol.
    pub fn insert(&mut self, items: &[[u8; L]]) {
        assert!(self.cells.is_empty(), "insertion after a coded symbol");
        let r = self.riblt;
        for (x, a) in items.iter().zip(r.addends(items)) {
            self.local_set
                .push(*x, r.group.neg_addend(&a), r.mapping(x));
        }
    }

    /// Receive A's next coded symbol: subtract B, apply the recovered items,
    /// and queue the result if it is decodable.
    pub fn receive(&mut self, mut c: Cell<G, L>) {
        let r = self.riblt;
        let g = &r.group;
        self.local_set.apply(g, &mut c);
        self.remote.apply(g, &mut c);
        self.local.apply(g, &mut c);
        let i = self.cells.len();
        match c.count {
            0 if g.is_identity(&c.sum) => self.decodable.push((i, None)),
            1 | -1 => {
                if let Some(a) = pure(r, &c) {
                    self.decodable.push((i, Some(a)));
                }
            }
            _ => {}
        }
        self.cells.push(c);
    }

    /// Apply the recovered item `x` with signed addend `a` and count change
    /// `sign` to the received symbols it maps to, queueing those that become
    /// pure, and return its mapping at the next symbol to be received.
    ///
    /// A symbol queued with count +1 or -1 holds one item, so a later change
    /// can only remove that item and leave count 0; queueing only counts of
    /// +1 and -1 therefore queues no symbol twice, as upstream argues.
    fn recover(&mut self, x: &[u8; L], a: &G::Addend, sign: i64) -> Mapping<P> {
        let r = self.riblt;
        let mut m = r.mapping(x);
        while m.index() < self.cells.len() as u64 {
            let i = m.index() as usize;
            let c = &mut self.cells[i];
            update(&r.group, c, x, a, sign);
            if c.count.abs() == 1
                && let Some(b) = pure(r, c)
            {
                self.decodable.push((i, Some(b)));
            }
            m.advance();
        }
        m
    }

    /// Decode the queued symbols, and those that recovered items make pure.
    pub fn try_decode(&mut self) {
        let mut k = 0;
        while k < self.decodable.len() {
            let (i, a) = self.decodable[k];
            k += 1;
            let c = self.cells[i];
            match (c.count, a) {
                // emptied by an item recovered from another symbol
                (0, _) => {}
                (1 | -1, Some(a)) => {
                    // An item of A \ B leaves the symbols; one of B \ A,
                    // already subtracted with B, returns to them.
                    let g = &self.riblt.group;
                    let signed = if c.count > 0 { g.neg_addend(&a) } else { a };
                    let m = self.recover(&c.key, &signed, -c.count);
                    let w = if c.count > 0 {
                        &mut self.remote
                    } else {
                        &mut self.local
                    };
                    w.push(c.key, signed, m);
                }
                // Only a checksum collision admits an impure symbol whose
                // count later leaves -1..=1.
                _ => panic!("decodable symbol {i} has count {}", c.count),
            }
            self.decoded += 1;
        }
        self.decodable.clear();
    }

    /// Whether every symbol received so far has been decoded.
    pub fn decoded(&self) -> bool {
        self.decoded == self.cells.len()
    }

    /// The number of coded symbols received.
    pub fn received(&self) -> usize {
        self.cells.len()
    }

    /// The items recovered so far: A \ B, and B \ A.
    pub fn difference(&self) -> Difference<L> {
        let keys = |w: &Window<G, P, L>| w.items.iter().map(|(x, _)| *x).collect();
        (keys(&self.remote), keys(&self.local))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curve::binary127;
    use crate::field::gf2_127;
    use crate::group::Group;
    use crate::hash::Salted;
    use crate::riblt::{ChaCha8, Mcg64};
    use proptest::prelude::*;

    /// XOR of 64-bit salted digests: a checksum without group structure.
    #[derive(Clone, Copy, Debug)]
    struct Xor;

    impl Group for Xor {
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

    impl Accumulate for Xor {
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

    impl Negate for Xor {
        fn neg(&self, a: &u64) -> u64 {
            *a
        }
        fn neg_addend(&self, a: &u64) -> u64 {
            *a
        }
    }

    impl HashToCurve for Xor {
        fn hash(&self, h: &Salted, msg: &[u8]) -> u64 {
            u64::from_le_bytes(h.digest(msg, 0)[..8].try_into().unwrap())
        }
    }

    fn bin() -> binary127::Curve {
        binary127::Curve::new(gf2_127::from_u128(
            0x1234_5678_9abc_def0_1234_5678_9abc_def1,
        ))
    }

    /// Go's test symbols: a little-endian counter in 64 bytes.
    fn items(from: u64, n: u64) -> Vec<[u8; 64]> {
        (from..from + n)
            .map(|i| {
                let mut x = [0; 64];
                x[..8].copy_from_slice(&i.to_le_bytes());
                x
            })
            .collect()
    }

    /// Stream A = common + only_a to a decoder holding B = common + only_b
    /// until it decodes; return the sorted difference and the symbols used.
    fn reconcile<G: HashToCurve + Negate, P: Prng>(
        r: &Riblt<G, P>,
        common: u64,
        only_a: u64,
        only_b: u64,
    ) -> (Difference<64>, usize) {
        let (c, a, b) = (
            items(0, common),
            items(common, only_a),
            items(common + only_a, only_b),
        );
        let mut enc = Encoder::new(r);
        let mut dec = Decoder::new(r);
        enc.insert(&[a, c.clone()].concat());
        dec.insert(&[b, c].concat());
        loop {
            dec.receive(enc.produce());
            dec.try_decode();
            if dec.decoded() {
                break;
            }
            assert!(dec.received() < 1 << 20, "no decode");
        }
        let (mut a, mut b) = dec.difference();
        a.sort();
        b.sort();
        ((a, b), dec.received())
    }

    fn check<G: HashToCurve + Negate, P: Prng>(
        r: &Riblt<G, P>,
        common: u64,
        only_a: u64,
        only_b: u64,
    ) -> usize {
        let (got, n) = reconcile(r, common, only_a, only_b);
        let (mut a, mut b) = (items(common, only_a), items(common + only_a, only_b));
        a.sort();
        b.sort();
        let want = (a, b);
        assert_eq!(got, want, "common {common}, a {only_a}, b {only_b}");
        n
    }

    /// Go's BenchmarkEncodeAndDecode sets: d/2 items on each side, d common.
    fn decodes_go_differences<G: HashToCurve + Negate>(g: G, ds: &[u64]) {
        let r = Riblt::new(g, &[9; 32]);
        let r64 = Riblt::<G, Mcg64>::with_prng(g, &[9; 32]);
        let rc = Riblt::<G, ChaCha8>::with_prng(g, &[9; 32]);
        for &d in ds {
            check(&r, d, d / 2, d - d / 2);
            check(&r64, d, d / 2, d - d / 2);
            check(&rc, d, d / 2, d - d / 2);
        }
    }

    #[test]
    fn decodes_exactly_the_difference() {
        decodes_go_differences(Xor, &[1, 2, 3, 10, 20, 40, 100, 1000]);
        decodes_go_differences(bin(), &[1, 2, 10, 40, 100]);
    }

    #[test]
    fn decodes_one_sided_differences() {
        let r = Riblt::new(bin(), &[4; 32]);
        for (a, b) in [(0, 1), (1, 0), (0, 30), (30, 0)] {
            check(&r, 20, a, b);
        }
    }

    /// Equal sets, empty or not, decode from the first symbol: its count is
    /// 0 and its checksum the identity.
    #[test]
    fn equal_sets_decode_after_one_symbol() {
        let r = Riblt::new(bin(), &[4; 32]);
        assert_eq!(check(&r, 0, 0, 0), 1);
        assert_eq!(check(&r, 50, 0, 0), 1);
    }

    /// The overhead approaches upstream's 1.35 symbols per difference.
    #[test]
    fn symbols_per_difference_near_upstream() {
        let r = Riblt::new(Xor, &[2; 32]);
        let d = 10_000;
        let n = check(&r, d, d / 2, d / 2);
        let ratio = n as f64 / d as f64;
        assert!((1.25..1.5).contains(&ratio), "{ratio}");
    }

    /// Coded symbol i is cell i of the finite simulator for any m > i.
    #[test]
    fn coded_symbols_are_the_finite_cells() {
        let r = Riblt::new(bin(), &[6; 32]);
        let g = r.group;
        let xs = items(0, 200);
        let mut enc = Encoder::new(&r);
        enc.insert(&xs);
        let symbols: Vec<_> = (0..300).map(|_| enc.produce()).collect();
        let mut cells = r.cells(300);
        r.encode(&mut cells, &xs, 1);
        for (s, c) in symbols.iter().zip(&cells) {
            assert_eq!((s.key, s.count), (c.key, c.count));
            assert!(g.is_identity(&g.add_affine(&s.sum, &g.neg(&g.to_affine(&c.sum)))));
        }
    }

    #[test]
    #[should_panic(expected = "insertion after a coded symbol")]
    fn encoder_rejects_insertion_after_a_symbol() {
        let r = Riblt::new(Xor, &[0; 32]);
        let mut enc = Encoder::<_, _, 64>::new(&r);
        enc.insert(&items(0, 3));
        enc.produce();
        enc.insert(&items(3, 1));
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(16))]
        #[test]
        fn decodes_random_differences(
            salt in any::<[u8; 32]>(),
            common in 0u64..60,
            a in 0u64..40,
            b in 0u64..40,
        ) {
            check(&Riblt::new(Xor, &salt), common, a, b);
            check(&Riblt::new(bin(), &salt), common, a, b);
        }
    }
}
