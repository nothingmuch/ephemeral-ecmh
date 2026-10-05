//! Maps from digest halves to curve points, and the drivers that apply
//! them. A map belongs to a curve model over a field, not to a group: it
//! needs only the model's equation in field elements, and every group
//! built on that model uses it. Which subgroup or quotient the result is
//! read in, and with what distribution, is the group's hash
//! (`group::HashToCurve`).
//!
//! `Lift` separates candidate construction from point recovery. Try-and-increment
//! examines the 128-bit halves of `Salted::digest(msg, ctr)` in the order
//! (0, low), (0, high), (1, low), (1, high), and so on, returning the first
//! successfully lifted point. `Map` maps a single digest half without retries.
//!
//! Maps by curve model, over the fields this crate instantiates:
//!
//! | model                     | fields                          | T&I | Elligator 2 | SSWU | Pornin |
//! |---------------------------|---------------------------------|-----|-------------|------|--------|
//! | Edwards, a = 1            | F_p (127)                       | yes | yes         | (1)  | (3)    |
//! | short Weierstrass, a = -3 | F_p (127)                       | yes | (2)         | yes  | (3)    |
//! | binary, a = 1             | GF(2^127)                       | yes | (4)         | (4)  | yes    |
//!
//! 1. SSWU applies to y^2 = x^3 + Ax + B with AB != 0. These curves have a
//!    point of order 2, hence a Montgomery model, to which Elligator 2
//!    applies directly; RFC 9380 (§6.7.1) uses it for such curves. SSWU
//!    would need a rational map to a Weierstrass model on top of the same
//!    exponentiation.
//! 2. Elligator 2 needs a point of order 2. The Weierstrass families have
//!    prime order.
//! 3. Pornin's map is a characteristic-2 construction.
//! 4. RFC 9380's Elligator 2 and SSWU are odd-characteristic maps. The
//!    characteristic-2 counterpart of their Shallue–van de Woestijne
//!    ancestor is binary Elligator squared (Aranha et al., SAC 2014);
//!    Pornin's map is the deterministic map evaluated here.
//!
//! Elligator 2 is `Elligator2` over the `Montgomery` trait; SSWU is
//! `curve::weier::Sswu`; Pornin's map is `curve::binary`'s `Map`, with
//! per-modulus constants (`binary::Pornin`). A single map is not uniform;
//! the group hashes use try-and-increment.

use crate::field::batch::{Invert, invert};
use crate::hash::{Salted, halves};

mod elligator2;

pub use elligator2::{Elligator2, Montgomery};

/// Candidate construction and partial lifting from 128-bit digest halves to `P`.
pub trait Lift<P> {
    type Candidate: Copy;
    /// Construct a candidate from `c`, or return `None` to skip this digest half.
    fn candidate(&self, c: u128) -> Option<Self::Candidate>;
    /// Recover a point from an admitted candidate, or return `None` if lifting fails.
    fn lift(&self, c: Self::Candidate) -> Option<P>;
}

/// A lift that exposes one denominator per candidate for shared inversion.
pub trait BatchLift<P>: Lift<P> {
    type F: Invert;
    /// Return the nonzero denominator for a candidate admitted by `candidate`.
    fn denominator(&self, c: Self::Candidate) -> Self::F;
    /// Return the same result as `lift(c)`, given `inv = 1 / denominator(c)`.
    fn lift_inv(&self, c: Self::Candidate, inv: Self::F) -> Option<P>;
}

/// A deterministic, total map from 128-bit digest halves to points of `P`.
/// Each implementation defines its field reduction and any sign-bit convention.
/// Use of an RFC 9380 map-to-curve operation does not by itself implement an
/// RFC 9380 hash-to-curve suite or imply a uniform distribution on the group.
pub trait Map<P> {
    fn map(&self, c: u128) -> P;
}

/// Apply `map` to the low 128-bit half of `h.digest(msg, 0)`.
pub fn map1<P, M: Map<P>>(m: &M, h: &Salted, msg: &[u8]) -> P {
    m.map(halves(&h.digest(msg, 0))[0])
}

/// Return the first point accepted by the candidate and lift operations in
/// digest order. The supplied `Salted` value determines the hash namespace.
pub fn try_and_increment<P, L: Lift<P>>(l: &L, h: &Salted, msg: &[u8]) -> P {
    for ctr in 0.. {
        for c in halves(&h.digest(msg, ctr)) {
            if let Some(p) = l.candidate(c).and_then(|c| l.lift(c)) {
                return p;
            }
        }
    }
    unreachable!()
}

/// Apply `try_and_increment` independently to each message, preserving input order.
/// Each round shares one inversion across the current pending candidates.
pub fn try_and_increment_batch<P: Clone, L: BatchLift<P>>(
    l: &L,
    h: &Salted,
    msgs: &[&[u8]],
) -> Vec<P> {
    /// The current digest half for one message: half `k` at counter `ctr`.
    struct Try {
        item: usize,
        ctr: u32,
        halves: [u128; 2],
        k: usize,
    }
    impl Try {
        /// Advance to the first remaining digest half admitted by `candidate`.
        fn seek<P, L: Lift<P>>(&mut self, l: &L, h: &Salted, msg: &[u8]) -> L::Candidate {
            loop {
                if self.k == 2 {
                    self.ctr += 1;
                    self.halves = halves(&h.digest(msg, self.ctr));
                    self.k = 0;
                }
                if let Some(c) = l.candidate(self.halves[self.k]) {
                    return c;
                }
                self.k += 1;
            }
        }
    }
    let mut out = vec![None; msgs.len()];
    let mut pending: Vec<(Try, L::Candidate)> = (0..msgs.len())
        .map(|item| {
            let mut t = Try {
                item,
                ctr: 0,
                halves: halves(&h.digest(msgs[item], 0)),
                k: 0,
            };
            let c = t.seek(l, h, msgs[item]);
            (t, c)
        })
        .collect();
    while !pending.is_empty() {
        let mut inv: Vec<L::F> = pending.iter().map(|&(_, c)| l.denominator(c)).collect();
        invert(&mut inv);
        let mut next = Vec::with_capacity(pending.len());
        for ((mut t, c), inv) in pending.into_iter().zip(inv) {
            if let Some(p) = l.lift_inv(c, inv) {
                out[t.item] = Some(p);
                continue;
            }
            t.k += 1;
            let c = t.seek(l, h, msgs[t.item]);
            next.push((t, c));
        }
        pending = next;
    }
    out.into_iter().map(Option::unwrap).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::Binary;
    use crate::field::gf2_127::{Gf, eq};
    use proptest::prelude::*;

    /// Reject digest halves whose low two bits are zero, then require bit 8
    /// for lifting. This exercises separate candidate and lift rejection paths.
    struct Mock;

    fn passes(c: u128) -> bool {
        c & 3 != 0
    }

    impl Lift<u128> for Mock {
        type Candidate = u128;
        fn candidate(&self, c: u128) -> Option<u128> {
            passes(c).then_some(c)
        }
        fn lift(&self, c: u128) -> Option<u128> {
            assert!(passes(c), "a skipped half reached lift");
            (c >> 8 & 1 == 1).then_some(c)
        }
    }

    impl BatchLift<u128> for Mock {
        type F = Gf;
        fn denominator(&self, c: u128) -> Gf {
            assert!(passes(c), "a skipped half reached denominator");
            Gf::new(c >> 64 | 1)
        }
        fn lift_inv(&self, c: u128, inv: Gf) -> Option<u128> {
            assert!(eq(inv * self.denominator(c), Gf::ONE));
            self.lift(c)
        }
    }

    fn reference(h: &Salted, msg: &[u8]) -> u128 {
        (0..)
            .flat_map(|ctr| halves(&h.digest(msg, ctr)))
            .find(|&c| passes(c) && c >> 8 & 1 == 1)
            .unwrap()
    }

    proptest! {
        #[test]
        fn drivers_match_the_reference(
            salt in any::<[u8; 32]>(),
            // Short messages over a small alphabet, so batches repeat some.
            msgs in prop::collection::vec(prop::collection::vec(0u8..4, 0..3), 0..40),
        ) {
            let h = Salted::new(b"test", &salt);
            let refs: Vec<&[u8]> = msgs.iter().map(|m| m.as_slice()).collect();
            let want: Vec<u128> = refs.iter().map(|m| reference(&h, m)).collect();
            let single: Vec<u128> = refs.iter().map(|m| try_and_increment(&Mock, &h, m)).collect();
            prop_assert_eq!(&single, &want);
            prop_assert_eq!(try_and_increment_batch(&Mock, &h, &refs), want);
        }
    }

    #[test]
    fn empty_batch() {
        let h = Salted::new(b"test", &[0; 32]);
        assert!(try_and_increment_batch(&Mock, &h, &[]).is_empty());
    }
}
