//! Salted candidate streams: 32 bytes per message and counter.
//!
//! `Salted::new` is tagged, salted SHA-256 with a per-salt midstate:
//! `SHA256(SHA256(tag) || SHA256(tag) || salt || 0^32 || msg || ctr_le32)`.
//! The first 128 bytes are exactly two compression blocks, so each salt
//! costs two compressions once, and a message of up to 51 bytes costs one.
//!
//! `Salted::projected` divides that work by lifetime. An item's ID, its
//! unsalted SHA-256 (`id`), outlives every salt; a salt then needs only a
//! keyed universal hash of the 32-byte ID. The ID's 128-bit halves m1, m2
//! are field elements, and each output half is U_k(m1, m2) = m1 k^2 + m2 k
//! under its own key k, as Poly1305 evaluates a two-block message but with
//! neither clamping nor a one-time pad. For two IDs fixed before the salt,
//! U_k(a) - U_k(b) is a nonzero polynomial of degree at most 2 in k, so they
//! collide for at most 2 keys: probability at most 2^-126. After the salt, colliding IDs
//! cost a generic search, as they would for any 128-bit hash, provided
//! the IDs are hash outputs: an adversary choosing m1, m2 directly solves
//! for collisions under a known key.
//!
//! The field is a matter of cost. Items reach the curve independently and
//! are added only as points, so linearity over a curve's own field yields
//! no relation among points short of a symmetry of the map, which needs
//! IDs agreeing up to sign in about 254 bits (docs/problem.md, Adversary).
//! GF(2^127) is cheapest with carry-less multiplication, F_(2^127 - 1)
//! without; F_(2^130 - 5) is Poly1305's. Each half is close to
//! uniform; the two of a counter are not independent: given one, the
//! other takes about 2^-2 of its values in F_(2^130 - 5), and in the
//! 127-bit fields each carries one of the ID's two leftover top bits,
//! which makes the pair a bijection of the ID (in F_(2^127 - 1), up to
//! the all-ones pattern, which reads as 0).
//! Later counters are fresh keys over the same ID, so their halves are
//! linear in counter 0's: a try-and-increment that reaches them draws
//! candidates correlated with the ones it rejected.
//!
//! The keys come from the salted SHA-256 under the tag, two per counter,
//! nonzero and distinct.

use core::ops::{Add, Mul};

use sha2::{Digest, Sha256};

use crate::field::fp127::Fp;
use crate::field::gf2_127::{self, Gf, MASK127};

/// The field a projection evaluates in. F_(2^130 - 5) and F_(2^127 - 1)
/// multiply with integer multiplication on every target, GF(2^127) by
/// carry-less multiplication where the build has it.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Field {
    Fp130,
    Fp127,
    Gf2_127,
}

#[derive(Clone)]
pub struct Salted(Mode);

#[derive(Clone)]
enum Mode {
    Sha256(Sha256),
    Fp130(Box<Projection<Fe130>>),
    Fp127(Box<Projection<Fp>>),
    Gf2_127(Box<Projection<Gf>>),
}

/// Counters whose keys are kept; try-and-increment rarely passes them.
const KEPT: usize = 4;

trait Element: Copy + Add<Output = Self> + Mul<Output = Self> {
    /// The element a 128-bit half reads as, and the bits it leaves over.
    fn split(v: u128) -> (Self, u128);
    /// The low 128 bits of the canonical representative.
    fn bits(self) -> u128;
    /// a b + c d.
    fn dot(a: Self, b: Self, c: Self, d: Self) -> Self {
        a * b + c * d
    }
}

impl Element for Gf {
    fn split(v: u128) -> (Self, u128) {
        (gf2_127::from_u128(v), v & !MASK127)
    }
    fn bits(self) -> u128 {
        gf2_127::to_u128(self)
    }
}

impl Element for Fp {
    fn split(v: u128) -> (Self, u128) {
        (Fp::new(v & MASK127), v & !MASK127)
    }
    fn bits(self) -> u128 {
        self.value()
    }
}

const M44: u64 = (1 << 44) - 1;
const M42: u64 = (1 << 42) - 1;

/// F_(2^130 - 5) in poly1305-donna-64's radix: limbs of 44, 44 and 42
/// bits, partly reduced, so a product's terms stay below 2^92. 2^132 is
/// 20 modulo p, which folds the high products.
#[derive(Clone, Copy, Debug)]
struct Fe130([u64; 3]);

impl Fe130 {
    fn from_u128(v: u128) -> Self {
        Self([v as u64 & M44, (v >> 44) as u64 & M44, (v >> 88) as u64])
    }

    fn carry(d: [u128; 3]) -> Self {
        let mut h = [0u64; 3];
        let d1 = d[1] + (d[0] >> 44);
        let d2 = d[2] + (d1 >> 44);
        h[0] = d[0] as u64 & M44;
        h[1] = d1 as u64 & M44;
        h[2] = d2 as u64 & M42;
        let h0 = h[0] + (d2 >> 42) as u64 * 5;
        h[0] = h0 & M44;
        h[1] += h0 >> 44;
        Self(h)
    }

    /// The canonical representative, in [0, p).
    fn canonical(self) -> [u64; 3] {
        let Self([mut h0, mut h1, mut h2]) = Self::carry(self.0.map(u128::from));
        h1 += h0 >> 44;
        h0 &= M44;
        h2 += h1 >> 44;
        h1 &= M44;
        h0 += (h2 >> 42) * 5;
        h2 &= M42;
        h1 += h0 >> 44;
        h0 &= M44;
        // h < 2^130 now; h - p = h + 5 - 2^130 when that does not borrow
        let g0 = h0 + 5;
        let g1 = h1 + (g0 >> 44);
        let g2 = h2 + (g1 >> 44);
        if g2 >> 42 != 0 {
            [g0 & M44, g1 & M44, g2 & M42]
        } else {
            [h0, h1 & M44, h2 + (h1 >> 44)]
        }
    }
}

impl Add for Fe130 {
    type Output = Self;
    fn add(self, r: Self) -> Self {
        Self([0, 1, 2].map(|i| self.0[i] + r.0[i]))
    }
}

impl Fe130 {
    /// The unreduced product terms of h r.
    fn terms(h: Self, r: Self) -> [u128; 3] {
        let [h0, h1, h2] = h.0.map(u128::from);
        let [r0, r1, r2] = r.0.map(u128::from);
        let (s1, s2) = (r1 * 20, r2 * 20);
        [
            h0 * r0 + h1 * s2 + h2 * s1,
            h0 * r1 + h1 * r0 + h2 * s2,
            h0 * r2 + h1 * r1 + h2 * r0,
        ]
    }
}

impl Mul for Fe130 {
    type Output = Self;
    fn mul(self, r: Self) -> Self {
        Self::carry(Self::terms(self, r))
    }
}

impl Element for Fe130 {
    fn split(v: u128) -> (Self, u128) {
        (Self::from_u128(v), 0)
    }
    fn bits(self) -> u128 {
        let [h0, h1, h2] = self.canonical().map(u128::from);
        h0 | h1 << 44 | h2 << 88
    }
    /// One carry for both products: their terms sum below 2^94.
    fn dot(a: Self, b: Self, c: Self, d: Self) -> Self {
        let (x, y) = (Self::terms(a, b), Self::terms(c, d));
        Self::carry([x[0] + y[0], x[1] + y[1], x[2] + y[2]])
    }
}

/// Each counter's two keys, as (k^2, k).
#[derive(Clone)]
struct Projection<F> {
    kept: [[(F, F); 2]; KEPT],
    keys: Sha256,
}

impl<F: Element> Projection<F> {
    fn new(keys: Sha256) -> Self {
        let kept = core::array::from_fn(|ctr| Self::derive(&keys, ctr as u32));
        Self { kept, keys }
    }

    /// A counter's keys: the halves of the salted SHA-256 of its bytes,
    /// rederived with a retry counter appended while a key is 0 or the
    /// two are equal, either of which would collapse the projection.
    fn derive(keys: &Sha256, ctr: u32) -> [(F, F); 2] {
        for retry in 0u32.. {
            let mut h = keys.clone();
            h.update(ctr.to_le_bytes());
            if retry > 0 {
                h.update(retry.to_le_bytes());
            }
            let d: [u8; 32] = h.finalize().into();
            let [k0, k1] = halves(&d).map(|k| F::split(k).0);
            let [b0, b1] = [k0, k1].map(F::bits);
            if b0 != 0 && b1 != 0 && b0 != b1 {
                return [(k0 * k0, k0), (k1 * k1, k1)];
            }
        }
        unreachable!()
    }

    fn keys(&self, ctr: u32) -> [(F, F); 2] {
        match self.kept.get(ctr as usize) {
            Some(k) => *k,
            None => Self::derive(&self.keys, ctr),
        }
    }

    fn digest(&self, id: &[u8], ctr: u32) -> [u8; 32] {
        let keys = self.keys(ctr);
        let [(m1, t1), (m2, t2)] = split(id);
        let mut out = [0; 32];
        for ((o, &(kk, k)), top) in out
            .as_chunks_mut::<16>()
            .0
            .iter_mut()
            .zip(&keys)
            .zip([t1, t2])
        {
            o.copy_from_slice(&(F::dot(m1, kk, m2, k).bits() | top).to_le_bytes());
        }
        out
    }

    fn half(&self, id: &[u8], ctr: u32, i: usize) -> u128 {
        let (kk, k) = self.keys(ctr)[i];
        let [(m1, t1), (m2, t2)] = split(id);
        F::dot(m1, kk, m2, k).bits() | [t1, t2][i]
    }
}

fn split<F: Element>(id: &[u8]) -> [(F, u128); 2] {
    let id: &[u8; 32] = id.try_into().expect("a projection hashes 32-byte IDs");
    halves(id).map(F::split)
}

fn midstate(tag: &[u8], salt: &[u8; 32]) -> Sha256 {
    let t = Sha256::digest(tag);
    let mut h = Sha256::new();
    h.update(t);
    h.update(t);
    h.update(salt);
    h.update([0u8; 32]);
    h
}

impl Salted {
    pub fn new(tag: &[u8], salt: &[u8; 32]) -> Self {
        Self(Mode::Sha256(midstate(tag, salt)))
    }

    /// The projection of 32-byte IDs under `tag` and `salt`, in `field`.
    pub fn projected(tag: &[u8], salt: &[u8; 32], field: Field) -> Self {
        let keys = midstate(tag, salt);
        Self(match field {
            Field::Fp130 => Mode::Fp130(Box::new(Projection::new(keys))),
            Field::Fp127 => Mode::Fp127(Box::new(Projection::new(keys))),
            Field::Gf2_127 => Mode::Gf2_127(Box::new(Projection::new(keys))),
        })
    }

    /// The candidate bytes of `msg` at `ctr`; under a projection, `msg`
    /// is an ID and must be 32 bytes.
    pub fn digest(&self, msg: &[u8], ctr: u32) -> [u8; 32] {
        match &self.0 {
            Mode::Sha256(s) => {
                let mut h = s.clone();
                h.update(msg);
                h.update(ctr.to_le_bytes());
                h.finalize().into()
            }
            Mode::Fp130(p) => p.digest(msg, ctr),
            Mode::Fp127(p) => p.digest(msg, ctr),
            Mode::Gf2_127(p) => p.digest(msg, ctr),
        }
    }

    /// Half `i` of `digest(msg, ctr)`. A projection evaluates only that
    /// half's polynomial, so a map that reads one half costs one product
    /// pair; the salted SHA-256 computes the whole digest regardless.
    pub fn half(&self, msg: &[u8], ctr: u32, i: usize) -> u128 {
        match &self.0 {
            Mode::Sha256(_) => halves(&self.digest(msg, ctr))[i],
            Mode::Fp130(p) => p.half(msg, ctr, i),
            Mode::Fp127(p) => p.half(msg, ctr, i),
            Mode::Gf2_127(p) => p.half(msg, ctr, i),
        }
    }
}

/// An item's ID: its SHA-256, without tag or salt.
pub fn id(item: &[u8]) -> [u8; 32] {
    Sha256::digest(item).into()
}

/// The two 128-bit little-endian halves of a digest.
pub fn halves(d: &[u8; 32]) -> [u128; 2] {
    let (lo, hi) = d.split_at(16);
    [
        u128::from_le_bytes(lo.try_into().unwrap()),
        u128::from_le_bytes(hi.try_into().unwrap()),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use poly1305::Poly1305;
    use poly1305::universal_hash::{KeyInit, UniversalHash};
    use proptest::prelude::*;

    const ALL: [Field; 3] = [Field::Fp130, Field::Fp127, Field::Gf2_127];

    fn keys(salt: &[u8; 32], ctr: u32) -> [u128; 2] {
        halves(&Salted::new(b"t", salt).digest(&[], ctr))
    }

    /// Poly1305 under key (r, s = 0) of the two blocks m1, m2 is
    /// (m1 + 2^128) r^2 + (m2 + 2^128) r mod p, read modulo 2^128.
    fn poly1305(r: u128, id: &[u8; 32]) -> u128 {
        let mut key = [0u8; 32];
        key[..16].copy_from_slice(&r.to_le_bytes());
        let mut p = Poly1305::new(&key.into());
        p.update_padded(id);
        u128::from_le_bytes(p.finalize().into())
    }

    /// Poly1305's clamping: the keys it evaluates under.
    fn clamp(r: u128) -> u128 {
        r & 0x0ffffffc_0ffffffc_0ffffffc_0fffffff
    }

    proptest! {
        #[test]
        fn midstate_matches_one_shot(
            tag in prop::collection::vec(any::<u8>(), 0..80),
            salt in any::<[u8; 32]>(),
            msg in prop::collection::vec(any::<u8>(), 0..200),
            ctr in any::<u32>(),
        ) {
            let t = Sha256::digest(&tag);
            let one_shot: [u8; 32] = Sha256::new()
                .chain_update(t)
                .chain_update(t)
                .chain_update(salt)
                .chain_update([0u8; 32])
                .chain_update(&msg)
                .chain_update(ctr.to_le_bytes())
                .finalize()
                .into();
            prop_assert_eq!(Salted::new(&tag, &salt).digest(&msg, ctr), one_shot);
        }

        /// The F_(2^130 - 5) arithmetic against Poly1305's, under clamped
        /// keys: U_k(m1, m2) + 2^128 (k^2 + k) is Poly1305's tag.
        #[test]
        fn fp130_is_poly1305s_field(r in any::<u128>(), id in any::<[u8; 32]>()) {
            let k = Fe130::from_u128(clamp(r));
            let [m1, m2] = halves(&id).map(Fe130::from_u128);
            let pad = Fe130([0, 0, 1 << 40]);
            let tag = m1 * (k * k) + m2 * k + pad * (k * k + k);
            prop_assert_eq!(tag.bits(), poly1305(clamp(r), &id));
        }

        /// Products of unreduced sums stay in range: (a + b) c canonicalizes
        /// to ac + bc for elements as `digest` forms them.
        #[test]
        fn fp130_distributes(a in any::<u128>(), b in any::<u128>(), c in any::<u128>()) {
            let [a, b, c] = [a, b, c].map(Fe130::from_u128);
            let c = c * c;
            prop_assert_eq!(((a + b) * c).bits(), (a * c + b * c).bits());
            prop_assert_eq!(((a + b) * c).canonical(), (a * c + b * c).canonical());
        }

        /// Each half is m1 k^2 + m2 k under the counter's keys, the halves
        /// of the salted SHA-256 of the empty message, with the 127-bit
        /// fields' leftover top bits carried through.
        #[test]
        fn projection_matches_its_definition(
            salt in any::<[u8; 32]>(),
            id in any::<[u8; 32]>(),
            ctr in 0u32..2 * KEPT as u32,
        ) {
            let ks = keys(&salt, ctr);
            let [h1, h2] = halves(&id);
            let gf = gf2_127::from_u128;
            let want: Vec<u128> = ks.iter().zip([h1, h2]).map(|(&k, top)| {
                let k = gf(k);
                gf2_127::to_u128(gf(h1) * k * k + gf(h2) * k) | (top & !MASK127)
            }).collect();
            let got = halves(&Salted::projected(b"t", &salt, Field::Gf2_127).digest(&id, ctr));
            prop_assert_eq!(got.to_vec(), want);
            let want: Vec<u128> = ks.iter().zip([h1, h2]).map(|(&k, top)| {
                let k = Fp::new(k & MASK127);
                let (m1, m2) = (Fp::new(h1 & MASK127), Fp::new(h2 & MASK127));
                (m1 * k * k + m2 * k).value() | (top & !MASK127)
            }).collect();
            let got = halves(&Salted::projected(b"t", &salt, Field::Fp127).digest(&id, ctr));
            prop_assert_eq!(got.to_vec(), want);
            let want: Vec<u128> = ks.iter().map(|&k| {
                let k = Fe130::from_u128(k);
                let (m1, m2) = (Fe130::from_u128(h1), Fe130::from_u128(h2));
                (m1 * (k * k) + m2 * k).bits()
            }).collect();
            let got = halves(&Salted::projected(b"t", &salt, Field::Fp130).digest(&id, ctr));
            prop_assert_eq!(got.to_vec(), want);
        }

        #[test]
        fn half_is_the_digests_half(
            salt in any::<[u8; 32]>(),
            id in any::<[u8; 32]>(),
            ctr in 0u32..2 * KEPT as u32,
        ) {
            let modes = ALL.map(|f| Salted::projected(b"t", &salt, f));
            for s in modes.iter().chain([&Salted::new(b"t", &salt)]) {
                let d = halves(&s.digest(&id, ctr));
                prop_assert_eq!([0, 1].map(|i| s.half(&id, ctr, i)), d);
            }
        }

        #[test]
        fn distinct_ids_project_apart(
            salt in any::<[u8; 32]>(),
            a in any::<[u8; 32]>(),
            b in any::<[u8; 32]>(),
        ) {
            prop_assume!(a != b);
            for field in ALL {
                let s = Salted::projected(b"t", &salt, field);
                prop_assert_ne!(s.digest(&a, 0), s.digest(&b, 0));
            }
        }
    }

    /// p - 1 and the values just past p, where the final subtraction acts.
    #[test]
    fn fp130_canonical_at_the_modulus() {
        let p = [M44 - 4, M44, M42];
        assert_eq!(Fe130(p).canonical(), [0, 0, 0]);
        assert_eq!(Fe130([M44 - 5, M44, M42]).canonical(), [M44 - 5, M44, M42]);
        assert_eq!(Fe130([M44, M44, M42]).canonical(), [4, 0, 0]);
    }

    #[test]
    fn projection_rejects_other_lengths() {
        for field in ALL {
            let s = Salted::projected(b"t", &[0; 32], field);
            assert!(std::panic::catch_unwind(|| s.digest(&[0; 36], 0)).is_err());
        }
    }
}
