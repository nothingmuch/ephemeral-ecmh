//! Prover adapters for the 109- and 122-bit binary families (`select109`,
//! `select122`), and their AGM counts.
//!
//! In all three, Tr(a2) = 1 makes #E = 2 mod 4, so 2E is the subgroup of
//! odd order #E/2, and an accepted curve has #E = 2r. Witnesses are
//! multiples of hashed points, which have Tr(x) = 1 and so lie in 2E: a
//! hashed point's order divides #E/2, and need not equal it. The families
//! differ in what a verifier takes for a composite r:
//!
//! - `Binary109`: the smallest prime factor of r below 2^10, else what
//!   the `Factor` supplies, else #E itself (the even marker). Sage's
//!   labels (sage/kat109.sage) are the smallest prime factor without a
//!   bound, so reproducing them takes a `Factor`.
//! - `Dense122`, `Gls122`: the smallest prime factor of r below 2^16,
//!   else #E itself (the even-order marker), as sage/kat122.sage labels
//!   them. The verdict trial-divides to 2^16, so `NoFactor` gives Sage's
//!   labels and any `Factor` a valid certificate with others.

use super::*;
use crate::curve::{binary, binary109};
use crate::curvegen::criteria::AdmissibleR;
use crate::curvegen::select::Arithmetic;
use crate::curvegen::{select109, select122};
use select122::{Dense122, Gls122};
use std::sync::OnceLock;

/// A point of order l on a binary curve of order n = 2 * odd, from hashed
/// points of the odd part.
fn witness_point<M: binary::Model>(
    c: &binary::Curve<M>,
    h: &Salted,
    j: u32,
    n: u128,
    l: u128,
) -> binary::Affine<M> {
    torsion_point(
        n / 2,
        l,
        |i| c.hash_to_curve(h, &msg(j, i)),
        |p, k| c.to_affine(&c.mul(&c.from_affine(p), k)),
        binary::Affine::is_identity,
    )
}

/// The 14-byte encoding at the prover's 16: zero-padded, as `Rejection`
/// holds it (`certificate109` narrows it back).
impl Arithmetic<16> for binary109::Curve {
    type Affine = binary109::Affine;
    type Point = binary109::Point;

    fn decode(&self, enc: [u8; 16]) -> Option<Self::Affine> {
        let (enc, pad) = enc.split_at(binary109::ENCODED_LEN);
        (pad == [0; 2]).then(|| self.decode(enc.try_into().unwrap()))?
    }
    fn hash(&self, salt: &Salted, msg: &[u8]) -> Self::Affine {
        self.hash_to_curve(salt, msg)
    }
    fn lift(&self, p: &Self::Affine) -> Option<Self::Point> {
        (!p.is_identity()).then(|| self.from_affine(p))
    }
    fn mul(&self, p: &Self::Point, k: u128) -> Self::Point {
        binary::Curve::mul(self, p, k)
    }
    fn is_identity(&self, p: &Self::Point) -> bool {
        binary::Point::is_identity(p)
    }
}

/// `select109::Binary109`, whose certificate holds 14-byte witnesses, at
/// the prover's 16-byte width. #E = 2r.
pub struct Binary109;

type Narrow = select109::Binary109;
const NARROW: usize = binary109::ENCODED_LEN;

impl Criteria for Binary109 {
    type Curve = binary109::Curve;
    const TAG: &'static [u8] = <Narrow as Criteria>::TAG;
    const R: AdmissibleR = <Narrow as Criteria>::R;
    fn candidate(seed: &[u8; 32], j: u32) -> Option<Self::Curve> {
        Narrow::candidate(seed, j)
    }
}

impl Policy<16> for Binary109 {
    const CLEAR: u128 = <Narrow as Policy<NARROW>>::CLEAR;
    const ENTRY_PER_INDEX: bool = <Narrow as Policy<NARROW>>::ENTRY_PER_INDEX;

    fn marker(c: &Self::Curve, p: &binary109::Affine, l: u128) -> Option<bool> {
        Narrow::marker(c, p, l)
    }
}

impl Family for Binary109 {
    fn witness(c: &Self::Curve, h: &Salted, j: u32, n: u128, l: u128) -> [u8; 16] {
        let mut w = [0; 16];
        w[..NARROW].copy_from_slice(&witness_point(c, h, j, n, l).encode());
        w
    }

    /// Any hashed point: not O, and of odd order dividing n.
    fn order_witness(c: &Self::Curve, h: &Salted, j: u32, _: u128) -> Option<[u8; 16]> {
        let mut w = [0; 16];
        w[..NARROW].copy_from_slice(&c.hash_to_curve(h, &msg(j, 0)).encode());
        Some(w)
    }
}

/// A `prove::<Binary109>` certificate in `select109`'s 14-byte encoding.
pub fn certificate109(c: Certificate) -> select109::Certificate {
    let narrow = |(l, w): Rejection| {
        let (enc, pad) = w.split_at(NARROW);
        assert_eq!(pad, [0; 2], "not a binary109 encoding");
        (l, enc.try_into().unwrap())
    };
    select109::Certificate {
        index: c.index,
        r: c.r,
        rejections: c.rejections.into_iter().map(narrow).collect(),
    }
}

/// sage/kat122.sage's trial bound: labels below it are prime factors.
const TRIAL122: u32 = 1 << 16;

/// The odd primes below `TRIAL122`.
fn primes122() -> &'static [u32] {
    static P: OnceLock<Vec<u32>> = OnceLock::new();
    P.get_or_init(|| {
        let mut composite = vec![false; TRIAL122 as usize];
        (3..TRIAL122)
            .step_by(2)
            .filter(|&l| {
                let prime = !composite[l as usize];
                if prime {
                    (l * l..TRIAL122)
                        .step_by(2 * l as usize)
                        .for_each(|k| composite[k as usize] = true);
                }
                prime
            })
            .collect()
    })
}

/// m mod l in 64-bit arithmetic, from m's halves: l < 2^16 keeps every
/// product below 2^32.
fn rem(m: u128, l: u32) -> u64 {
    let l = u64::from(l);
    let (hi, lo) = ((m >> 64) as u64 % l, m as u64 % l);
    let wrap = (u64::MAX % l + 1) % l;
    (hi * wrap + lo) % l
}

/// r = n/2 with sage/kat122.sage's labels: its smallest prime factor below
/// 2^16, else `Composite` (for a rejection by order) or acceptance, which
/// requires the embedding degree to exceed `EMBEDDING_MIN`.
fn verdict122<F: Criteria>(n: u128) -> Verdict {
    let m = n / 2;
    if let Some(&l) = primes122().iter().find(|&&l| rem(m, l) == 0) {
        return Verdict::Reject(l.into());
    }
    if !is_prime(m) {
        return Verdict::Composite(m);
    }
    if let Err(reason) = F::check_order(n) {
        return Verdict::Inadmissible(reason);
    }
    Verdict::Accept(m)
}

macro_rules! binary122_family {
    ($name:ident) => {
        impl Family for $name {
            fn verdict(n: u128) -> Verdict {
                verdict122::<Self>(n)
            }

            fn witness(c: &Self::Curve, h: &Salted, j: u32, n: u128, l: u128) -> [u8; 16] {
                witness_point(c, h, j, n, l).encode()
            }

            /// Any hashed point: not O, and of odd order dividing n.
            fn order_witness(c: &Self::Curve, h: &Salted, j: u32, _: u128) -> Option<[u8; 16]> {
                Some(c.hash_to_curve(h, &msg(j, 0)).encode())
            }
        }
    };
}

binary122_family!(Dense122);
binary122_family!(Gls122);

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn primes_below_the_trial_bound() {
        let p = primes122();
        assert_eq!((p.len(), p[0], *p.last().unwrap()), (6541, 3, 65521));
    }

    /// The 16-byte adapter decodes a native encoding under a zero pad,
    /// nothing under a nonzero one, and `certificate109` strips the pad.
    #[test]
    fn padded_decoding_requires_a_zero_pad() {
        let seed = [7; 32];
        let c = (0..).find_map(|j| select109::candidate(&seed, j)).unwrap();
        let p = c.hash_to_curve(&Salted::new(b"test", &seed), b"pad");
        let mut w = [0; 16];
        w[..NARROW].copy_from_slice(&p.encode());
        let enc = p.encode();
        assert_eq!(
            Arithmetic::<16>::decode(&c, w).map(|a| a.encode()),
            Some(enc)
        );
        for i in [NARROW, NARROW + 1] {
            let mut bad = w;
            bad[i] = 1;
            assert!(Arithmetic::<16>::decode(&c, bad).is_none());
        }
        let cert = Certificate {
            index: 1,
            r: 0,
            rejections: vec![(3, w)],
        };
        assert_eq!(certificate109(cert).rejections, [(3, enc)]);
    }

    proptest! {
        #[test]
        fn rem_matches_u128(m in any::<u128>(), l in 1u32..TRIAL122) {
            prop_assert_eq!(u128::from(rem(m, l)), m % u128::from(l));
        }
    }
}
