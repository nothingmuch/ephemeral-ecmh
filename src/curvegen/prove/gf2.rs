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
use crate::curvegen::select109;

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

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

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

    proptest! {}
}
