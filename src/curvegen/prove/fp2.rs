//! Witnesses for the quadratic-field policies in `select_fp2`.

use super::*;
use crate::curve::{edwards61x2, twisted, weier61x2};
use crate::curvegen::select_fp2::{
    Edwards61x2, Twisted61x2, Twisted64x2, TwistedGoldilocks2, Weier61x2,
};
use crate::curvegen::{poly, sieve};
use crate::field::Packed;

/// In the quotient `E/<T>` of order n/2, label 8 is a point of order 4.
fn quotient_witness<F: Packed>(
    c: &twisted::Curve<F>,
    h: &Salted,
    j: u32,
    n: u128,
    l: u128,
) -> [u8; 16] {
    let (base, target) = if l == 8 { (2, 4) } else { (l, l) };
    prime_power_point(
        n / 2,
        base,
        target,
        |i| c.hash_to_curve(h, &msg(j, i)),
        |p, k| c.to_affine(&c.mul(&c.from_affine(p), k)),
        twisted::Affine::is_identity,
    )
    .encode()
}

macro_rules! twisted_witness {
    ($($name:ident),*) => {$(
        impl Family for $name {
            fn witness(c: &Self::Curve, h: &Salted, j: u32, n: u128, l: u128) -> [u8; 16] {
                quotient_witness(c, h, j, n, l)
            }
        }
    )*};
}

twisted_witness!(Twisted61x2, Twisted64x2, TwistedGoldilocks2);

impl Family for Edwards61x2 {
    fn witness(c: &Self::Curve, h: &Salted, j: u32, n: u128, l: u128) -> [u8; 16] {
        torsion_point(
            n,
            l,
            |i| c.hash_to_curve(h, &msg(j, i)),
            |p, k| c.to_affine(&c.mul(&c.from_affine(p), k)),
            edwards61x2::Affine::is_identity,
        )
        .encode()
    }
}

impl Family for Weier61x2 {
    /// The affine chord law is complete before the order is known:
    /// rejected candidates can have rational two-torsion.
    fn witness(c: &Self::Curve, h: &Salted, j: u32, n: u128, l: u128) -> [u8; 16] {
        torsion_point(
            n,
            l,
            |i| c.hash_to_curve(h, &msg(j, i)),
            |p, k| affine_mul(c, p, k),
            weier61x2::Affine::is_identity,
        )
        .encode()
    }
}

/// No quadratic family has a `quick_reject`, so the sieve tries 8 itself.
impl Sieve<edwards61x2::Curve> for SmallL {
    fn reject(&mut self, c: &edwards61x2::Curve) -> Option<Rejection> {
        sieve::edwards(c, self.0)
    }
}

impl Sieve<weier61x2::Curve> for SmallL {
    fn reject(&mut self, c: &weier61x2::Curve) -> Option<Rejection> {
        sieve::weier(c, self.0)
    }
}

impl<F: poly::Field + Packed> Sieve<twisted::Curve<F>> for SmallL {
    fn reject(&mut self, c: &twisted::Curve<F>) -> Option<Rejection> {
        sieve::twisted(c, self.0)
    }
}
