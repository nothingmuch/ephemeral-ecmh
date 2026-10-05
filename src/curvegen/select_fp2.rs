//! The policies over the quadratic fields: `select`'s scheme with a tag
//! and a candidate filter per family, and the bounds for q = p^2.
//!
//! The Edwards curves have order 4r and the twisted codecs represent their
//! two-torsion quotients, of order 2r; the Weierstrass comparison has
//! order r.
//!
//! Candidates whose j-invariant belongs to the base field are excluded as
//! a conservative family policy. This is not a claim that every such curve
//! has composite order, nor a concrete discrete-logarithm security proof.

use crate::curve::{
    edwards61x2, twisted, twisted_goldilocks2, twisted61x2, twisted64x2, weier61x2,
};
use crate::curvegen::criteria::{AdmissibleR, Criteria};
use crate::curvegen::select::{Certificate, Error, Policy, eight_torsion, half, two_torsion};
use crate::field::batch::Invert;
use crate::field::{Packed, fp61x2, fp64x2, goldilocks2};

pub const TAG_TWISTED61X2: &[u8] = b"ephemeral-ecmh/curve/twisted61x2";
pub const TAG_TWISTED64X2: &[u8] = b"ephemeral-ecmh/curve/twisted64x2";
pub const TAG_TWISTED_GOLDILOCKS2: &[u8] = b"ephemeral-ecmh/curve/twisted-goldilocks2";
pub const TAG_EDWARDS61X2: &[u8] = b"ephemeral-ecmh/curve/edwards61x2";
pub const TAG_WEIER61X2: &[u8] = b"ephemeral-ecmh/curve/weier61x2";

/// Over F_{p^2} the Hasse interval is [(p - 1)^2, (p + 1)^2] exactly,
/// and 2 sqrt q = 2p.
const fn admissible(p: u128, cofactor: u128) -> AdmissibleR {
    AdmissibleR::new(p * p, 2 * p, cofactor)
}

pub(super) const R61: AdmissibleR = admissible(fp61x2::P as u128, 4);
pub(super) const R64: AdmissibleR = admissible(fp64x2::P as u128, 4);
pub(super) const R_GOLDILOCKS: AdmissibleR = admissible(goldilocks2::P as u128, 4);
pub(super) const R_WEIER61: AdmissibleR = admissible(fp61x2::P as u128, 1);

/// Parameter sampling reduces coefficients; point hashing instead rejects
/// noncanonical encodings. The two operations define different distributions.
fn reduce_parameter<F: Packed>(v: u128) -> F {
    F::reduce(v & (u128::MAX >> (128 - F::BITS)))
}

fn parameter<F: Packed>(tag: &[u8], seed: &[u8; 32], index: u32) -> F {
    reduce_parameter(half(tag, seed, index))
}

/// These three packings put the base coefficient below the extension one.
fn outside_base_field<F: Packed>(j: F) -> bool {
    j.pack() >> (F::BITS / 2) != 0
}

/// j for y^2 = x^3 + a2*x^2 + a4*x on a nonsingular model.
fn cubic_j<F: Packed>(a2: F, a4: F) -> F {
    let a = a2.square();
    let t = a - F::reduce(3) * a4;
    F::reduce(256) * t.square() * t * (a4.square() * (a - F::reduce(4) * a4)).inv()
}

fn twisted_candidate<F: Packed>(
    tag: &[u8],
    seed: &[u8; 32],
    index: u32,
) -> Option<twisted::Curve<F>> {
    let c = twisted::Curve::new(parameter::<F>(tag, seed, index))?;
    let half = (F::ONE + F::ONE).inv();
    let a2 = (c.d - F::ONE) * half;
    let a4 = ((c.d + F::ONE) * half.square()).square();
    outside_base_field(cubic_j(a2, a4)).then_some(c)
}

fn weier_curve(b: fp61x2::Fq) -> Option<weier61x2::Curve> {
    let c = weier61x2::Curve::new(b)?;
    let j = fp61x2::Fq::reduce(6912) * (fp61x2::Fq::reduce(4) - b.square()).inv();
    outside_base_field(j).then_some(c)
}

/// A twisted family: -x^2 + y^2 = 1 + d x^2 y^2 over F_{p^2} of order 4r,
/// represented in `E/<T>`, where the hashed point is cleared by 2 and the
/// l = 8 marker is a point of order 4.
macro_rules! twisted_policy {
    ($name:ident, $module:ident, $tag:ident, $r:ident, $candidate:ident, $verify:ident) => {
        pub struct $name;

        impl Criteria for $name {
            type Curve = $module::Curve;
            const TAG: &'static [u8] = $tag;
            const R: AdmissibleR = $r;
            fn candidate(seed: &[u8; 32], j: u32) -> Option<Self::Curve> {
                twisted_candidate(Self::TAG, seed, j)
            }
        }

        impl Policy<16> for $name {
            const CLEAR: u128 = 2;
            const ENTRY_PER_INDEX: bool = false;

            fn marker(c: &Self::Curve, p: &$module::Affine, l: u128) -> Option<bool> {
                eight_torsion(Self::CLEAR, c, p, l)
            }
        }

        pub fn $candidate(seed: &[u8; 32], j: u32) -> Option<$module::Curve> {
            $name::candidate(seed, j)
        }

        pub fn $verify(seed: &[u8; 32], cert: &Certificate) -> Result<$module::Curve, Error> {
            crate::curvegen::select::verify::<16, $name>(seed, cert)
        }
    };
}

twisted_policy!(
    Twisted61x2,
    twisted61x2,
    TAG_TWISTED61X2,
    R61,
    twisted61x2_candidate,
    verify_twisted61x2
);
twisted_policy!(
    Twisted64x2,
    twisted64x2,
    TAG_TWISTED64X2,
    R64,
    twisted64x2_candidate,
    verify_twisted64x2
);
twisted_policy!(
    TwistedGoldilocks2,
    twisted_goldilocks2,
    TAG_TWISTED_GOLDILOCKS2,
    R_GOLDILOCKS,
    twisted_goldilocks2_candidate,
    verify_twisted_goldilocks2
);

/// Edwards x^2 + y^2 = 1 + d x^2 y^2 over F_{p^2}, p = 2^61 - 1, #E = 4r.
pub struct Edwards61x2;

impl Criteria for Edwards61x2 {
    type Curve = edwards61x2::Curve;
    const TAG: &'static [u8] = TAG_EDWARDS61X2;
    const R: AdmissibleR = R61;
    fn candidate(seed: &[u8; 32], j: u32) -> Option<Self::Curve> {
        let c = edwards61x2::Curve::new(parameter(Self::TAG, seed, j))?;
        outside_base_field(cubic_j(c.a2, c.a4)).then_some(c)
    }
}

impl Policy<16> for Edwards61x2 {
    const CLEAR: u128 = 4;
    const ENTRY_PER_INDEX: bool = false;

    fn marker(c: &Self::Curve, p: &edwards61x2::Affine, l: u128) -> Option<bool> {
        eight_torsion(Self::CLEAR, c, p, l)
    }
}

/// y^2 = x^3 - 3x + b over F_{p^2}, p = 2^61 - 1, #E = r. (q = p^2 is
/// composite, so the anomalous exclusion never decides a verdict here.)
pub struct Weier61x2;

impl Criteria for Weier61x2 {
    type Curve = weier61x2::Curve;
    const TAG: &'static [u8] = TAG_WEIER61X2;
    const R: AdmissibleR = R_WEIER61;
    fn candidate(seed: &[u8; 32], j: u32) -> Option<Self::Curve> {
        weier_curve(parameter(Self::TAG, seed, j))
    }
}

impl Policy<16> for Weier61x2 {
    const CLEAR: u128 = 1;
    const ENTRY_PER_INDEX: bool = false;

    fn marker(_: &Self::Curve, p: &weier61x2::Affine, l: u128) -> Option<bool> {
        two_torsion(p, l)
    }
}

pub fn edwards61x2_candidate(seed: &[u8; 32], j: u32) -> Option<edwards61x2::Curve> {
    Edwards61x2::candidate(seed, j)
}

pub fn weier61x2_candidate(seed: &[u8; 32], j: u32) -> Option<weier61x2::Curve> {
    Weier61x2::candidate(seed, j)
}

pub fn verify_edwards61x2(
    seed: &[u8; 32],
    cert: &Certificate,
) -> Result<edwards61x2::Curve, Error> {
    crate::curvegen::select::verify::<16, Edwards61x2>(seed, cert)
}

pub fn verify_weier61x2(seed: &[u8; 32], cert: &Certificate) -> Result<weier61x2::Curve, Error> {
    crate::curvegen::select::verify::<16, Weier61x2>(seed, cert)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::{Field, fp61x2, fp64x2, goldilocks2};

    #[test]
    fn quadratic_bounds_are_exact() {
        for (p, r) in [
            (fp61x2::P as u128, R61),
            (fp64x2::P as u128, R64),
            (goldilocks2::P as u128, R_GOLDILOCKS),
        ] {
            // Both Hasse endpoints are multiples of 4, so division is exact.
            assert_eq!(r.q, p * p);
            assert_eq!(r.lo, (p - 1) * (p - 1) / 4);
            assert_eq!(r.hi, (p + 1) * (p + 1) / 4);
            assert!(r.hasse_contains((p - 1) * (p - 1)) && r.hasse_contains((p + 1) * (p + 1)));
            assert!(!r.hasse_contains((p - 1) * (p - 1) - 1));
            assert!(r.admits(r.lo) && r.admits(r.hi));
            assert!(!r.admits(r.lo - 1) && !r.admits(r.hi + 1) && !r.admits(u128::MAX));
            // r = q is out of the window at cofactor 4
            assert!(r.hi < r.q);
        }
        assert_eq!(R_WEIER61.lo, (fp61x2::P as u128 - 1).pow(2));
        assert_eq!(R_WEIER61.hi, (fp61x2::P as u128 + 1).pow(2));
        // The window contains q at cofactor 1, though q = p^2 is composite.
        assert!(R_WEIER61.admits(R_WEIER61.q));
    }

    fn reduction<F: Packed>(p: u128) {
        let w = F::BITS / 2;
        let mask = u128::MAX >> (128 - F::BITS);
        for v in [0, 1, p, p << w, mask, u128::MAX] {
            let masked = v & mask;
            let a = masked & (u128::MAX >> (128 - w));
            let b = masked >> w;
            let want = (a % p) | ((b % p) << w);
            assert_eq!(reduce_parameter::<F>(v).pack(), want);
        }
    }

    #[test]
    fn candidate_reduction_is_coefficient_wise() {
        reduction::<fp61x2::Fq>(fp61x2::P as u128);
        reduction::<fp64x2::Fq>(fp64x2::P as u128);
        reduction::<goldilocks2::Fq>(goldilocks2::P as u128);
    }

    fn pow<F: Field>(mut a: F, mut n: u128) -> F {
        let mut r = F::ONE;
        while n != 0 {
            if n & 1 != 0 {
                r = r * a;
            }
            a = a.square();
            n >>= 1;
        }
        r
    }

    fn subfield<F: Packed>(p: u128) {
        for v in [0, 1, 2, 1 << (F::BITS / 2), 1 | 1 << (F::BITS / 2)] {
            let x = F::reduce(v);
            assert_eq!(outside_base_field(x), !pow(x, p).equals(x));
        }
    }

    #[test]
    fn structural_filter_uses_j_not_the_curve_coefficient() {
        subfield::<fp61x2::Fq>(fp61x2::P as u128);
        subfield::<fp64x2::Fq>(fp64x2::P as u128);
        subfield::<goldilocks2::Fq>(goldilocks2::P as u128);
        type F = fp61x2::Fq;
        let i = F::unpack(1 << 61).unwrap();
        assert!(outside_base_field(i));
        // b=i is outside F_p, but b^2=-1 and j=6912/5 are in F_p.
        assert!(weier61x2::Curve::new(i).is_some());
        assert!(weier_curve(i).is_none());
        assert!(weier_curve(F::ONE + i).is_some());
    }
}
