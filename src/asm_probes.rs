//! Assembly probes: every field operation and every family's principal
//! group operations, each as its own non-inlined function, so the
//! compiler's output for each can be read and compared across targets
//! (`nix run .#asm-report`).
//!
//! In real callers these operations inline, so a probe shows the
//! arithmetic itself plus a call boundary: its argument and return moves
//! do not occur in the inlined code. Each probe is
//! `asm_probes::<field or family>::<op>`, which `asm-report` finds by its
//! (demangled) path.
//!
//! Only built with the `asm-probes` feature.

use crate::group::{Accumulate, Decode, Group, Negate};

/// Probes for one field: `mul`, `square`, `add`, `sub`, `neg`, `invert` and
/// `sqrt`, each a function of the field's own type, with the operations
/// that vary between fields' APIs passed in as closures.
macro_rules! field {
    ($m:ident, $F:ty, neg: $neg:expr, sub: $sub:expr, square: $square:expr, invert: $invert:expr, sqrt: $sqrt:expr $(,)?) => {
        pub mod $m {
            #[allow(unused_imports)]
            use super::*;
            type F = $F;
            #[inline(never)]
            pub fn mul(a: F, b: F) -> F {
                a * b
            }
            #[inline(never)]
            pub fn square(a: F) -> F {
                ($square)(a)
            }
            #[inline(never)]
            pub fn add(a: F, b: F) -> F {
                a + b
            }
            #[inline(never)]
            pub fn sub(a: F, b: F) -> F {
                ($sub)(a, b)
            }
            #[inline(never)]
            pub fn neg(a: F) -> F {
                ($neg)(a)
            }
            #[inline(never)]
            pub fn invert(a: F) -> F {
                ($invert)(a)
            }
            #[inline(never)]
            pub fn sqrt(a: F) -> Option<F> {
                ($sqrt)(a)
            }
        }
    };
}

/// `field!` for odd characteristic.
macro_rules! prime_field {
    ($m:ident, $F:ty, square: $square:expr, invert: $invert:expr, sqrt: $sqrt:expr $(,)?) => {
        field!($m, $F, neg: |a: $F| -a, sub: |a: $F, b: $F| a - b,
            square: $square, invert: $invert, sqrt: $sqrt);
    };
}

/// `field!` for characteristic 2: negation is the identity, subtraction
/// addition.
macro_rules! binary_field {
    ($m:ident, $F:ty) => {
        field!($m, $F, neg: |a: $F| a, sub: |a: $F, b: $F| a + b,
            square: |a: $F| a.square(), invert: |a: $F| a.invert(), sqrt: |a: $F| Some(a.sqrt()));
    };
}

pub mod field {
    binary_field!(gf2_127, crate::field::gf2_127::Gf);
    binary_field!(gf2_109, crate::field::gf2_109::Gf);
    binary_field!(gf2_122, crate::field::gf2_122::Gf);
    binary_field!(gf2_61, crate::field::gf2_122::Gf61);

    pub mod halftrace {
        #[inline(never)]
        pub fn gf2_127(a: crate::field::gf2_127::Gf) -> crate::field::gf2_127::Gf {
            crate::field::gf2_127::halftrace8(a)
        }
        #[inline(never)]
        pub fn gf2_109(a: crate::field::gf2_109::Gf) -> crate::field::gf2_109::Gf {
            crate::field::gf2_109::halftrace8(a)
        }
        #[inline(never)]
        pub fn gf2_122(a: crate::field::gf2_122::Gf) -> crate::field::gf2_122::Gf {
            a.qsolve()
        }
    }

    pub mod prime {
        use p3_field::PrimeCharacteristicRing;

        prime_field!(fp127, crate::field::fp127::Fp, square: |a: crate::field::fp127::Fp| a.square(),
            invert: |a: crate::field::fp127::Fp| a.invert(), sqrt: |a: crate::field::fp127::Fp| a.sqrt());
        prime_field!(fp107, crate::field::fp107::Fp, square: |a: crate::field::fp107::Fp| a.square(),
            invert: |a: crate::field::fp107::Fp| a.invert(), sqrt: |a: crate::field::fp107::Fp| a.sqrt());
        prime_field!(fp128, crate::field::fp128::Fp, square: |a: crate::field::fp128::Fp| a.square(),
            invert: |a: crate::field::fp128::Fp| a.invert(), sqrt: |a: crate::field::fp128::Fp| a.sqrt());
        prime_field!(fp61, crate::field::fp61x2::Fp, square: |a: crate::field::fp61x2::Fp| a.square(),
            invert: |a: crate::field::fp61x2::Fp| a.invert(), sqrt: |a: crate::field::fp61x2::Fp| a.sqrt());
        prime_field!(fp61x2, crate::field::fp61x2::Fq, square: |a: crate::field::fp61x2::Fq| a.square(),
            invert: |a: crate::field::fp61x2::Fq| a.invert(), sqrt: |a: crate::field::fp61x2::Fq| a.sqrt());
        prime_field!(fp64, crate::field::fp64x2::Fp, square: |a: crate::field::fp64x2::Fp| a.square(),
            invert: |a: crate::field::fp64x2::Fp| a.invert(), sqrt: |a: crate::field::fp64x2::Fp| a.sqrt());
        prime_field!(fp64x2, crate::field::fp64x2::Fq, square: |a: crate::field::fp64x2::Fq| a.square(),
            invert: |a: crate::field::fp64x2::Fq| a.invert(), sqrt: |a: crate::field::fp64x2::Fq| a.sqrt());
        prime_field!(goldilocks, crate::field::goldilocks2::Fp, square: |a: crate::field::goldilocks2::Fp| a.square(),
            invert: crate::field::goldilocks2::invert_base, sqrt: crate::field::goldilocks2::sqrt_base);
        prime_field!(goldilocks2, crate::field::goldilocks2::Fq, square: |a: crate::field::goldilocks2::Fq| a.square(),
            invert: crate::field::goldilocks2::invert, sqrt: crate::field::goldilocks2::sqrt);
    }
}

/// Probes for one family through the `group` traits, as RIBLT drives it:
/// `add` (accumulator += addend), `sub`, `equals` (the purity check of a
/// cell against an addend), `prepare` and `decode`, called
/// through the traits, as some curves have inherent methods of the same
/// names.
macro_rules! family {
    ($m:ident, $G:ty) => {
        pub mod $m {
            use super::*;
            type G = $G;
            #[inline(never)]
            pub fn add(
                g: &G,
                p: &<G as Group>::Point,
                a: &<G as Accumulate>::Addend,
            ) -> <G as Group>::Point {
                <G as Accumulate>::add(g, p, a)
            }
            #[inline(never)]
            pub fn sub(
                g: &G,
                p: &<G as Group>::Point,
                a: &<G as Accumulate>::Addend,
            ) -> <G as Group>::Point {
                <G as Accumulate>::add(g, p, &<G as Negate>::neg_addend(g, a))
            }
            #[inline(never)]
            pub fn equals(g: &G, p: &<G as Group>::Point, a: &<G as Accumulate>::Addend) -> bool {
                <G as Negate>::equals_addend(g, p, a)
            }
            #[inline(never)]
            pub fn prepare(g: &G, a: &<G as Group>::Affine) -> <G as Accumulate>::Addend {
                <G as Accumulate>::prepare(g, a)
            }
            #[inline(never)]
            pub fn decode(
                g: &G,
                e: &<G as crate::group::Encode>::Encoding,
            ) -> Option<<G as Group>::Affine> {
                <G as Decode>::decode(g, e)
            }
        }
    };
}

pub mod group {
    use super::*;

    family!(binary127, crate::curve::binary127::Curve);
    family!(
        binary127_lambda,
        crate::curve::binary::lambda::Curve<crate::curve::binary127::M127>
    );
    family!(
        binary127_w,
        crate::curve::binary::wcodec::Curve<crate::curve::binary127::M127>
    );
    family!(binary109, crate::curve::binary109::Curve);
    family!(
        binary109_lambda,
        crate::curve::binary::lambda::Curve<crate::curve::binary109::M109>
    );
    family!(
        binary109_w,
        crate::curve::binary::wcodec::Curve<crate::curve::binary109::M109>
    );
    family!(binary122, crate::curve::binary122::Dense);
    family!(
        binary122_lambda,
        crate::curve::binary::lambda::Curve<crate::curve::binary122::M122>
    );
    family!(binary122_gls, crate::curve::binary122::Gls);
    family!(
        binary122_gls_lambda,
        crate::curve::binary::lambda::Curve<crate::curve::binary122::M122Gls>
    );
    family!(
        binary122_w,
        crate::curve::binary::wcodec::Curve<crate::curve::binary122::M122>
    );
    family!(
        binary122_gls_w,
        crate::curve::binary::wcodec::Curve<crate::curve::binary122::M122Gls>
    );
    family!(edwards127, crate::curve::edwards127::Curve);
    family!(edwards107, crate::curve::edwards107::Curve);
    family!(twisted128, crate::curve::twisted128::Curve);
    family!(weier127, crate::curve::weier127::OddCurve);
    family!(
        weier127_jacobian,
        crate::curve::weier::jacobian::Curve<crate::field::fp127::Fp>
    );
    family!(weier107, crate::curve::weier107::OddCurve);
    family!(
        weier107_jacobian,
        crate::curve::weier::jacobian::Curve<crate::field::fp107::Fp>
    );
    family!(
        binary127_u,
        crate::curve::binary::unscaled::Curve<crate::curve::binary127::M127>
    );
    family!(
        binary109_u,
        crate::curve::binary::unscaled::Curve<crate::curve::binary109::M109>
    );
    family!(
        binary122_u,
        crate::curve::binary::unscaled::Curve<crate::curve::binary122::M122>
    );
    family!(
        binary122_gls_u,
        crate::curve::binary::unscaled::Curve<crate::curve::binary122::M122Gls>
    );
    family!(edwards61x2, crate::curve::edwards61x2::Curve);
    family!(weier61x2, crate::curve::weier61x2::OddCurve);
    family!(
        weier61x2_jacobian,
        crate::curve::weier::jacobian::Curve<crate::field::fp61x2::Fq>
    );
    family!(twisted61x2, crate::curve::twisted61x2::Curve);
    family!(twisted64x2, crate::curve::twisted64x2::Curve);
    family!(
        twisted_goldilocks2,
        crate::curve::twisted_goldilocks2::Curve
    );
}
