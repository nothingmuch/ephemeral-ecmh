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

    pub mod halftrace {
        #[inline(never)]
        pub fn gf2_127(a: crate::field::gf2_127::Gf) -> crate::field::gf2_127::Gf {
            crate::field::gf2_127::halftrace8(a)
        }
    }

    pub mod prime {
        prime_field!(fp127, crate::field::fp127::Fp, square: |a: crate::field::fp127::Fp| a.square(),
            invert: |a: crate::field::fp127::Fp| a.invert(), sqrt: |a: crate::field::fp127::Fp| a.sqrt());
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
    family!(
        binary127_u,
        crate::curve::binary::unscaled::Curve<crate::curve::binary127::M127>
    );
}
