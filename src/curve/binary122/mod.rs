//! `binary`'s curves over GF(2^122) = GF(2^61)\[u\]/(u^2 + u + 1): GLS254's
//! field shape at half the size, in two families.
//!
//! - Field: `field::gf2_122`, the tower over F_2\[z\]/(z^61 + z^23 +
//!   z^15 + z^5 + 1).
//! - Encoding: 16 bytes, `to_u128(x) | Tr(y) << 127`. Tr(x) = 1 is bit 64
//!   (Tr_122(x) = Tr_61(x1), which is bit 0 of a canonical x1), so every
//!   point but O has it set.
//!
//! - `map`: `binary::Pornin`'s map for both families, crrl's GLS254
//!   map at half the size.
//!
//! The families differ only in the type B is drawn from:
//!
//! - Dense ([`M122`], `binary.122`): B random in GF(2^122), so
//!   beta = B^(1/4) is a full element and each m_beta costs a
//!   multiplication, as on `binary127`'s random curves. rho is about 2^60.3.
//! - GLS ([`M122Gls`], `binary.122-gls`): B = beta^4 with beta random in
//!   GF(2^61). With a = u that is the quadratic twist of a curve over
//!   GF(2^61), of order (q - 1)^2 + t^2 (q = 2^61), so #E = 2r still
//!   happens. m_beta is a GF(2^122) x GF(2^61) product, 2 M61 against
//!   M = 3 M61. The GLS endomorphism psi exists, and gains an attacker
//!   sqrt 2: rho about 2^59.8.
//!
//! # Why a = u
//!
//! m = 122 is even, so Tr(1) = 0 and a = 1 would be isomorphic to a = 0,
//! whose orders are 4 * odd at best. Tr(u) = Tr_61(1) = 1 (see
//! `Gf::trace`), so a = u gives #E = 2 * odd, and the doubles, E\[r\] = 2E,
//! are the points with Tr(x) = Tr(a) = 1. Pornin's formulas then carry
//! a^2 = u^2 = u + 1 and 1 + a = u^2, both multiplications by u^2:
//! (a0 + a1 u)(1 + u) = (a0 + a1) + a0 u, one addition ([`mul_u2`]).
//!
//! # Costs per addition
//!
//! M = GF(2^122) multiplication = 3 M61; S = squaring = 2 S61.
//!
//! | family | add (extended + extended) | m_beta | total |
//! |---|---|---|---|
//! | `binary.122` | 8M + 2S + 2 m_beta | 1M | 10M + 2S |
//! | `binary.122-gls` | 8M + 2S + 2 m_beta | 2 M61 = 2/3 M | 9.33M + 2S (28 M61 + 4 S61) |
//! | `binary-lambda.122`, `-gls` | 8M + 2S (`binary::lambda`) | none | 8M + 2S |
//!
//! `binary127`'s add is also 8M + 2S + 2 m_beta, with its M being 4
//! carry-less products and a trinomial fold (see `field::gf2_122`'s
//! Backends). The a^2 term (E = a^2 T1 T2) costs one addition here; for
//! Pornin's GLS254 and for a = 1 it requires no operation.

use crate::curve::binary::{self, Constant};
use crate::field::gf2_122::{Gf, Gf61};

mod map;

/// x u^2 = x (u + 1): (a0 + a1 u)(1 + u) = (a0 + a1) + a0 u.
#[inline(always)]
pub fn mul_u2(x: Gf) -> Gf {
    let [a0, a1] = x.limbs();
    Gf::w64le(a0 ^ a1, a0)
}

/// A constant in GF(2^61): its product with an element of GF(2^122) is
/// 2 M61.
impl Constant<Gf> for Gf61 {
    #[inline(always)]
    fn gf(self) -> Gf {
        Gf::from_base(self)
    }
    #[inline(always)]
    fn mul(self, x: Gf) -> Gf {
        x.mul_base(self)
    }
    fn sqrt(self) -> Self {
        Gf61::sqrt(self)
    }
    fn recip(self) -> Self {
        crate::field::batch::Invert::inv(self)
    }
}

/// a = u over GF(2^122), B drawn from `$k`.
macro_rules! model {
    ($(#[$doc:meta])* $name:ident, $k:ty) => {
        $(#[$doc])*
        #[derive(Clone, Copy, Debug)]
        pub struct $name;

        impl binary::Model for $name {
            type F = Gf;
            type K = $k;
            const A: Gf = Gf::U;
            #[inline(always)]
            fn mul_a(x: Gf) -> Gf {
                x.mul_u()
            }
            /// u^2 = u + 1
            #[inline(always)]
            fn mul_a2(x: Gf) -> Gf {
                mul_u2(x)
            }
            /// 1 + u = u^2
            #[inline(always)]
            fn mul_1a(x: Gf) -> Gf {
                mul_u2(x)
            }
            const SIGN: u32 = 127;
            type Bytes = [u8; 16];
        }
    };
}
model!(
    /// `binary.122`: B dense in GF(2^122).
    M122,
    Gf
);
model!(
    /// `binary.122-gls`: B = beta^4, beta in GF(2^61).
    M122Gls,
    Gf61
);

pub type Dense = binary::Curve<M122>;
pub type Gls = binary::Curve<M122Gls>;

pub fn candidate(c: u128) -> (Gf, u32) {
    binary::candidate::<M122>(c)
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::field::gf2_122::gf2_61::{self, MASK61};
    use crate::field::gf2_122::{MASK122, eq, from_u128, to_u128};
    use proptest::prelude::*;

    pub(crate) mod dense {
        crate::curve::binary::tests::suite!(gf2_122, crate::curve::binary122::M122);
    }

    pub(crate) mod gls {
        use crate::field::gf2_122::gf2_61::{self, MASK61};
        crate::curve::binary::tests::suite!(
            gf2_122,
            crate::curve::binary122::M122Gls,
            any::<u64>()
                .prop_map(|v| v & MASK61)
                .prop_filter("beta != 0", |&v| v != 0)
                .prop_map(|v| Curve::new(gf2_61::from_u64(v).xsquare(2)))
        );
    }

    #[test]
    fn mul_u2_is_mul_by_u_squared() {
        let x = from_u128(0x1234_5678_9abc_def0_0fed_cba9_8765_4321 & MASK122);
        assert!(eq(mul_u2(x), x * Gf::U * Gf::U));
        assert!(eq(mul_u2(x), x.mul_u().mul_u()));
    }

    #[test]
    fn constants_agree() {
        let beta = gf2_61::from_u64(0x0123_4567_89ab_cdef & MASK61);
        let (g, d) = (
            Gls::new(beta.xsquare(2)),
            Dense::new(Gf::from_base(beta).xsquare(2)),
        );
        assert!(eq(g.beta.gf(), d.beta) && eq(g.b.gf(), d.b) && eq(g.big_b.gf(), d.big_b));
        let x = from_u128(0xfeed_face_cafe_beef_0bad_f00d_dead_c0de & MASK122);
        assert!(eq(g.beta.mul(x), d.beta.mul(x)));
    }

    proptest! {
        #[test]
        fn candidate_sets_bit_64(v in any::<u128>()) {
            // Tr(x) = x_64 here, so fixing the trace is setting bit 64
            let (x, sign) = candidate(v);
            prop_assert_eq!(to_u128(x), (v | 1 << 64) & MASK122);
            prop_assert_eq!(sign as u128, v >> 127);
        }

        #[test]
        fn encodings_set_bit_64((_c, ps) in dense::curve_and_points(1)) {
            prop_assert_eq!(ps[0].encode()[8] & 1, 1);
        }

        #[test]
        fn gls_matches_dense_on_the_same_curve((c, ps) in gls::curve_and_points(2)) {
            let d = Dense::new(c.big_b.gf());
            let dense = |p: &binary::Affine<M122Gls>| binary::Affine::<M122> { x: p.x, y: p.y };
            let (p, q) = (ps[0], ps[1]);
            let a = c.add(&c.from_affine(&p), &c.from_affine(&q));
            let b = d.add(&d.from_affine(&dense(&p)), &d.from_affine(&dense(&q)));
            prop_assert_eq!(dense(&c.to_affine(&a)), d.to_affine(&b));
        }
    }
}
