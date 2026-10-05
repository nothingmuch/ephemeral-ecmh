//! `binary::Pornin`'s map over GF(2^122) = GF(2^61)\[u\]/(u^2 + u + 1),
//! a = u, for both families: crrl's `gls254::map_to_curve` (GF(2^127)\[u\],
//! a = u) at half the size, with this base modulus's constants.
//!
//! Tr_122(a0 + a1 u) = Tr_61(a1), and Tr_61(v) = v_0 on a canonical v for
//! f = z^61 + z^23 + z^15 + z^5 + 1 (`field::gf2_122`), so, as in crrl:
//!
//! - k = z^2 lies in GF(2^61), so Tr(k) = 0.
//! - Tr(c) = c1_0, bit 64 of the canonical c.
//! - Tr(c^2/z^2) = Tr(c/z) = Tr_61(c1/z) = c1_1, bit 65: c1/z is c1 >> 1
//!   plus c1_0 z^-1 = c1_0 (z^60 + z^22 + z^14 + z^4), of trace 0.
//!
//! So the map sets bit 64 and clears bit 65, and c^2/k is (c/z)^2, both
//! halves shifted and one squaring. m = 122 is even and Tr(1) = 0, so
//! the trace does not tell w from w + 1; adding 1 flips bit 0 of w0,
//! which picks the root as in crrl. Tr_61(w0) is that bit on any
//! representative. The input keeps 120 bits of c, and the sign at bit 127
//! (sage/pornin_map.sage checks each step).

use super::{M122, M122Gls};
use crate::curve::binary::Pornin;
use crate::field::gf2_122::{Gf, MASK122, from_u128, to_u128};

/// f with its z^61 term: w/z = (w + w_0 f) >> 1 in GF(2^61).
const F61: u64 = 1 << 61 | 1 << 23 | 1 << 15 | 1 << 5 | 1;

/// w/z for a canonical w in GF(2^61).
#[inline(always)]
fn div_z(w: u64) -> u64 {
    (w ^ (w & 1).wrapping_neg() & F61) >> 1
}

macro_rules! pornin {
    ($($m:ty),*) => {$(
        impl Pornin for $m {
            #[inline(always)]
            fn k() -> Gf {
                Gf::w64le(4, 0)
            }
            /// c1_0 = 1 and c1_1 = 0: bits 64 and 65.
            #[inline(always)]
            fn force(v: u128) -> Gf {
                from_u128((v & MASK122 & !(1 << 65)) | 1 << 64)
            }
            /// (c/z)^2.
            #[inline(always)]
            fn sq_over_k(c: Gf) -> Gf {
                let v = to_u128(c);
                Gf::w64le(div_z(v as u64), div_z((v >> 64) as u64)).square()
            }
            #[inline(always)]
            fn root_bit(w: Gf) -> u32 {
                w.parts().0.trace()
            }
        }
    )*};
}
pornin!(M122, M122Gls);

#[cfg(test)]
mod tests {
    mod dense {
        crate::curve::binary::map_suite!(
            crate::curve::binary122::M122,
            crate::curve::binary122::tests::dense::curve()
        );

        #[test]
        fn matches_sage() {
            // sage/pornin_map.sage 122: (input, x, y) on B = 0x1f3c8a712e94d0b61c47a389e5f26d1b
            const VECTORS: [(u128, u128, u128); 8] = [
                (
                    0x9e3779b97f4a7c15f39cc0605cedc835,
                    0x8e9cc9487f7c57f02296ec81c609cc7,
                    0x985d6da184f3ae902b1c395d96a8545,
                ),
                (
                    0x3c6ef372fe94f82be73980c0b9db906a,
                    0xae0170d959b4290cba06a32cb31632,
                    0x1b137984e514c8820e6bf8a0475f0c02,
                ),
                (
                    0xdaa66d2c7ddf7441dad6412116c9589f,
                    0x113c1b834702ddf0b38178005660b8b,
                    0x1e50204eafaee92f0ea75301304ecbb9,
                ),
                (
                    0x78dde6e5fd29f057ce73018173b720d4,
                    0xf108b4ef38f311d061592a490fb9a9e,
                    0x6c47ad54ac244ff0cdd2b274f8c1ba9,
                ),
                (
                    0x1715609f7c746c6dc20fc1e1d0a4e909,
                    0xf27164ce1925cb7028ab0f89886abb3,
                    0x19d357c5406c3af31cffa0fc61031d57,
                ),
                (
                    0xb54cda58fbbee883b5ac82422d92b13e,
                    0x1feb612ce1b9f611052b9e0983425142,
                    0x42a0d81e69cb896031cb73981d303f5,
                ),
                (
                    0x538454127b096499a94942a28a807973,
                    0x1e5ea56214d7b31306053d667247c99d,
                    0x19d84986bd6e71d108fd39af5c151b85,
                ),
                (
                    0xf1bbcdcbfa53e0af9ce60302e76e41a8,
                    0x10aaecd031dd3add1979ca96e3337400,
                    0x1cf3eca5aaa04d8309cbf147eeb9c99f,
                ),
            ];
            let c = Curve::new(crate::field::gf2_122::from_u128(
                0x5f3c_8a71_2e94_d0b6_1c47_a389_e5f2_6d1b,
            ));
            crate::curve::binary::check_map_vectors(&c, &VECTORS);
        }
    }

    mod gls {
        crate::curve::binary::map_suite!(
            crate::curve::binary122::M122Gls,
            crate::curve::binary122::tests::gls::curve()
        );

        #[test]
        fn matches_sage() {
            // sage/pornin_map.sage 122-gls: (input, x, y) on B = 0xa7ed0c2882f395b
            const VECTORS: [(u128, u128, u128); 8] = [
                (
                    0x9e3779b97f4a7c15f39cc0605cedc835,
                    0x162108390d0424590e1cf5da500ba997,
                    0x1ec75edefed63d2f050ab20dadf2e9b6,
                ),
                (
                    0x3c6ef372fe94f82be73980c0b9db906a,
                    0x12d01d16d2a68d07179b19b8626744f5,
                    0x13bee55cf9cc30c61e9b5e8284126e12,
                ),
                (
                    0xdaa66d2c7ddf7441dad6412116c9589f,
                    0x7847bb85e4f541d17631f02e951e130,
                    0x1734792f0d91706a03041268075c13e4,
                ),
                (
                    0x78dde6e5fd29f057ce73018173b720d4,
                    0x63a3bcbe085515703b225e701e8355c,
                    0x18f31e0a4442eb3d02a1c959bcd82da6,
                ),
                (
                    0x1715609f7c746c6dc20fc1e1d0a4e909,
                    0x13bfdbe4012b2e5f066ff81f41094561,
                    0x7126da0e30dd49306ee35e212a5d5cd,
                ),
                (
                    0xb54cda58fbbee883b5ac82422d92b13e,
                    0x8fddc1f7cbf218302ca19a0b70ebe28,
                    0x1ae720d945dbab7700e8403d97d934fa,
                ),
                (
                    0x538454127b096499a94942a28a807973,
                    0x76f9201ba7ba95f13c6b22228751c47,
                    0xc516bb3a5a5baa00700f566a0d638bb,
                ),
                (
                    0xf1bbcdcbfa53e0af9ce60302e76e41a8,
                    0x5d7c0e60a49a9110ad1f97186e6f519,
                    0x1f5a8b63adf56b7d06ba557cf6d5760d,
                ),
            ];
            let beta = crate::field::gf2_122::gf2_61::from_u64(0x1c47_a389_e5f2_6d1b);
            crate::curve::binary::check_map_vectors(&Curve::new(beta.xsquare(2)), &VECTORS);
        }
    }
}
