//! `binary::Pornin`'s constants for f = z^109 + z^5 + z^4 + z^2 + 1,
//! a = 1: the GF(2^127) constants carry over, with two parities where
//! that modulus has single bits (sage/pornin_map.sage checks each step).
//!
//! The map needs Tr(k) = 0 and, on the forced c, Tr(c) = Tr(a) = 1 and
//! Tr(c^2/k) = 0. Take k = z^2 again:
//!
//! - Tr(z^2) = Tr(z) = 0, since Tr(z^i) = 1 for i < 109 exactly at
//!   i = 0, 105, 107 (`Gf::trace`). So Tr(c) = c_0 + c_105 + c_107.
//! - Tr(c^2/z^2) = Tr(c/z). With z^-1 = z^108 + z^4 + z^3 + z (from
//!   z (z^108 + z^4 + z^3 + z) = f + 1), c/z is c >> 1 plus c_0 z^-1,
//!   and Tr(z^-1) = 0, so Tr(c/z) = c_1 + c_106 + c_108.
//!
//! Bit 0 enters only the first form and bit 1 only the second, so setting
//! c_0 = 1 + c_105 + c_107 and c_1 = c_106 + c_108 forces both, as
//! setting c_0 = 1 and c_1 = 0 does for z^127 + z^63 + 1. c^2/k is then
//! (c/z)^2, a shift on the canonical c and a squaring. 109 is odd, so
//! Tr(1) = 1 and the trace tells w from w + 1, as over GF(2^127). The input
//! keeps 107 bits of c, and the sign at bit 109 (`M109::SIGN`).

use super::M109;
use crate::curve::binary::Pornin;
use crate::field::gf2_109::{Gf, MASK109, from_u128, to_u128};

/// f, with its z^109 term: c/z = (c + c_0 f) >> 1.
const F: u128 = 1 << 109 | 0b11_0101;

impl Pornin for M109 {
    #[inline(always)]
    fn k() -> Gf {
        Gf::w64le(4, 0)
    }
    /// c_0 = 1 + c_105 + c_107 and c_1 = c_106 + c_108.
    #[inline(always)]
    fn force(v: u128) -> Gf {
        let v = v & MASK109;
        let c0 = !(v >> 105 ^ v >> 107) & 1;
        let c1 = (v >> 106 ^ v >> 108) & 1;
        from_u128(v & !3 | c0 | c1 << 1)
    }
    /// (c/z)^2.
    #[inline(always)]
    fn sq_over_k(c: Gf) -> Gf {
        let v = to_u128(c);
        from_u128((v ^ (v & 1).wrapping_neg() & F) >> 1).square()
    }
    #[inline(always)]
    fn root_bit(w: Gf) -> u32 {
        w.trace()
    }
}

#[cfg(test)]
mod tests {
    crate::curve::binary::map_suite!(
        crate::curve::binary109::M109,
        crate::curve::binary109::tests::curve()
    );

    #[test]
    fn matches_sage() {
        // sage/pornin_map.sage 109: (input, x, y) on B = 0xa712e94d0b61c47a389e5f26d1b
        const VECTORS: [(u128, u128, u128); 8] = [
            (
                0x9e3779b97f4a7c15f39cc0605cedc835,
                0x1ef58619c78ff98dba7b22c91cc5,
                0x63d4cf2840369e32e1eec36a019,
            ),
            (
                0x3c6ef372fe94f82be73980c0b9db906a,
                0x1f79527c9b3f3b2753085273eb83,
                0xe7ae3ddb05b2c96c74d19b3559a,
            ),
            (
                0xdaa66d2c7ddf7441dad6412116c9589f,
                0x72f011b521296c8c5c2eeffa090,
                0x8871382e0baf00c99b175162b7a,
            ),
            (
                0x78dde6e5fd29f057ce73018173b720d4,
                0x1394f726ab2d8bba1eb7e30c9138,
                0x7399cb3f9cbc592e118ce9a7103,
            ),
            (
                0x1715609f7c746c6dc20fc1e1d0a4e909,
                0x1e729d40f6837953ac8669cbea97,
                0x15fb440086c7e7a02ab930c064b1,
            ),
            (
                0xb54cda58fbbee883b5ac82422d92b13e,
                0x10be68747863a8f1dd20961c21d7,
                0x1df344819b7a8bc381abf8d72f64,
            ),
            (
                0x538454127b096499a94942a28a807973,
                0xb598d0399cccd663d2a5758d5ef,
                0x41a1dc904e35b326c6b57cc14be,
            ),
            (
                0xf1bbcdcbfa53e0af9ce60302e76e41a8,
                0x1abe26662428dabe8c51bea97f15,
                0x956d9910a741fe3a4e43c261529,
            ),
        ];
        let c = Curve::new(crate::field::gf2_109::from_u128(
            0x5f3c_8a71_2e94_d0b6_1c47_a389_e5f2_6d1b,
        ));
        crate::curve::binary::check_map_vectors(&c, &VECTORS);
    }
}
