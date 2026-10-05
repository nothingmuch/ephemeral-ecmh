//! `binary::Pornin`'s constants for z^127 + z^63 + 1, a = 1: crrl's own,
//! from GLS254, whose base field this is.
//!
//! - k = z^2: Tr(z^2) = Tr(z) = 0.
//! - Tr(v) = v_0 (Tr(z^i) = 1 for i < 127 only at i = 0), so Tr(c) = 1 is
//!   c_0 = 1.
//! - Tr(c^2/z^2) = Tr(c/z) = c_1: c/z is c >> 1 plus c_0 z^-1 =
//!   c_0 (z^126 + z^62), of trace 0. So Tr(c/z) = 0 is c_1 = 0.
//! - Tr(1) = 1, so the trace tells w from w + 1; on w's canonical form it
//!   is bit 0, crrl's choice for GLS254.

use super::M127;
use crate::curve::binary::Pornin;
use crate::field::gf2_127::{Gf, from_u128};

impl Pornin for M127 {
    #[inline(always)]
    fn k() -> Gf {
        Gf::w64le(4, 0)
    }
    /// c_0 = 1 and c_1 = 0.
    #[inline(always)]
    fn force(v: u128) -> Gf {
        from_u128((v & !2) | 1)
    }
    #[inline(always)]
    fn sq_over_k(c: Gf) -> Gf {
        c.square().div_z2()
    }
    #[inline(always)]
    fn root_bit(w: Gf) -> u32 {
        w.trace()
    }
}

#[cfg(test)]
mod tests {
    crate::curve::binary::map_suite!(
        crate::curve::binary127::M127,
        crate::curve::binary127::tests::curve()
    );

    // Check sage/pornin_map.sage against the Rust map pinned by tests/kat.rs:
    // this grounds the script also used for the 109- and 122-bit fixtures.
    #[test]
    fn matches_sage() {
        // sage/pornin_map.sage 127: (input, x, y) on B = 0x5f3c8a712e94d0b61c47a389e5f26d1b
        const VECTORS: [(u128, u128, u128); 8] = [
            (
                0x9e3779b97f4a7c15f39cc0605cedc835,
                0x142f5a1e102997516f358fc8dd15d379,
                0x3fb6222d1b865470278039c22ef1389b,
            ),
            (
                0x3c6ef372fe94f82be73980c0b9db906a,
                0x190b722bc8ef266961267838e7e25f27,
                0x4f54a0000a0ff4246e73d92dab3c8701,
            ),
            (
                0xdaa66d2c7ddf7441dad6412116c9589f,
                0x4d33b1898c3a593d9c5d20b1f25f8957,
                0x1b677091e800ce9a5fed274ea736fcc7,
            ),
            (
                0x78dde6e5fd29f057ce73018173b720d4,
                0x116843bc5b3652ef5342d54330a4ea73,
                0x53d4738db8087aa99835732dba32a188,
            ),
            (
                0x1715609f7c746c6dc20fc1e1d0a4e909,
                0x29701cb63a09e436a9318c0b0b19fd65,
                0x4783b01f43184a9926f0b5041eef1c29,
            ),
            (
                0xb54cda58fbbee883b5ac82422d92b13e,
                0x457a15091c5a410ef25cc9d0c27c25f3,
                0x7e614aa21940a1a7c2875d0d0c54a1,
            ),
            (
                0x538454127b096499a94942a28a807973,
                0x3a264d7e1c910068837d2a6acc88ef6b,
                0x1f60a68a68d19194a79fb264af2cdff7,
            ),
            (
                0xf1bbcdcbfa53e0af9ce60302e76e41a8,
                0x1f132f76ba8a7f32ca74e4f884a07169,
                0x24c9191ac06847967ccd6a78d1abe621,
            ),
        ];
        let c = Curve::new(crate::field::gf2_127::from_u128(
            0x5f3c_8a71_2e94_d0b6_1c47_a389_e5f2_6d1b,
        ));
        crate::curve::binary::check_map_vectors(&c, &VECTORS);
    }
}
