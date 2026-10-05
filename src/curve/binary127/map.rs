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
}
