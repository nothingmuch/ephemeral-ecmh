//! aarch64 backend: PMULL (the crypto extension), elements in one NEON register.

use super::tables::RHO;
use core::arch::aarch64::*;
use core::mem::transmute;

#[derive(Clone, Copy)]
pub struct Gf(uint64x2_t);

impl Gf {
    /// Low limb first; any 128 bits (reduced lazily).
    pub const fn w64le(x0: u64, x1: u64) -> Self {
        unsafe { Self(transmute::<[u64; 2], uint64x2_t>([x0, x1])) }
    }

    /// The 128-bit representative, not reduced.
    #[inline(always)]
    pub fn limbs(self) -> [u64; 2] {
        unsafe { transmute(self.0) }
    }
}

#[inline(always)]
pub fn add(a: Gf, b: Gf) -> Gf {
    unsafe { Gf(veorq_u64(a.0, b.0)) }
}

#[inline(always)]
unsafe fn pmull_lo(a: uint64x2_t, b: uint64x2_t) -> uint64x2_t {
    unsafe { vreinterpretq_u64_p128(vmull_p64(vgetq_lane_u64(a, 0), vgetq_lane_u64(b, 0))) }
}

#[inline(always)]
unsafe fn pmull_hi(a: uint64x2_t, b: uint64x2_t) -> uint64x2_t {
    unsafe {
        vreinterpretq_u64_p128(vmull_high_p64(
            vreinterpretq_p64_u64(a),
            vreinterpretq_p64_u64(b),
        ))
    }
}

/// d0 + d1 z^128 mod z^19 f: fold d1 onto rho twice (see the parent module).
#[inline(always)]
unsafe fn reduce(d0: uint64x2_t, d1: uint64x2_t) -> Gf {
    unsafe {
        let rho = vdupq_n_u64(RHO);
        let z = vdupq_n_u64(0);
        // c = d1_hi rho; its high word is floor(d1 rho / z^128)
        let c = pmull_hi(d1, rho);
        let g = veorq_u64(d1, vextq_u64(c, c, 1));
        let e = pmull_lo(g, rho);
        Gf(veorq_u64(veorq_u64(d0, e), vextq_u64(z, c, 1)))
    }
}

/// 4 PMULL for the product, 2 for the reduction.
#[inline(always)]
pub fn mul(a: Gf, b: Gf) -> Gf {
    unsafe {
        let (a, b) = (a.0, b.0);
        let bx = vextq_u64(b, b, 1);
        let c0 = pmull_lo(a, b);
        let c2 = pmull_hi(a, b);
        let c1 = veorq_u64(pmull_lo(a, bx), pmull_hi(a, bx));
        let z = vdupq_n_u64(0);
        reduce(
            veorq_u64(c0, vextq_u64(z, c1, 1)),
            veorq_u64(c2, vextq_u64(c1, z, 1)),
        )
    }
}

#[inline(always)]
pub fn square(a: Gf) -> Gf {
    unsafe { reduce(pmull_lo(a.0, a.0), pmull_hi(a.0, a.0)) }
}
