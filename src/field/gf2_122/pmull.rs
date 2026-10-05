//! aarch64 backend: PMULL (the crypto extension). A GF(2^122) element is
//! one NEON register, a0 in lane 0 and a1 in lane 1.

use super::tables::RHO;
use core::arch::aarch64::*;
use core::mem::transmute;

#[derive(Clone, Copy)]
pub struct Gf(uint64x2_t);

impl Gf {
    /// a0 + a1 u, low word first; any 64-bit representatives.
    pub const fn w64le(x0: u64, x1: u64) -> Self {
        unsafe { Self(transmute::<[u64; 2], uint64x2_t>([x0, x1])) }
    }

    /// The two 64-bit representatives, not reduced.
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

/// Lane 0 of the result is c = lo + hi z^64 mod z^3 f (lane 1 is junk):
/// x = hi rho has T in its high word, then lo + low((hi + T) rho), 2 PMULL.
#[inline(always)]
unsafe fn fold(c: uint64x2_t, rho: uint64x2_t) -> uint64x2_t {
    unsafe {
        let x = pmull_hi(c, rho);
        // lane 1: hi + T
        let e = pmull_hi(veorq_u64(c, x), rho);
        veorq_u64(c, e)
    }
}

#[inline(always)]
pub fn mul61(a: u64, b: u64) -> u64 {
    unsafe {
        let c = vreinterpretq_u64_p128(vmull_p64(a, b));
        vgetq_lane_u64(fold(c, vdupq_n_u64(RHO)), 0)
    }
}

#[inline(always)]
pub fn square61(a: u64) -> u64 {
    mul61(a, a)
}

/// 3 PMULL for the Karatsuba products, 4 for the two folds.
#[inline(always)]
pub fn mul(a: Gf, b: Gf) -> Gf {
    unsafe {
        let (a, b) = (a.0, b.0);
        let rho = vdupq_n_u64(RHO);
        let p00 = pmull_lo(a, b);
        let p11 = pmull_hi(a, b);
        // lane 0: a0 + a1
        let sa = veorq_u64(a, vextq_u64(a, a, 1));
        let sb = veorq_u64(b, vextq_u64(b, b, 1));
        let pm = pmull_lo(sa, sb);
        let c0 = fold(veorq_u64(p00, p11), rho);
        let c1 = fold(veorq_u64(pm, p00), rho);
        Gf(vzip1q_u64(c0, c1))
    }
}

/// (a0 + a1)^2 and a1^2: 2 PMULL, 4 for the folds.
#[inline(always)]
pub fn square(a: Gf) -> Gf {
    unsafe {
        let rho = vdupq_n_u64(RHO);
        // (a0 + a1, a1)
        let v = veorq_u64(a.0, vextq_u64(a.0, vdupq_n_u64(0), 1));
        let c0 = fold(pmull_lo(v, v), rho);
        let c1 = fold(pmull_hi(v, v), rho);
        Gf(vzip1q_u64(c0, c1))
    }
}
