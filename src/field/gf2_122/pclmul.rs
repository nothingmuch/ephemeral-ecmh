//! x86_64 backend: PCLMULQDQ. A GF(2^122) element is one SSE register, a0
//! in the low quadword and a1 in the high one. Everything else is SSE2.

use super::tables::RHO;
use core::arch::x86_64::*;
use core::mem::transmute;

#[derive(Clone, Copy)]
pub struct Gf(__m128i);

impl Gf {
    /// a0 + a1 u, low word first; any 64-bit representatives.
    pub const fn w64le(x0: u64, x1: u64) -> Self {
        unsafe { Self(transmute::<[u64; 2], __m128i>([x0, x1])) }
    }

    /// The two 64-bit representatives, not reduced.
    #[inline(always)]
    pub fn limbs(self) -> [u64; 2] {
        unsafe { transmute(self.0) }
    }
}

#[inline(always)]
pub fn add(a: Gf, b: Gf) -> Gf {
    unsafe { Gf(_mm_xor_si128(a.0, b.0)) }
}

/// The low quadword of the result is c = lo + hi z^64 mod z^3 f (the high
/// one is junk): x = hi rho has T in its high word, then
/// lo + low((hi + T) rho), 2 PCLMULQDQ.
#[inline(always)]
unsafe fn fold(c: __m128i, rho: __m128i) -> __m128i {
    unsafe {
        let x = _mm_clmulepi64_si128(c, rho, 0x01);
        // high quadword: hi + T
        let e = _mm_clmulepi64_si128(_mm_xor_si128(c, x), rho, 0x01);
        _mm_xor_si128(c, e)
    }
}

#[inline(always)]
pub fn mul61(a: u64, b: u64) -> u64 {
    unsafe {
        let c = _mm_clmulepi64_si128(_mm_cvtsi64_si128(a as i64), _mm_cvtsi64_si128(b as i64), 0);
        _mm_cvtsi128_si64(fold(c, _mm_set1_epi64x(RHO as i64))) as u64
    }
}

#[inline(always)]
pub fn square61(a: u64) -> u64 {
    mul61(a, a)
}

/// 3 PCLMULQDQ for the Karatsuba products, 4 for the two folds.
#[inline(always)]
pub fn mul(a: Gf, b: Gf) -> Gf {
    unsafe {
        let (a, b) = (a.0, b.0);
        let rho = _mm_set1_epi64x(RHO as i64);
        let p00 = _mm_clmulepi64_si128(a, b, 0x00);
        let p11 = _mm_clmulepi64_si128(a, b, 0x11);
        // low quadword: a0 + a1
        let sa = _mm_xor_si128(a, _mm_shuffle_epi32(a, 0x4E));
        let sb = _mm_xor_si128(b, _mm_shuffle_epi32(b, 0x4E));
        let pm = _mm_clmulepi64_si128(sa, sb, 0x00);
        let c0 = fold(_mm_xor_si128(p00, p11), rho);
        let c1 = fold(_mm_xor_si128(pm, p00), rho);
        Gf(_mm_unpacklo_epi64(c0, c1))
    }
}

/// (a0 + a1)^2 and a1^2: 2 PCLMULQDQ, 4 for the folds.
#[inline(always)]
pub fn square(a: Gf) -> Gf {
    unsafe {
        let rho = _mm_set1_epi64x(RHO as i64);
        // (a0 + a1, a1)
        let v = _mm_xor_si128(a.0, _mm_srli_si128(a.0, 8));
        let c0 = fold(_mm_clmulepi64_si128(v, v, 0x00), rho);
        let c1 = fold(_mm_clmulepi64_si128(v, v, 0x11), rho);
        Gf(_mm_unpacklo_epi64(c0, c1))
    }
}
