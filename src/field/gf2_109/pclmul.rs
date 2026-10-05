//! x86_64 backend: PCLMULQDQ, elements in one SSE register.

use super::tables::RHO;
use core::arch::x86_64::*;
use core::mem::transmute;

#[derive(Clone, Copy)]
pub struct Gf(__m128i);

impl Gf {
    /// Low limb first; any 128 bits (reduced lazily).
    pub const fn w64le(x0: u64, x1: u64) -> Self {
        unsafe { Self(transmute::<[u64; 2], __m128i>([x0, x1])) }
    }

    /// The 128-bit representative, not reduced.
    #[inline(always)]
    pub fn limbs(self) -> [u64; 2] {
        unsafe { transmute(self.0) }
    }
}

#[inline(always)]
pub fn add(a: Gf, b: Gf) -> Gf {
    unsafe { Gf(_mm_xor_si128(a.0, b.0)) }
}

/// d0 + d1 z^128 mod z^19 f: fold d1 onto rho twice (see the parent module).
#[inline(always)]
unsafe fn reduce(d0: __m128i, d1: __m128i) -> Gf {
    unsafe {
        let rho = _mm_set1_epi64x(RHO as i64);
        // c = d1_hi rho; its high word is floor(d1 rho / z^128)
        let c = _mm_clmulepi64_si128(d1, rho, 0x11);
        let g = _mm_xor_si128(d1, _mm_bsrli_si128(c, 8));
        let e = _mm_clmulepi64_si128(g, rho, 0x00);
        Gf(_mm_xor_si128(_mm_xor_si128(d0, e), _mm_bslli_si128(c, 8)))
    }
}

/// 4 PCLMULQDQ for the product, 2 for the reduction.
#[inline(always)]
pub fn mul(a: Gf, b: Gf) -> Gf {
    unsafe {
        let (a, b) = (a.0, b.0);
        let c0 = _mm_clmulepi64_si128(a, b, 0x00);
        let c1 = _mm_xor_si128(
            _mm_clmulepi64_si128(a, b, 0x01),
            _mm_clmulepi64_si128(a, b, 0x10),
        );
        let c2 = _mm_clmulepi64_si128(a, b, 0x11);
        reduce(
            _mm_xor_si128(c0, _mm_bslli_si128(c1, 8)),
            _mm_xor_si128(c2, _mm_bsrli_si128(c1, 8)),
        )
    }
}

#[inline(always)]
pub fn square(a: Gf) -> Gf {
    unsafe {
        reduce(
            _mm_clmulepi64_si128(a.0, a.0, 0x00),
            _mm_clmulepi64_si128(a.0, a.0, 0x11),
        )
    }
}
