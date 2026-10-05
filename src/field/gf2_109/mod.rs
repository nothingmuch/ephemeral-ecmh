//! GF(2^109) = F_2\[z\]/(f), f = z^109 + z^5 + z^4 + z^2 + 1: the field
//! of `binary109`, the 14-byte binary field.
//!
//! - Modulus: the pentanomial f, so a degree-109 extension of F_2
//!   directly, with no subfield tower.
//! - Representation: `Gf`, 128 bits reduced lazily (below); `to_u128`
//!   gives the canonical 109-bit polynomial.
//! - Backend, chosen at compile time by target feature, all three with the
//!   same `Gf` API: `pmull` (aarch64 with `aes`: PMULL on a NEON register),
//!   `pclmul` (x86_64 with `pclmulqdq`: PCLMULQDQ on an SSE register), else
//!   `soft` (carry-less products from integer multiplications). Tests build
//!   `soft` on every target too, and check it against a bit-serial reference.
//!
//! There is no irreducible trinomial of degree 109, and f is the lowest
//! pentanomial (sage/gf2_109.sage). Elements are kept in 128 bits and
//! reduced only modulo z^19 f = z^128 + rho, rho = z^24 + z^23 + z^21 + z^19,
//! as crrl does for GF(2^127) with z f. A 256-bit product d0 + d1 z^128 then
//! folds with two word-aligned carry-less multiplications by rho: with
//! c = d1_hi rho, whose high word T = floor(d1 rho / z^128) has 24 bits,
//! d1 rho = (d1 + T) rho mod z^128 (mod z^19 f), as deg(T rho) < 128.
//! Any f with its middle terms below z^45 costs the same.
//!
//! `to_u128` normalizes; `eq`, `is_zero` and encodings go through it.

use core::fmt;
use core::ops::{Add, AddAssign, Div, Mul, MulAssign};

mod tables;

#[cfg(all(target_arch = "aarch64", target_feature = "aes"))]
mod pmull;
#[cfg(all(target_arch = "aarch64", target_feature = "aes"))]
use pmull as backend;

#[cfg(all(target_arch = "x86_64", target_feature = "pclmulqdq"))]
mod pclmul;
#[cfg(all(target_arch = "x86_64", target_feature = "pclmulqdq"))]
use pclmul as backend;

#[cfg(any(
    test,
    not(any(
        all(target_arch = "aarch64", target_feature = "aes"),
        all(target_arch = "x86_64", target_feature = "pclmulqdq"),
    ))
))]
mod soft;
#[cfg(not(any(
    all(target_arch = "aarch64", target_feature = "aes"),
    all(target_arch = "x86_64", target_feature = "pclmulqdq"),
)))]
use soft as backend;

pub use backend::Gf;

pub const MASK109: u128 = (1 << 109) - 1;

/// Canonical 109-bit integer form (bit i = coefficient of z^i).
#[inline]
pub fn to_u128(x: Gf) -> u128 {
    let [lo, hi] = x.limbs();
    // z^109 t = (z^5 + z^4 + z^2 + 1) t, and t (19 bits) times that fits in lo
    let t = hi >> 45;
    let lo = lo ^ t ^ t << 2 ^ t << 4 ^ t << 5;
    lo as u128 | ((hi & ((1 << 45) - 1)) as u128) << 64
}

/// Bits above 108 are dropped (not reduced), so the map is onto and exact.
#[inline]
pub fn from_u128(v: u128) -> Gf {
    let v = v & MASK109;
    Gf::w64le(v as u64, (v >> 64) as u64)
}

#[inline]
pub fn is_zero(x: Gf) -> bool {
    to_u128(x) == 0
}

#[inline]
pub fn eq(a: Gf, b: Gf) -> bool {
    is_zero(a + b)
}

/// The even-indexed bits of x, packed into the low half.
fn squeeze(x: u64) -> u64 {
    let mut x = x & 0x5555_5555_5555_5555;
    x = (x | x >> 1) & 0x3333_3333_3333_3333;
    x = (x | x >> 2) & 0x0F0F_0F0F_0F0F_0F0F;
    x = (x | x >> 4) & 0x00FF_00FF_00FF_00FF;
    x = (x | x >> 8) & 0x0000_FFFF_0000_FFFF;
    (x | x >> 16) & 0x0000_0000_FFFF_FFFF
}

impl Gf {
    pub const ZERO: Self = Self::w64le(0, 0);
    pub const ONE: Self = Self::w64le(1, 0);

    #[inline(always)]
    pub fn square(self) -> Self {
        backend::square(self)
    }

    pub fn xsquare(self, n: u32) -> Self {
        (0..n).fold(self, |x, _| x.square())
    }

    /// Itoh-Tsujii: 1/a = (a^2)^(2^108 - 1), 8M + 108S. Maps 0 to 0.
    pub fn invert(self) -> Self {
        let a1 = self.square();
        let a2 = a1 * a1.square();
        let a3 = a1 * a2.square();
        let a6 = a3 * a3.xsquare(3);
        let a12 = a6 * a6.xsquare(6);
        let a24 = a12 * a12.xsquare(12);
        let a27 = a3 * a24.xsquare(3);
        let a54 = a27 * a27.xsquare(27);
        a54 * a54.xsquare(54)
    }

    /// a = ae + z ao with ae, ao squares of their squeezed halves, so
    /// sqrt(a) = sqrt(ae) + sqrt(z) sqrt(ao): 1M (sqrt(z) is dense).
    pub fn sqrt(self) -> Self {
        let [lo, hi] = self.limbs();
        let e = squeeze(lo) | squeeze(hi) << 32;
        let o = squeeze(lo >> 1) | squeeze(hi >> 1) << 32;
        let sz = tables::SQRT_Z;
        Self::w64le(e, 0) + Self::w64le(o, 0) * Self::w64le(sz as u64, (sz >> 64) as u64)
    }

    /// Tr(z^i) = 1 for i < 128 exactly at i = 0, 105, 107, 109, so this
    /// needs no normalization.
    pub fn trace(self) -> u32 {
        let [lo, hi] = self.limbs();
        ((lo ^ hi >> 41 ^ hi >> 43 ^ hi >> 45) & 1) as u32
    }

    /// H(a) = sum of a^(4^i), i = 0..=54: H(v)^2 + H(v) = v + Tr(v).
    ///
    /// Odd-indexed bits go through a table. Even ones are a square,
    /// H(s^2) = H(s) + s + Tr(s), and s is half as long: repeat until s has
    /// one bit, then H(1) = 1 (55 terms). With e the sum of the s, the
    /// Tr(s) sum to Tr(e) = e_0 (every s is shorter than z^105).
    pub fn halftrace(self) -> Self {
        let v = to_u128(self);
        let (lo, hi) = (v as u64, (v >> 64) as u64);
        let (mut odd, odd_hi) = (lo, hi);
        let mut s = squeeze(lo) | squeeze(hi) << 32;
        let mut e = s;
        for _ in 0..6 {
            odd ^= s;
            s = squeeze(s);
            e ^= s;
        }
        // s is now the top bit, H(s) = s
        let mut d = [e ^ (e & 1) ^ s, 0];
        for (i, &h) in tables::HALFTRACE.iter().enumerate() {
            let w = if i < 32 { odd } else { odd_hi };
            let m = (w >> (2 * (i % 32) + 1) & 1).wrapping_neg();
            d[0] ^= m & h as u64;
            d[1] ^= m & (h >> 64) as u64;
        }
        Self::w64le(d[0], d[1])
    }
}

/// H(z^i) for every i < 109: odd i from the table, even i = 2k by
/// H(z^(2k)) = H(z^k) + z^k + Tr(z^k), where Tr(z^k) = 0 as 0 < k < 105.
const fn monomial_halftraces() -> [u128; 109] {
    let mut h = [0; 109];
    h[0] = 1;
    let mut i = 1;
    while i < 109 {
        h[i] = if i % 2 == 1 {
            tables::HALFTRACE[i / 2]
        } else {
            h[i / 2] ^ 1 << (i / 2)
        };
        i += 1;
    }
    h
}

static HT8: [[u128; 256]; 14] = super::window_tables(&monomial_halftraces());

/// The half-trace by byte tables: 14 lookups into 56 KiB, where
/// `Gf::halftrace` folds the even bits and walks 54 odd-bit rows.
#[inline(always)]
pub fn halftrace8(x: Gf) -> Gf {
    let v = to_u128(x);
    let mut h = 0;
    for (j, t) in HT8.iter().enumerate() {
        h ^= t[(v >> (8 * j)) as u8 as usize];
    }
    from_u128(h)
}

impl fmt::Debug for Gf {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Gf({:#030x})", to_u128(*self))
    }
}

impl Add for Gf {
    type Output = Self;
    #[inline(always)]
    fn add(self, o: Self) -> Self {
        backend::add(self, o)
    }
}

impl AddAssign for Gf {
    #[inline(always)]
    fn add_assign(&mut self, o: Self) {
        *self = *self + o;
    }
}

impl Mul for Gf {
    type Output = Self;
    #[inline(always)]
    fn mul(self, o: Self) -> Self {
        backend::mul(self, o)
    }
}

impl MulAssign for Gf {
    #[inline(always)]
    fn mul_assign(&mut self, o: Self) {
        *self = *self * o;
    }
}

impl Div for Gf {
    type Output = Self;
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, o: Self) -> Self {
        self * o.invert()
    }
}

impl crate::field::batch::Invert for Gf {
    const ONE: Self = Gf::ONE;
    fn inv(self) -> Self {
        self.invert()
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use proptest::prelude::*;

    /// Bit-serial reference multiply, mod f, independent of the backends.
    pub fn mul_ref(a: u128, b: u128) -> u128 {
        let (mut acc, mut a) = (0u128, a & MASK109);
        for i in 0..109 {
            if b >> i & 1 == 1 {
                acc ^= a;
            }
            // a *= z, then z^109 = z^5 + z^4 + z^2 + 1
            a <<= 1;
            if a >> 109 & 1 == 1 {
                a ^= 1 << 109 | 0b110101;
            }
        }
        acc
    }

    fn pow_ref(a: u128, e: u32) -> u128 {
        // a^(2^e)
        (0..e).fold(a, |x, _| mul_ref(x, x))
    }

    pub fn fe() -> impl Strategy<Value = u128> {
        any::<u128>().prop_map(|v| v & MASK109)
    }

    /// Any 128-bit representative, as products leave them.
    fn raw() -> impl Strategy<Value = Gf> {
        any::<u128>().prop_map(|v| Gf::w64le(v as u64, (v >> 64) as u64))
    }

    /// v mod f for any 128-bit v, bit by bit.
    fn reduce_ref(mut v: u128) -> u128 {
        for i in (109..128).rev() {
            if v >> i & 1 == 1 {
                v ^= 1 << i | 0b110101 << (i - 109);
            }
        }
        v
    }

    proptest! {
        #[test]
        fn roundtrip(a in fe()) {
            prop_assert_eq!(to_u128(from_u128(a)), a);
        }

        #[test]
        fn normalization_reduces_mod_f(v in any::<u128>()) {
            prop_assert_eq!(to_u128(Gf::w64le(v as u64, (v >> 64) as u64)), reduce_ref(v));
        }

        #[test]
        fn mul_matches_reference(a in raw(), b in raw()) {
            prop_assert_eq!(to_u128(a * b), mul_ref(to_u128(a), to_u128(b)));
        }

        #[test]
        fn square_matches_mul(a in raw()) {
            let v = to_u128(a);
            prop_assert_eq!(to_u128(a.square()), mul_ref(v, v));
        }

        #[test]
        fn soft_backend_matches_reference(a in any::<u128>(), b in any::<u128>()) {
            let s = |v: u128| soft::Gf::w64le(v as u64, (v >> 64) as u64);
            let n = |x: soft::Gf| to_u128(Gf::w64le(x.limbs()[0], x.limbs()[1]));
            let (ra, rb) = (reduce_ref(a), reduce_ref(b));
            prop_assert_eq!(n(soft::mul(s(a), s(b))), mul_ref(ra, rb));
            prop_assert_eq!(n(soft::square(s(a))), mul_ref(ra, ra));
            prop_assert_eq!(n(soft::add(s(a), s(b))), ra ^ rb);
        }

        #[test]
        fn soft_clmul_matches_shift_and_add(a in any::<u64>(), b in any::<u64>()) {
            let want = (0..64).filter(|i| b >> i & 1 == 1).fold(0u128, |acc, i| acc ^ (a as u128) << i);
            prop_assert_eq!(soft::clmul(a, b), want);
        }

        #[test]
        fn invert(a in raw()) {
            let y = a.invert();
            if is_zero(a) {
                prop_assert!(is_zero(y));
            } else {
                prop_assert_eq!(mul_ref(to_u128(a), to_u128(y)), 1);
            }
        }

        #[test]
        fn sqrt(a in raw()) {
            prop_assert!(eq(a.sqrt().square(), a));
        }

        #[test]
        fn trace_matches_definition(a in raw()) {
            // Tr(v) = sum of v^(2^i), i < 109
            let v = to_u128(a);
            let t = (0..109).fold((0, v), |(t, x), _| (t ^ x, mul_ref(x, x))).0;
            prop_assert!(t <= 1);
            prop_assert_eq!(a.trace() as u128, t);
        }

        #[test]
        fn halftrace_matches_definition(a in raw()) {
            let v = to_u128(a);
            let want = (0..55).fold((0, v), |(h, x), _| (h ^ x, pow_ref(x, 2))).0;
            prop_assert_eq!(to_u128(a.halftrace()), want);
            let h = a.halftrace();
            prop_assert!(eq(h.square() + h + a, Gf::w64le(a.trace() as u64, 0)));
        }

        #[test]
        fn halftrace_tables_match_rows(a in raw()) {
            prop_assert_eq!(to_u128(halftrace8(a)), to_u128(a.halftrace()));
        }

        #[test]
        fn distributive(a in raw(), b in raw(), c in raw()) {
            prop_assert!(eq(a * (b + c), a * b + a * c));
        }

        #[test]
        fn batch_invert_matches_single(v in prop::collection::vec(fe().prop_filter("nonzero", |&x| x != 0), 0..40)) {
            let v: Vec<Gf> = v.into_iter().map(from_u128).collect();
            let mut w = v.clone();
            crate::field::batch::invert(&mut w);
            for (x, y) in v.iter().zip(&w) {
                prop_assert!(eq(x.invert(), *y));
            }
        }
    }

    #[test]
    fn tables_match_reference() {
        use super::tables::*;
        // z^128 = z^19 (z^5 + z^4 + z^2 + 1) mod z^19 f
        assert_eq!(RHO, 0b110101 << 19);
        assert_eq!(mul_ref(SQRT_Z, SQRT_Z), 2);
        for (i, &h) in HALFTRACE.iter().enumerate() {
            let v = 1 << (2 * i + 1);
            let want = (0..55).fold((0, v), |(h, x), _| (h ^ x, pow_ref(x, 2))).0;
            assert_eq!(h, want, "H(z^{})", 2 * i + 1);
        }
    }
}
