//! GF(2^61) = F_2\[z\]/(f), f = z^61 + z^23 + z^15 + z^5 + 1: the base
//! field of the tower (see the parent module for why this f).
//!
//! An element is one u64, reduced only modulo z^3 f = z^64 + rho,
//! rho = z^26 + z^18 + z^8 + z^3: any 64 bits represent their residue, so
//! a 128-bit product lo + hi z^64 folds onto rho at word boundaries. With
//! T = floor(hi rho / z^64) (26 bits), hi rho = (hi + T) rho mod z^64
//! (mod z^3 f), as deg(T rho) < 64: two folds, each 4 shifts or one
//! carry-less multiplication (the backends' `mul61`).
//!
//! `to_u64` normalizes; `eq`, `is_zero` and the tower's encoding go through
//! it. Trace, square root and half-trace take any representative.

use super::backend;
use super::tables::HALFTRACE;
use core::fmt;
use core::ops::{Add, AddAssign, Div, Mul, MulAssign};

pub const MASK61: u64 = (1 << 61) - 1;

/// An element of GF(2^61), any 64-bit representative mod z^3 f.
#[derive(Clone, Copy)]
pub struct Gf61(u64);

/// Canonical 61-bit form (bit i = coefficient of z^i).
#[inline(always)]
pub fn to_u64(x: Gf61) -> u64 {
    // z^61 t = (z^23 + z^15 + z^5 + 1) t, and t < z^3 times that fits
    let (v, t) = (x.0, x.0 >> 61);
    v & MASK61 ^ t ^ t << 5 ^ t << 15 ^ t << 23
}

/// Bits above 60 are dropped (not reduced), so the map is onto and exact.
pub const fn from_u64(v: u64) -> Gf61 {
    Gf61(v & MASK61)
}

#[inline(always)]
pub fn is_zero(x: Gf61) -> bool {
    to_u64(x) == 0
}

#[inline(always)]
pub fn eq(a: Gf61, b: Gf61) -> bool {
    is_zero(a + b)
}

/// The even-indexed bits of x in the low half, the odd ones in the high half.
#[inline(always)]
fn unshuffle(mut x: u64) -> u64 {
    // Hacker's Delight 7-2: swap bit groups inward, 2 bits at a time up to 32
    for (m, s) in [
        (0x2222_2222_2222_2222, 1),
        (0x0C0C_0C0C_0C0C_0C0C, 2),
        (0x00F0_00F0_00F0_00F0, 4),
        (0x0000_FF00_0000_FF00, 8),
        (0x0000_0000_FFFF_0000, 16),
    ] {
        let t = (x ^ x >> s) & m;
        x ^= t ^ t << s;
    }
    x
}

/// Half-trace byte tables: `HT[j][b]` = H(sum of z^(8j + k) over the bits k of b).
static HT: [[u64; 256]; 8] = byte_tables(&HALFTRACE);

const fn byte_tables(h: &[u64; 64]) -> [[u64; 256]; 8] {
    let mut t = [[0; 256]; 8];
    let mut j = 0;
    while j < 8 {
        let mut b = 1;
        while b < 256 {
            // b without its lowest bit is already done
            t[j][b] = t[j][b & (b - 1)] ^ h[8 * j + b.trailing_zeros() as usize];
            b += 1;
        }
        j += 1;
    }
    t
}

impl Gf61 {
    pub const ZERO: Self = Self(0);
    pub const ONE: Self = Self(1);

    /// Any 64 bits, reduced lazily (mod z^3 f).
    #[inline(always)]
    pub const fn w64(v: u64) -> Self {
        Self(v)
    }

    /// The 64-bit representative, not reduced.
    #[inline(always)]
    pub const fn limb(self) -> u64 {
        self.0
    }

    #[inline(always)]
    pub fn square(self) -> Self {
        Self(backend::square61(self.0))
    }

    pub fn xsquare(self, n: u32) -> Self {
        (0..n).fold(self, |x, _| x.square())
    }

    /// Itoh-Tsujii: 1/a = (a^(2^60 - 1))^2, 7M + 60S. Maps 0 to 0.
    pub fn invert(self) -> Self {
        let a2 = self * self.square();
        let a3 = self * a2.square();
        let a6 = a3 * a3.xsquare(3);
        let a12 = a6 * a6.xsquare(6);
        let a15 = a3 * a12.xsquare(3);
        let a30 = a15 * a15.xsquare(15);
        let a60 = a30 * a30.xsquare(30);
        a60.square()
    }

    /// a = e^2 + z o^2 with e, o the even and odd bits squeezed, so
    /// sqrt(a) = e + sqrt(z) o, and sqrt(z) = z^31 + z^12 + z^8 + z^3
    /// (z = z^62 + z^24 + z^16 + z^6 mod f): 4 shifts, no multiplication.
    /// Both halves have 32 bits, so the result fits in 63.
    #[inline(always)]
    pub fn sqrt(self) -> Self {
        let v = unshuffle(self.0);
        let (e, o) = (v & 0xFFFF_FFFF, v >> 32);
        Self(e ^ o << 31 ^ o << 12 ^ o << 8 ^ o << 3)
    }

    /// Tr(z^i) = 1 for i < 64 exactly at i = 0 and 61, so this needs no
    /// normalization; on a canonical element it is bit 0.
    #[inline(always)]
    pub fn trace(self) -> u32 {
        ((self.0 ^ self.0 >> 61) & 1) as u32
    }

    /// H(a) = sum of a^(4^i), i = 0..=30: H(v)^2 + H(v) = v + Tr(v). Linear,
    /// so one byte-table lookup per byte of the representative (16 KiB of
    /// tables). The result is canonical.
    #[inline(always)]
    pub fn halftrace(self) -> Self {
        let v = self.0;
        let mut h = 0;
        for (j, t) in HT.iter().enumerate() {
            h ^= t[(v >> (8 * j)) as u8 as usize];
        }
        Self(h)
    }
}

impl fmt::Debug for Gf61 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Gf61({:#018x})", to_u64(*self))
    }
}

impl Add for Gf61 {
    type Output = Self;
    #[inline(always)]
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn add(self, o: Self) -> Self {
        Self(self.0 ^ o.0)
    }
}

impl AddAssign for Gf61 {
    #[inline(always)]
    fn add_assign(&mut self, o: Self) {
        *self = *self + o;
    }
}

impl Mul for Gf61 {
    type Output = Self;
    #[inline(always)]
    fn mul(self, o: Self) -> Self {
        Self(backend::mul61(self.0, o.0))
    }
}

impl MulAssign for Gf61 {
    #[inline(always)]
    fn mul_assign(&mut self, o: Self) {
        *self = *self * o;
    }
}

impl Div for Gf61 {
    type Output = Self;
    #[allow(clippy::suspicious_arithmetic_impl)]
    fn div(self, o: Self) -> Self {
        self * o.invert()
    }
}

impl crate::field::batch::Invert for Gf61 {
    const ONE: Self = Gf61::ONE;
    fn inv(self) -> Self {
        self.invert()
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::super::soft;
    use super::*;
    use proptest::prelude::*;

    /// f - z^61
    const G: u64 = 1 << 23 | 1 << 15 | 1 << 5 | 1;

    /// Bit-serial reference multiply, mod f, independent of the backends.
    pub fn mul_ref(a: u64, b: u64) -> u64 {
        let (mut acc, mut a, b) = (0u64, reduce_ref(a), reduce_ref(b));
        for i in 0..61 {
            if b >> i & 1 == 1 {
                acc ^= a;
            }
            // a *= z, then z^61 = z^23 + z^15 + z^5 + 1
            a <<= 1;
            if a >> 61 & 1 == 1 {
                a ^= 1 << 61 | G;
            }
        }
        acc
    }

    /// v mod f for any 64-bit v, bit by bit.
    pub fn reduce_ref(mut v: u64) -> u64 {
        for i in (61..64).rev() {
            if v >> i & 1 == 1 {
                v ^= 1 << i | G << (i - 61);
            }
        }
        v
    }

    /// a^(2^e)
    pub fn frob_ref(a: u64, e: u32) -> u64 {
        (0..e).fold(reduce_ref(a), |x, _| mul_ref(x, x))
    }

    pub fn halftrace_ref(a: u64) -> u64 {
        (0..31)
            .fold((0, reduce_ref(a)), |(h, x), _| (h ^ x, frob_ref(x, 2)))
            .0
    }

    /// Tr(v) = sum of v^(2^i), i < 61
    pub fn trace_ref(a: u64) -> u64 {
        let t = (0..61)
            .fold((0, reduce_ref(a)), |(t, x), _| (t ^ x, mul_ref(x, x)))
            .0;
        assert!(t <= 1);
        t
    }

    pub fn fe() -> impl Strategy<Value = u64> {
        any::<u64>().prop_map(|v| v & MASK61)
    }

    /// Any 64-bit representative, as products leave them.
    pub fn raw() -> impl Strategy<Value = Gf61> {
        any::<u64>().prop_map(Gf61::w64)
    }

    proptest! {
        #[test]
        fn roundtrip(a in fe()) {
            prop_assert_eq!(to_u64(from_u64(a)), a);
        }

        #[test]
        fn normalization_reduces_mod_f(v in any::<u64>()) {
            prop_assert_eq!(to_u64(Gf61::w64(v)), reduce_ref(v));
        }

        #[test]
        fn mul_matches_reference(a in raw(), b in raw()) {
            prop_assert_eq!(to_u64(a * b), mul_ref(a.limb(), b.limb()));
        }

        #[test]
        fn square_matches_mul(a in raw()) {
            prop_assert_eq!(to_u64(a.square()), mul_ref(a.limb(), a.limb()));
        }

        #[test]
        fn soft_backend_matches_reference(a in any::<u64>(), b in any::<u64>()) {
            prop_assert_eq!(reduce_ref(soft::mul61(a, b)), mul_ref(a, b));
            prop_assert_eq!(reduce_ref(soft::square61(a)), mul_ref(a, a));
        }

        #[test]
        fn soft_clmul_matches_shift_and_add(a in any::<u64>(), b in any::<u64>()) {
            let want = (0..64).filter(|i| b >> i & 1 == 1).fold(0u128, |acc, i| acc ^ (a as u128) << i);
            prop_assert_eq!(soft::clmul(a, b), want);
            prop_assert_eq!(soft::clmul(a, a), soft::spread(a));
        }

        #[test]
        fn invert(a in raw()) {
            let y = a.invert();
            if is_zero(a) {
                prop_assert!(is_zero(y));
            } else {
                prop_assert_eq!(mul_ref(a.limb(), y.limb()), 1);
            }
        }

        #[test]
        fn sqrt(a in raw()) {
            let s = a.sqrt();
            prop_assert!(eq(s.square(), a));
            prop_assert_eq!(to_u64(s), frob_ref(a.limb(), 60));
        }

        #[test]
        fn trace_matches_definition(a in raw()) {
            prop_assert_eq!(a.trace() as u64, trace_ref(a.limb()));
            prop_assert_eq!(a.trace() as u64, to_u64(a) & 1);
        }

        #[test]
        fn halftrace_matches_definition(a in raw()) {
            let h = a.halftrace();
            prop_assert_eq!(h.limb(), halftrace_ref(a.limb()));
            prop_assert!(eq(h.square() + h + a, Gf61::w64(a.trace() as u64)));
        }

        #[test]
        fn field_axioms(a in raw(), b in raw(), c in raw()) {
            prop_assert!(eq(a * (b + c), a * b + a * c));
            prop_assert!(eq((a * b) * c, a * (b * c)));
            prop_assert!(eq(a * b, b * a));
            prop_assert!(eq(a * Gf61::ONE, a) && is_zero(a * Gf61::ZERO) && is_zero(a + a));
            prop_assert!(eq((a + b).square(), a.square() + b.square()));
        }

        #[test]
        fn batch_invert_matches_single(v in prop::collection::vec(fe().prop_filter("nonzero", |&x| x != 0), 0..40)) {
            let v: Vec<Gf61> = v.into_iter().map(from_u64).collect();
            let mut w = v.clone();
            crate::field::batch::invert(&mut w);
            for (x, y) in v.iter().zip(&w) {
                prop_assert!(eq(x.invert(), *y));
            }
        }
    }

    #[test]
    fn tables_match_reference() {
        use super::super::tables::*;
        // z^64 = z^3 (z^23 + z^15 + z^5 + 1) mod z^3 f
        assert_eq!(RHO, G << 3);
        for (i, &h) in HALFTRACE.iter().enumerate() {
            assert_eq!(h, halftrace_ref(1 << i), "H(z^{i})");
        }
        assert_eq!(
            to_u64(Gf61::w64(2).sqrt()),
            1 << 31 | 1 << 12 | 1 << 8 | 1 << 3
        );
    }

    #[test]
    fn known_answers() {
        use super::super::kats::{MUL61, UNARY61};
        for &(a, b, c) in MUL61 {
            assert_eq!(to_u64(from_u64(a) * from_u64(b)), c);
        }
        for &(a, inv, sqrt, h, tr) in UNARY61 {
            let x = from_u64(a);
            assert_eq!(to_u64(x.invert()), inv, "1/{a:#x}");
            assert_eq!(to_u64(x.sqrt()), sqrt, "sqrt {a:#x}");
            assert_eq!(to_u64(x.halftrace()), h, "H({a:#x})");
            assert_eq!(x.trace(), tr, "Tr {a:#x}");
        }
    }
}
