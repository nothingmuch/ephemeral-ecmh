//! F_p, p = 2^127 - 1: the field of `edwards127` and `weier127`.
//!
//! - Modulus: the Mersenne prime 2^127 - 1.
//! - Representation: `Fp`, one canonical `u128` in [0, p).
//! - Backend: portable Rust, the same on every target. A product is four
//!   64x64 -> 128 multiplications; the Mersenne modulus makes reduction two
//!   shifts and an add (and one conditional subtraction), and makes sqrt a
//!   plain 125-fold squaring since (p + 1)/4 = 2^125.

use core::ops::{Add, Mul, Neg, Sub};

pub const P: u128 = (1 << 127) - 1;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Fp(u128);

impl Fp {
    pub const ZERO: Self = Self(0);
    pub const ONE: Self = Self(1);

    /// Reduces any u128 to the canonical representative in [0, p).
    #[inline]
    pub const fn new(v: u128) -> Self {
        Self::reduce(v)
    }

    #[inline]
    pub const fn value(self) -> u128 {
        self.0
    }

    #[inline]
    const fn reduce(v: u128) -> Self {
        let s = (v & P) + (v >> 127);
        Self(if s >= P { s - P } else { s })
    }

    #[inline]
    pub fn is_zero(self) -> bool {
        self.0 == 0
    }

    #[inline]
    pub fn square(self) -> Self {
        self * self
    }

    pub fn xsquare(self, n: u32) -> Self {
        (0..n).fold(self, |x, _| x.square())
    }

    /// x^(2^125 - 1), shared by inversion and square-root ratios.
    fn pow_2k_minus_1(x1: Self) -> Self {
        let x2 = x1.xsquare(1) * x1;
        let x4 = x2.xsquare(2) * x2;
        let x8 = x4.xsquare(4) * x4;
        let x16 = x8.xsquare(8) * x8;
        let x32 = x16.xsquare(16) * x16;
        let x64 = x32.xsquare(32) * x32;
        let x96 = x64.xsquare(32) * x32;
        let x112 = x96.xsquare(16) * x16;
        let x120 = x112.xsquare(8) * x8;
        let x124 = x120.xsquare(4) * x4;
        x124.xsquare(1) * x1
    }

    /// x^(p-2) = x^(4 * (2^125 - 1) + 1); maps 0 to 0.
    pub fn invert(self) -> Self {
        Self::pow_2k_minus_1(self).xsquare(2) * self
    }

    /// x^((p-3)/4) = x^(2^125 - 1): with it, s = n (n d)^((p-3)/4) satisfies
    /// s^2 d = n chi(n d), a square root of n/d or a non-residue witness
    /// with no separate inversion.
    pub fn pow_p34(self) -> Self {
        Self::pow_2k_minus_1(self)
    }

    /// A square root of n/d, if d != 0 and n/d is a square: n (n d)^((p-3)/4),
    /// one exponentiation and no inversion.
    pub fn sqrt_ratio(n: Self, d: Self) -> Option<Self> {
        if d.is_zero() {
            return None;
        }
        let s = n * (n * d).pow_p34();
        (s.square() * d == n).then_some(s)
    }

    /// Quadratic character via the binary Jacobi symbol (variable time):
    /// roughly 2 * 127 shift-and-subtract steps instead of a 126-bit
    /// exponentiation. 0 counts as a square.
    pub fn is_square(self) -> bool {
        let (mut a, mut n, mut t) = (self.0, P, false);
        while a != 0 {
            let z = a.trailing_zeros();
            a >>= z;
            // (2/n) = -1 iff n = 3, 5 mod 8
            if z & 1 == 1 && matches!(n & 7, 3 | 5) {
                t = !t;
            }
            if a < n {
                // reciprocity: flip iff a = n = 3 mod 4
                if a & n & 3 == 3 {
                    t = !t;
                }
                core::mem::swap(&mut a, &mut n);
            }
            a -= n;
        }
        !t
    }

    /// Some square root if `self` is a square.
    pub fn sqrt(self) -> Option<Self> {
        let y = self.xsquare(125);
        (y.square() == self).then_some(y)
    }

    #[inline]
    pub fn is_odd(self) -> bool {
        self.0 & 1 == 1
    }
}

impl Add for Fp {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self::reduce(self.0 + rhs.0)
    }
}

impl Sub for Fp {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        // p - b, as b ^ p for canonical b: in the builds inspected, LLVM
        // rewrote a + (p - b) into (a - b) + p but leaves the xor alone
        Self::reduce(self.0 + (rhs.0 ^ P))
    }
}

impl Neg for Fp {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        // For canonical a, a ^ P = p - a.
        Self::reduce(self.0 ^ P)
    }
}

impl Mul for Fp {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        let (a0, a1) = (self.0 as u64 as u128, self.0 >> 64);
        let (b0, b1) = (rhs.0 as u64 as u128, rhs.0 >> 64);
        // a1, b1 < 2^63, so the middle sum cannot overflow.
        let mid = a0 * b1 + a1 * b0;
        let (lo, c) = (a0 * b0).overflowing_add(mid << 64);
        let hi = a1 * b1 + (mid >> 64) + c as u128;
        // a*b < 2^254: t = (a*b) >> 127 < 2^127, so t + low < 2^128.
        let t = (hi << 1) | (lo >> 127);
        Self::reduce(t + (lo & P))
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use proptest::prelude::*;

    /// Double-and-add reference, independent of the limb multiply.
    pub fn mul_ref(a: u128, b: u128) -> u128 {
        let mut acc = 0u128;
        for i in (0..127).rev() {
            acc = (acc << 1) % P;
            if b >> i & 1 == 1 {
                acc = (acc + a) % P;
            }
        }
        acc
    }

    pub fn fp() -> impl Strategy<Value = Fp> {
        prop_oneof![
            (0..P).prop_map(Fp),
            // Edges of the limb split and of the modulus.
            prop::sample::select(vec![
                0,
                1,
                2,
                P - 1,
                P - 2,
                1 << 64,
                (1 << 64) - 1,
                1 << 126
            ])
            .prop_map(Fp),
        ]
    }

    proptest! {
        #[test]
        fn canonical_after_ops(a in fp(), b in fp()) {
            for v in [a + b, a - b, a * b, -a] {
                prop_assert!(v.0 < P);
            }
        }

        #[test]
        fn new_reduces(v in any::<u128>()) {
            prop_assert_eq!(Fp::new(v).0, v % P);
        }

        #[test]
        fn mul_matches_reference(a in fp(), b in fp()) {
            prop_assert_eq!((a * b).0, mul_ref(a.0, b.0));
        }

        #[test]
        fn add_sub_neg(a in fp(), b in fp()) {
            prop_assert_eq!((a + b).0, (a.0 + b.0) % P);
            prop_assert_eq!(a - b + b, a);
            prop_assert_eq!(a + -a, Fp::ZERO);
        }

        #[test]
        fn distributive(a in fp(), b in fp(), c in fp()) {
            prop_assert_eq!(a * (b + c), a * b + a * c);
        }

        #[test]
        fn invert(a in fp()) {
            let i = a.invert();
            prop_assert_eq!(a * i, if a.is_zero() { Fp::ZERO } else { Fp::ONE });
        }

        #[test]
        fn is_square_matches_euler(a in fp()) {
            prop_assert_eq!(a.is_square(), a.sqrt().is_some());
            prop_assert!((a * a).is_square());
        }

        #[test]
        fn pow_p34_gives_sqrt_ratio(n in fp(), d in fp()) {
            prop_assume!(!d.is_zero());
            let s = n * (n * d).pow_p34();
            let expect = if (n * d.invert()).is_square() { n } else { -n };
            prop_assert_eq!(s.square() * d, expect);
        }

        #[test]
        fn sqrt_of_square(a in fp()) {
            let s = a.square().sqrt().unwrap();
            prop_assert!(s == a || s == -a);
        }

        #[test]
        fn sqrt_rejects_nonsquares(a in fp()) {
            // p = 3 mod 4: -1 is a non-square, so exactly one of ±a² is a square.
            prop_assume!(!a.is_zero());
            prop_assert!((-a.square()).sqrt().is_none());
        }
    }
}
