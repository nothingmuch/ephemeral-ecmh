//! Arithmetic in F_p for p = 2^128 - 275, with p = 5 (mod 8).
//!
//! Elements are canonical u128 values. Multiplication forms a 256-bit product
//! and folds its upper half using 2^128 = 275 (mod p). Unlike the Mersenne
//! fields, both the product assembly and the reduction require carry handling.
//! The parameter search in `sage/edwards_search.py` proves this modulus prime.
//! Arithmetic is variable-time and intended for public-input experiments.

use core::ops::{Add, Mul, Neg, Sub};

pub const P: u128 = u128::MAX - 274;
const C: u128 = 275;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Fp(u128);

impl Fp {
    pub const ZERO: Self = Self(0);
    pub const ONE: Self = Self(1);
    /// 2^((p-1)/4), one square root of -1 since 2 is a nonsquare.
    pub const SQRT_M1: Self = Self(293588884923084747205589728388446189287);

    #[inline]
    pub const fn new(value: u128) -> Self {
        Self(if value >= P { value - P } else { value })
    }

    #[inline]
    pub const fn value(self) -> u128 {
        self.0
    }

    #[inline]
    pub fn is_zero(self) -> bool {
        self.0 == 0
    }

    #[inline]
    pub fn is_odd(self) -> bool {
        self.0 & 1 == 1
    }

    #[inline]
    pub fn square(self) -> Self {
        self * self
    }

    pub fn xsquare(self, n: u32) -> Self {
        (0..n).fold(self, |x, _| x.square())
    }

    /// Return x^(2^119-1) and x^7, shared prefixes of the three exponents.
    fn power_prefix(self) -> (Self, Self) {
        let x2 = self.square() * self;
        let x3 = x2.square() * self;
        let x6 = x3.xsquare(3) * x3;
        let x7 = x6.square() * self;
        let x14 = x7.xsquare(7) * x7;
        let x28 = x14.xsquare(14) * x14;
        let x56 = x28.xsquare(28) * x28;
        let x112 = x56.xsquare(56) * x56;
        (x112.xsquare(7) * x7, x3)
    }

    /// x^(p-2); zero maps to zero.
    pub fn invert(self) -> Self {
        let (prefix, x7) = self.power_prefix();
        let x11 = x7 * self.xsquare(2);
        let x235 = x7.xsquare(5) * x11;
        // p-2 = 2^9 * (2^119-1) + 235.
        prefix.xsquare(9) * x235
    }

    /// x^((p-5)/8), the exponent used by the square-root ratio formula.
    pub fn pow_p58(self) -> Self {
        let (prefix, x7) = self.power_prefix();
        // (p-5)/8 = 2^6 * (2^119-1) + 29.
        prefix.xsquare(6) * (x7.xsquare(2) * self)
    }

    /// Some square root, or None for a nonsquare.
    pub fn sqrt(self) -> Option<Self> {
        let (prefix, x7) = self.power_prefix();
        let x30 = (x7.square() * self).square();
        // (p+3)/8 = 2^6 * (2^119-1) + 30.
        let root = prefix.xsquare(6) * x30;
        let square = root.square();
        if square == self {
            Some(root)
        } else if square == -self {
            Some(root * Self::SQRT_M1)
        } else {
            None
        }
    }

    /// A square root of n/d using one exponentiation and no inversion.
    /// A zero denominator is rejected, including the indeterminate 0/0.
    pub fn sqrt_ratio(n: Self, d: Self) -> Option<Self> {
        if d.is_zero() {
            return None;
        }
        let d3 = d.square() * d;
        let d7 = d3.square() * d;
        let root = n * d3 * (n * d7).pow_p58();
        let check = root.square() * d;
        if check == n {
            Some(root)
        } else if check == -n {
            Some(root * Self::SQRT_M1)
        } else {
            None
        }
    }

    /// Quadratic character by the variable-time binary Jacobi algorithm.
    /// Zero is a square.
    pub fn is_square(self) -> bool {
        let (mut a, mut n, mut negative) = (self.0, P, false);
        while a != 0 {
            let zeros = a.trailing_zeros();
            a >>= zeros;
            if zeros & 1 == 1 && matches!(n & 7, 3 | 5) {
                negative = !negative;
            }
            if a < n {
                if a & n & 3 == 3 {
                    negative = !negative;
                }
                core::mem::swap(&mut a, &mut n);
            }
            a -= n;
        }
        !negative
    }
}

impl Add for Fp {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        let (sum, carry) = self.0.overflowing_add(rhs.0);
        // With a carry, sum < 2^128-2C, so adding C cannot overflow.
        if carry { Self(sum + C) } else { Self::new(sum) }
    }
}

impl Sub for Fp {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        let (difference, borrow) = self.0.overflowing_sub(rhs.0);
        // A borrowed difference is at least C+1; subtracting C replaces
        // the implicit 2^128 by p.
        Self(if borrow { difference - C } else { difference })
    }
}

impl Neg for Fp {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        if self.is_zero() {
            self
        } else {
            Self(P - self.0)
        }
    }
}

impl Mul for Fp {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        let (a0, a1) = (self.0 as u64 as u128, self.0 >> 64);
        let (b0, b1) = (rhs.0 as u64 as u128, rhs.0 >> 64);
        let (mid, cross_carry) = (a0 * b1).overflowing_add(a1 * b0);
        let (lo, carry) = (a0 * b0).overflowing_add(mid << 64);
        let hi = a1 * b1 + (mid >> 64) + ((cross_carry as u128) << 64) + carry as u128;

        // hi*C has at most 137 bits; form it without a general wide product.
        let h0 = (hi as u64 as u128) * C;
        let h1 = (hi >> 64) * C + (h0 >> 64);
        let (folded, carry) = lo.overflowing_add((h1 << 64) | (h0 as u64 as u128));
        let top = (h1 >> 64) + carry as u128; // at most C
        let (folded, carry) = folded.overflowing_add(top * C);
        // If the second fold carries, its residue is below C^2. The final
        // C addition therefore cannot carry again.
        Self::new(if carry { folded + C } else { folded })
    }
}

impl crate::field::batch::Invert for Fp {
    const ONE: Self = Fp::ONE;
    fn inv(self) -> Self {
        self.invert()
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use proptest::prelude::*;

    fn add_ref(a: u128, b: u128) -> u128 {
        if a >= P - b { a - (P - b) } else { a + b }
    }

    fn mul_ref(a: u128, mut b: u128) -> u128 {
        let (mut acc, mut term) = (0, a);
        while b != 0 {
            if b & 1 == 1 {
                acc = add_ref(acc, term);
            }
            term = add_ref(term, term);
            b >>= 1;
        }
        acc
    }

    fn pow_ref(a: u128, mut exponent: u128) -> u128 {
        let (mut acc, mut term) = (1, a);
        while exponent != 0 {
            if exponent & 1 == 1 {
                acc = mul_ref(acc, term);
            }
            term = mul_ref(term, term);
            exponent >>= 1;
        }
        acc
    }

    /// Reduced u128 inputs mixed with explicit carry boundaries.
    pub fn fp() -> impl Strategy<Value = Fp> {
        prop_oneof![
            any::<u128>().prop_map(Fp::new),
            prop::sample::select(vec![
                0,
                1,
                2,
                274,
                275,
                P - 1,
                P - 2,
                1 << 64,
                (1 << 64) - 1,
                1 << 127
            ])
            .prop_map(Fp::new)
        ]
    }

    #[test]
    fn carry_boundaries_and_square_root_constant() {
        assert_eq!(Fp::new(u128::MAX).value(), 274);
        assert_eq!(Fp::new(P), Fp::ZERO);
        assert_eq!(Fp::new(P - 1) + Fp::new(P - 1), Fp::new(P - 2));
        assert_eq!(Fp::ZERO - Fp::ONE, Fp::new(P - 1));
        assert_eq!(Fp::new(P - 1).square(), Fp::ONE);
        assert_eq!(Fp::SQRT_M1.square(), -Fp::ONE);
        assert_eq!(Fp::new(2).sqrt(), None);
        assert_eq!(Fp::ZERO.invert(), Fp::ZERO);
        assert_eq!(Fp::sqrt_ratio(Fp::ONE, Fp::ZERO), None);
        assert_eq!(Fp::sqrt_ratio(Fp::ZERO, Fp::ZERO), None);
    }

    proptest! {
        #[test]
        fn arithmetic_matches_independent_reference(a in fp(), b in fp()) {
            prop_assert_eq!((a + b).value(), add_ref(a.value(), b.value()));
            prop_assert_eq!((a * b).value(), mul_ref(a.value(), b.value()));
            prop_assert_eq!(a - b + b, a);
            prop_assert_eq!(a + -a, Fp::ZERO);
            for value in [a + b, a - b, a * b, -a] {
                prop_assert!(value.value() < P);
            }
        }

        #[test]
        fn inversion_matches_reference(a in fp()) {
            prop_assert_eq!(a.invert().value(), pow_ref(a.value(), P - 2));
            prop_assert_eq!(a * a.invert(), if a.is_zero() { Fp::ZERO } else { Fp::ONE });
        }

        #[test]
        fn roots_and_character_agree(a in fp()) {
            let root = a.square().sqrt().unwrap();
            prop_assert!(root == a || root == -a);
            let euler = pow_ref(a.value(), (P - 1) / 2);
            prop_assert_eq!(a.is_square(), euler != P - 1);
            prop_assert_eq!(a.sqrt().is_some(), a.is_square());
        }

        #[test]
        fn square_root_ratio_matches_division(n in fp(), d in fp()) {
            prop_assume!(!d.is_zero());
            let result = Fp::sqrt_ratio(n, d);
            prop_assert_eq!(result.is_some(), (n * d.invert()).is_square());
            if let Some(root) = result {
                prop_assert_eq!(root.square() * d, n);
            }
        }

        #[test]
        fn batch_matches_single(v in prop::collection::vec(fp().prop_filter("nonzero", |x| !x.is_zero()), 0..40)) {
            let mut w = v.clone();
            crate::field::batch::invert(&mut w);
            for (x, y) in v.iter().zip(&w) {
                prop_assert_eq!(x.invert(), *y);
            }
        }
    }
}
