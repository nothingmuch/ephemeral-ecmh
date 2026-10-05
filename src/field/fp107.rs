//! F_p, p = 2^107 - 1: the field of `edwards107` and `weier107`, the
//! 14-byte counterpart of `fp127`.
//!
//! - Modulus: the Mersenne prime 2^107 - 1.
//! - Representation: `Fp`, one `u128` below 2^108 = 2p + 2, not necessarily
//!   canonical: 0, p and 2p all mean 0. The spare bits let every operation end
//!   with one Mersenne fold and skip the conditional subtraction `fp127`
//!   performs after each operation; `value`, `==`, parity and the Jacobi symbol
//!   canonicalize first.
//! - Backend: portable Rust, the same on every target. Multiplication
//!   splits into 54-bit limbs (both < 2^54 under that bound), so
//!   Karatsuba's middle product fits a u128: 3 64x64 products, not 4
//!   (`mul_schoolbook` keeps the 4-product version for comparison). sqrt
//!   is 105 squarings since (p + 1)/4 = 2^105.

use core::ops::{Add, Mul, Neg, Sub};

pub const P: u128 = (1 << 107) - 1;
pub const MASK107: u128 = P;
/// Encodings are 108 bits: a 107-bit coordinate and a sign bit.
pub const MASK108: u128 = (1 << 108) - 1;
/// Byte length of a point encoding.
pub const BYTES: usize = 14;
const M54: u64 = (1 << 54) - 1;
const M53: u128 = (1 << 53) - 1;

#[derive(Clone, Copy, Default)]
pub struct Fp(u128);

impl Fp {
    pub const ZERO: Self = Self(0);
    pub const ONE: Self = Self(1);

    /// Reduces any u128.
    #[inline]
    pub const fn new(v: u128) -> Self {
        Self::fold(v)
    }

    /// The canonical representative, in [0, p).
    #[inline]
    pub const fn value(self) -> u128 {
        // self.0 < 2^108, so s <= p + 1
        let s = (self.0 & P) + (self.0 >> 107);
        if s >= P { s - P } else { s }
    }

    /// v mod p, up to a multiple: below 2^108 for any v (v >> 107 < 2^21).
    #[inline]
    const fn fold(v: u128) -> Self {
        Self((v & P) + (v >> 107))
    }

    #[inline]
    pub fn is_zero(self) -> bool {
        self.value() == 0
    }

    #[inline]
    pub fn square(self) -> Self {
        let (a0, a1) = self.limbs();
        Self::combine(
            a0 as u128 * a0 as u128,
            (a0 << 1) as u128 * a1 as u128,
            a1 as u128 * a1 as u128,
        )
    }

    /// Schoolbook: 4 widening products to Karatsuba's 3 (`*`), but no
    /// pre-additions or subtractions on the path. Kept for the benches.
    pub fn mul_schoolbook(self, rhs: Self) -> Self {
        let ((a0, a1), (b0, b1)) = (self.limbs(), rhs.limbs());
        let w = |x: u64, y: u64| x as u128 * y as u128;
        Self::combine(w(a0, b0), w(a0, b1) + w(a1, b0), w(a1, b1))
    }

    pub fn xsquare(self, n: u32) -> Self {
        (0..n).fold(self, |x, _| x.square())
    }

    #[inline]
    fn limbs(self) -> (u64, u64) {
        (self.0 as u64 & M54, (self.0 >> 54) as u64)
    }

    /// lo + mid 2^54 + hi 2^108 with lo, hi < 2^108 and mid < 2^109.
    /// 2^107 = 1 and 2^108 = 2, so it is lo + 2 hi + (mid mod 2^53) 2^54 +
    /// mid >> 53 < 2^110, and one fold leaves it below 2^107 + 8.
    #[inline]
    fn combine(lo: u128, mid: u128, hi: u128) -> Self {
        Self::fold(lo + (hi << 1) + ((mid & M53) << 54) + (mid >> 53))
    }

    /// x^(2^105 - 1), shared by inversion and square-root ratios.
    fn pow_2k_minus_1(x1: Self) -> Self {
        let x2 = x1.xsquare(1) * x1;
        let x4 = x2.xsquare(2) * x2;
        let x8 = x4.xsquare(4) * x4;
        let x16 = x8.xsquare(8) * x8;
        let x32 = x16.xsquare(16) * x16;
        let x64 = x32.xsquare(32) * x32;
        let x96 = x64.xsquare(32) * x32;
        let x104 = x96.xsquare(8) * x8;
        x104.xsquare(1) * x1
    }

    /// x^(p-2) = x^(4 * (2^105 - 1) + 1); maps 0 to 0. 106 S + 10 M.
    pub fn invert(self) -> Self {
        Self::pow_2k_minus_1(self).xsquare(2) * self
    }

    /// x^((p-3)/4) = x^(2^105 - 1); see `fp127::Fp::pow_p34`.
    pub fn pow_p34(self) -> Self {
        Self::pow_2k_minus_1(self)
    }

    /// A square root of n/d, if d != 0 and n/d is a square; see
    /// `fp127::Fp::sqrt_ratio`.
    pub fn sqrt_ratio(n: Self, d: Self) -> Option<Self> {
        if d.is_zero() {
            return None;
        }
        let s = n * (n * d).pow_p34();
        (s.square() * d == n).then_some(s)
    }

    /// Quadratic character via the binary Jacobi symbol (variable time),
    /// as in `fp127`. 0 counts as a square.
    pub fn is_square(self) -> bool {
        let (mut a, mut n, mut t) = (self.value(), P, false);
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
        let y = self.xsquare(105);
        (y.square() == self).then_some(y)
    }

    #[inline]
    pub fn is_odd(self) -> bool {
        self.value() & 1 == 1
    }
}

impl PartialEq for Fp {
    #[inline]
    fn eq(&self, o: &Self) -> bool {
        self.value() == o.value()
    }
}

impl Eq for Fp {}

impl core::hash::Hash for Fp {
    fn hash<H: core::hash::Hasher>(&self, h: &mut H) {
        self.value().hash(h)
    }
}

impl core::fmt::Debug for Fp {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_tuple("Fp").field(&self.value()).finish()
    }
}

/// 4p = 2^109 - 4 exceeds every representative, so 4p - b doesn't wrap.
const P4: u128 = 4 * P;

impl Add for Fp {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self::fold(self.0 + rhs.0)
    }
}

impl Sub for Fp {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self::fold(self.0 + (P4 - rhs.0))
    }
}

impl Neg for Fp {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self::fold(P4 - self.0)
    }
}

impl Mul for Fp {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        let ((a0, a1), (b0, b1)) = (self.limbs(), rhs.limbs());
        let (lo, hi) = (a0 as u128 * b0 as u128, a1 as u128 * b1 as u128);
        // limb sums < 2^55: the product < 2^110 and it is at least lo + hi
        let mid = (a0 + a1) as u128 * (b0 + b1) as u128 - lo - hi;
        Self::combine(lo, mid, hi)
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

    /// Double-and-add reference, independent of the limb multiply.
    pub fn mul_ref(a: u128, b: u128) -> u128 {
        let mut acc = 0u128;
        for i in (0..107).rev() {
            acc = (acc << 1) % P;
            if b >> i & 1 == 1 {
                acc = (acc + a) % P;
            }
        }
        acc
    }

    /// Every representative below 2^108, not only canonical ones.
    pub fn fp() -> impl Strategy<Value = Fp> {
        prop_oneof![
            (0..P).prop_map(Fp),
            (0..1u128 << 108).prop_map(Fp),
            // Edges of the limb split, of the modulus and of the bound.
            prop::sample::select(vec![
                0,
                1,
                2,
                P - 1,
                P,
                P + 1,
                2 * P,
                2 * P + 1,
                (1 << 108) - 1,
                1 << 54,
                (1 << 54) - 1,
                1 << 106,
                1 << 107,
            ])
            .prop_map(Fp),
        ]
    }

    proptest! {
        #[test]
        fn bounded_after_ops(a in fp(), b in fp()) {
            for v in [a + b, a - b, a * b, a.square(), -a, Fp::new(a.0 << 20 | b.0)] {
                prop_assert!(v.0 < 1 << 108);
                prop_assert!(v.value() < P);
            }
        }

        #[test]
        fn new_reduces(v in any::<u128>()) {
            prop_assert_eq!(Fp::new(v).value(), v % P);
        }

        #[test]
        fn mul_matches_reference(a in fp(), b in fp()) {
            prop_assert_eq!((a * b).value(), mul_ref(a.value(), b.value()));
            prop_assert_eq!(a.square(), a * a);
            prop_assert_eq!(a.mul_schoolbook(b), a * b);
        }

        #[test]
        fn add_sub_neg(a in fp(), b in fp()) {
            prop_assert_eq!((a + b).value(), (a.value() + b.value()) % P);
            prop_assert_eq!((a - b).value(), (a.value() + P - b.value()) % P);
            prop_assert_eq!(a - b + b, a);
            prop_assert_eq!(a + -a, Fp::ZERO);
        }

        #[test]
        fn eq_is_mod_p(a in fp()) {
            prop_assert_eq!(Fp(a.value()), a);
            prop_assert_ne!(a + Fp::ONE, a);
            prop_assert_eq!(a.is_zero(), a.value() == 0);
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
