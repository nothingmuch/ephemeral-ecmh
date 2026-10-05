//! F_{p^2}, p = 2^64 - 59, i^2 = 2: `fp61x2`'s full-width counterpart, the
//! largest prime below 2^64, for about 2^128 points in the same two words.
//!
//! - Base field: `Fp`, one canonical u64 in [0, p). No bits are spare, so
//!   a sum carries out of its word, and a product reduces with 2^64 = 59:
//!   two folds by 59 and one conditional subtraction. Sums of up to three
//!   products reduce together, their carries folded with 2^128 = 59^2.
//! - Extension: `Fq` = a + b i. p = 5 mod 8, so -1 is a square and 2 is
//!   not: i^2 = 2 is irreducible, and multiplying by it is a doubling. A
//!   product is four base products and two reductions, a square three and
//!   two.
//! - Inversion: through the norm a^2 - 2b^2, inverted in F_p as x^(p - 2).
//! - Square root: the complex method, from Atkin's base square roots
//!   (p = 5 mod 8), with no inversion, and of a ratio n/d through
//!   n conj(d) / N(d), with none either.
//!
//! Curves: `curve::twisted64x2`.

use core::ops::{Add, Mul, Neg, Sub};

pub const P: u64 = 0u64.wrapping_sub(59);

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Fp(u64);

/// x + k 2^128 to [0, p), for k <= 2, with 2^64 = 59.
#[inline(always)]
fn reduce(x: u128, k: u64) -> Fp {
    // < 2^64 + 59 (2^64 + 2^65) < 2^72
    let t = (x as u64 as u128) + 59 * ((x >> 64) + ((k as u128) << 64));
    // < 2^64 + 59 * 2^8: at most one more carry, after which s is small
    let (s, c) = (t as u64).overflowing_add(59 * (t >> 64) as u64);
    let s = if c { s + 59 } else { s };
    Fp(if s >= P { s - P } else { s })
}

#[inline(always)]
fn add_wide(x: u128, y: u128) -> (u128, u64) {
    let (s, c) = x.overflowing_add(y);
    (s, c as u64)
}

impl Fp {
    pub const ZERO: Self = Self(0);
    pub const ONE: Self = Self(1);
    /// 1/2 = (p + 1)/2
    const HALF: Self = Self(P / 2 + 1);

    #[inline]
    pub const fn new(v: u64) -> Self {
        Self(if v >= P { v - P } else { v })
    }

    #[inline]
    pub const fn value(self) -> u64 {
        self.0
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

    /// x^(2^58 - 1), shared by inversion and square roots.
    fn pow_2_58_minus_1(x1: Self) -> Self {
        let x2 = x1.xsquare(1) * x1;
        let x4 = x2.xsquare(2) * x2;
        let x8 = x4.xsquare(4) * x4;
        let x16 = x8.xsquare(8) * x8;
        let x32 = x16.xsquare(16) * x16;
        let x48 = x32.xsquare(16) * x16;
        let x56 = x48.xsquare(8) * x8;
        x56.xsquare(2) * x2
    }

    /// x^(p - 2) = x^(64 (2^58 - 1) + 3); maps 0 to 0.
    pub fn invert(self) -> Self {
        Self::pow_2_58_minus_1(self).xsquare(6) * self.square() * self
    }

    /// x^((p - 5)/8) = x^(8 (2^58 - 1)).
    pub fn pow_p58(self) -> Self {
        Self::pow_2_58_minus_1(self).xsquare(3)
    }

    /// Atkin's square root for p = 5 mod 8, and its inverse for nonzero x:
    /// with u = (2x)^((p-5)/8) and v = 2x u^2, v^2 = -1 for a nonzero
    /// square, so x u (v - 1) is a root and u (v - 1) its inverse.
    /// Zero gives Some((0, 0)).
    fn sqrt_and_inverse(self) -> Option<(Self, Self)> {
        let x2 = self + self;
        let u = x2.pow_p58();
        let w = u * (x2 * u.square() - Self::ONE);
        let r = self * w;
        (r.square() == self).then_some((r, w))
    }

    /// Some square root if `self` is a square.
    pub fn sqrt(self) -> Option<Self> {
        self.sqrt_and_inverse().map(|(r, _)| r)
    }

    /// For n != 0, a root x of t/n if there is one, and for t != 0 its
    /// 1/(n x) = 1/sqrt(t n): x = t/sqrt(t n), from one Atkin root.
    fn sqrt_ratio(t: Self, n: Self) -> Option<(Self, Self)> {
        let (_, s) = (t * n).sqrt_and_inverse()?;
        Some((t * s, s))
    }
}

impl Add for Fp {
    type Output = Self;
    #[inline(always)]
    fn add(self, o: Self) -> Self {
        // with a carry, s + 59 = a + b - p, already canonical
        let (s, c) = self.0.overflowing_add(o.0);
        Self(if c {
            s + 59
        } else if s >= P {
            s - P
        } else {
            s
        })
    }
}

impl Sub for Fp {
    type Output = Self;
    #[inline(always)]
    fn sub(self, o: Self) -> Self {
        let (s, b) = self.0.overflowing_sub(o.0);
        Self(if b { s.wrapping_add(P) } else { s })
    }
}

impl Neg for Fp {
    type Output = Self;
    #[inline(always)]
    fn neg(self) -> Self {
        Self(if self.0 == 0 { 0 } else { P - self.0 })
    }
}

impl Mul for Fp {
    type Output = Self;
    #[inline(always)]
    fn mul(self, o: Self) -> Self {
        reduce(self.0 as u128 * o.0 as u128, 0)
    }
}

/// a + b i.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Fq {
    pub a: Fp,
    pub b: Fp,
}

impl Fq {
    pub const ZERO: Self = Self::new(Fp::ZERO, Fp::ZERO);
    pub const ONE: Self = Self::new(Fp::ONE, Fp::ZERO);

    #[inline]
    pub const fn new(a: Fp, b: Fp) -> Self {
        Self { a, b }
    }

    #[inline]
    pub fn is_zero(self) -> bool {
        self.a.is_zero() && self.b.is_zero()
    }

    /// (a + b i)(a - b i) = a^2 - 2b^2.
    #[inline]
    pub fn norm(self) -> Fp {
        self.a.square() - (self.b.square() + self.b.square())
    }

    /// (a^2 + 2b^2) + 2ab i: 3M.
    #[inline]
    pub fn square(self) -> Self {
        let (a, b) = (self.a.0 as u128, self.b.0 as u128);
        let (bb, ab) = (b * b, a * b);
        let (s, k1) = add_wide(a * a, bb);
        let (s, k2) = add_wide(s, bb);
        let (d, k) = add_wide(ab, ab);
        Self {
            a: reduce(s, k1 + k2),
            b: reduce(d, k),
        }
    }

    pub fn xsquare(self, n: u32) -> Self {
        (0..n).fold(self, |x, _| x.square())
    }

    /// (a - b i) / (a^2 - 2b^2); maps 0 to 0.
    pub fn invert(self) -> Self {
        let n = self.norm().invert();
        Self {
            a: self.a * n,
            b: -(self.b * n),
        }
    }

    /// Some square root if `self` is a square, equivalently if its norm
    /// is a square in F_p.
    pub fn sqrt(self) -> Option<Self> {
        Self::sqrt_over(self, Fp::ONE)
    }

    /// Some square root of n/d, if d != 0 and n/d is a square: of
    /// n conj(d) / N(d), with no inversion.
    pub fn sqrt_ratio(n: Self, d: Self) -> Option<Self> {
        let nd = d.norm();
        if nd.is_zero() {
            return None;
        }
        Self::sqrt_over(n * Self::new(d.a, -d.b), nd)
    }

    /// Some square root of w/n, n != 0 in F_p, by the complex method.
    ///
    /// (x + y i)^2 = a + b i is x^2 + 2y^2 = a, 2xy = b, and then the norm
    /// m = x^2 - 2y^2 is a square root of a^2 - 2b^2, so x^2 = (a + m)/2.
    /// With b != 0 the two signs of m give t = (a + m)/2 and (a - m)/2,
    /// whose product is b^2/2, a non-square: exactly one is a square, x^2.
    /// Then y = b/(2x). For a + b i = w/n: t is t'/n for t' = (w.a ± m')/2
    /// and m' a root of N(w), x = sqrt(t'/n), and y = w.b/(2 n x), which is
    /// w.b s/2 for s = 1/sqrt(t' n).
    fn sqrt_over(w: Self, n: Fp) -> Option<Self> {
        if w.b.is_zero() {
            // w.a/n, or w.a/(2n) as 2 is a non-square, is a square in F_p
            return Some(match Fp::sqrt_ratio(w.a, n) {
                Some((x, _)) => Self::new(x, Fp::ZERO),
                None => Self::new(Fp::ZERO, Fp::sqrt_ratio(w.a * Fp::HALF, n)?.0),
            });
        }
        let m = w.norm().sqrt()?;
        let (x, s) = Fp::sqrt_ratio((w.a + m) * Fp::HALF, n)
            .or_else(|| Fp::sqrt_ratio((w.a - m) * Fp::HALF, n))?;
        Some(Self::new(x, w.b * s * Fp::HALF))
    }
}

impl Add for Fq {
    type Output = Self;
    #[inline(always)]
    fn add(self, o: Self) -> Self {
        Self::new(self.a + o.a, self.b + o.b)
    }
}

impl Sub for Fq {
    type Output = Self;
    #[inline(always)]
    fn sub(self, o: Self) -> Self {
        Self::new(self.a - o.a, self.b - o.b)
    }
}

impl Neg for Fq {
    type Output = Self;
    #[inline(always)]
    fn neg(self) -> Self {
        Self::new(-self.a, -self.b)
    }
}

impl Mul for Fq {
    type Output = Self;
    /// (ac + 2bd) + (ad + bc) i: 4M, each half reduced once.
    #[inline(always)]
    fn mul(self, o: Self) -> Self {
        let (a, b) = (self.a.0 as u128, self.b.0 as u128);
        let (c, d) = (o.a.0 as u128, o.b.0 as u128);
        let bd = b * d;
        let (s, k1) = add_wide(a * c, bd);
        let (s, k2) = add_wide(s, bd);
        let (t, k) = add_wide(a * d, b * c);
        Self {
            a: reduce(s, k1 + k2),
            b: reduce(t, k),
        }
    }
}

impl crate::field::batch::Invert for Fq {
    const ONE: Self = Fq::ONE;
    fn inv(self) -> Self {
        self.invert()
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use proptest::prelude::*;

    pub fn fp() -> impl Strategy<Value = Fp> {
        prop_oneof![
            (0..P).prop_map(Fp),
            prop::sample::select(vec![0, 1, 2, 58, 59, 60, P - 1, P - 2, 1 << 63]).prop_map(Fp),
        ]
    }

    pub fn fq() -> impl Strategy<Value = Fq> {
        (fp(), fp()).prop_map(|(a, b)| Fq::new(a, b))
    }

    fn m(a: Fp, b: Fp) -> u64 {
        (a.0 as u128 * b.0 as u128 % P as u128) as u64
    }

    #[test]
    fn modulus() {
        assert_eq!(P as u128, (1 << 64) - 59);
        assert_eq!(P % 8, 5);
        assert_eq!(Fp::HALF + Fp::HALF, Fp::ONE);
    }

    proptest! {
        #[test]
        fn base_ops_match_reference(a in fp(), b in fp(), v in any::<u64>()) {
            prop_assert_eq!((a * b).0, m(a, b));
            prop_assert_eq!((a + b).0 as u128, (a.0 as u128 + b.0 as u128) % P as u128);
            prop_assert_eq!(a - b + b, a);
            prop_assert_eq!(Fp::new(v).0, v % P);
            prop_assert_eq!(a * a.invert(), if a.is_zero() { Fp::ZERO } else { Fp::ONE });
            match a.sqrt() {
                Some(s) => prop_assert_eq!(s.square(), a),
                // 2 is the non-square
                None => prop_assert!((a + a).sqrt().is_some()),
            }
        }

        #[test]
        fn reduce_matches_reference(x in any::<u128>(), k in 0u64..3) {
            // x + k 2^128 mod p, with 2^128 = 59^2
            let want = (x % P as u128 + k as u128 * 59 * 59) % P as u128;
            prop_assert_eq!(reduce(x, k).0 as u128, want);
        }

        #[test]
        fn mul_matches_reference(x in fq(), y in fq()) {
            let (a, b, c, d) = (x.a, x.b, y.a, y.b);
            let bd = Fp(m(b, d));
            let want = Fq::new(Fp(m(a, c)) + bd + bd, Fp(m(a, d)) + Fp(m(b, c)));
            prop_assert_eq!(x * y, want);
            prop_assert_eq!(x.square(), x * x);
            prop_assert_eq!(x.norm(), (x * Fq::new(x.a, -x.b)).a);
        }

        #[test]
        fn field_laws(x in fq(), y in fq(), z in fq()) {
            prop_assert_eq!(x * (y + z), x * y + x * z);
            prop_assert_eq!(x - y + y, x);
            prop_assert_eq!(x * x.invert(), if x.is_zero() { Fq::ZERO } else { Fq::ONE });
        }

        #[test]
        fn sqrt_exactly_on_squares(x in fq()) {
            let s = x.square().sqrt().unwrap();
            prop_assert!(s == x || s == -x);
            match x.sqrt() {
                Some(s) => prop_assert_eq!(s.square(), x),
                None => prop_assert!(x.norm().sqrt().is_none()),
            }
        }

        /// `sqrt_ratio` where n conj(d) is in F_p, here a N(d): a square,
        /// as every element of F_p is in F_{p^2}.
        #[test]
        fn sqrt_ratio_of_base_elements(a in fp(), d in fq()) {
            let n = Fq::new(a, Fp::ZERO) * d;
            match Fq::sqrt_ratio(n, d) {
                Some(s) => prop_assert_eq!(s.square() * d, n),
                None => prop_assert!(d.is_zero()),
            }
        }
    }

    #[test]
    fn i_squared_is_two() {
        let i = Fq::new(Fp::ZERO, Fp::ONE);
        let two = Fq::new(Fp(2), Fp::ZERO);
        assert_eq!(i * i, two);
        assert_eq!(i.square(), two);
        assert_eq!(two.sqrt().map(|s| s.square()), Some(two));
        // -1 is a square already in F_p
        assert!((-Fp::ONE).sqrt().is_some());
    }
}
