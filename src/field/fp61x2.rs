//! F_{p^2}, p = 2^61 - 1, i^2 = -1: a quadratic extension of a 61-bit
//! Mersenne prime, supporting curves of about 2^122 points in two
//! 64-bit words.
//!
//! - Base field: `Fp`, one canonical u64 in [0, p). The three bits above
//!   p leave room to multiply unreduced sums: anything below 2^124 reduces
//!   with two 61-bit folds and one conditional subtraction.
//! - Extension: `Fq` = a + b i, with i^2 = -1 irreducible as p = 3 mod 4.
//!   A product is four base products and two reductions, a square two
//!   and two.
//! - Inversion: through the norm a^2 + b^2, inverted in F_p as x^(p - 2).
//! - Square root: the complex method, from base square roots x^((p+1)/4)
//!   = x^(2^59), with no inversion, and of a ratio n/d through
//!   n conj(d) / N(d), with none either.
//!
//! Curves: `curve::edwards61x2`, `curve::weier61x2` and
//! `curve::twisted61x2`.

use core::ops::{Add, Mul, Neg, Sub};

pub const P: u64 = (1 << 61) - 1;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Fp(u64);

/// x < 2^124 to [0, p), with 2^61 = 1.
#[inline(always)]
const fn reduce(x: u128) -> Fp {
    let s = (x as u64 & P) + (x >> 61) as u64;
    let s = (s & P) + (s >> 61);
    Fp(if s >= P { s - P } else { s })
}

impl Fp {
    pub const ZERO: Self = Self(0);
    pub const ONE: Self = Self(1);
    /// 1/2 = (p + 1)/2
    const HALF: Self = Self(1 << 60);

    #[inline]
    pub const fn new(v: u64) -> Self {
        reduce(v as u128)
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

    /// x^((p - 3)/4) = x^(2^59 - 1).
    pub fn pow_p34(self) -> Self {
        let x1 = self;
        let x2 = x1.xsquare(1) * x1;
        let x4 = x2.xsquare(2) * x2;
        let x8 = x4.xsquare(4) * x4;
        let x16 = x8.xsquare(8) * x8;
        let x32 = x16.xsquare(16) * x16;
        let x48 = x32.xsquare(16) * x16;
        let x56 = x48.xsquare(8) * x8;
        let x58 = x56.xsquare(2) * x2;
        x58.xsquare(1) * x1
    }

    /// x^(p - 2) = x^(4 (2^59 - 1) + 1); maps 0 to 0.
    pub fn invert(self) -> Self {
        self.pow_p34().xsquare(2) * self
    }

    /// Some square root if `self` is a square.
    pub fn sqrt(self) -> Option<Self> {
        let y = self.xsquare(59);
        (y.square() == self).then_some(y)
    }

    /// For n != 0, a root x of t/n if there is one, and for t != 0 its
    /// 1/(n x) = 1/sqrt(t n): s = (t n)^((p - 3)/4) and x = t s. One
    /// exponentiation.
    fn sqrt_ratio(t: Self, n: Self) -> Option<(Self, Self)> {
        let s = (t * n).pow_p34();
        let x = t * s;
        (x.square() * n == t).then_some((x, s))
    }
}

impl Add for Fp {
    type Output = Self;
    #[inline(always)]
    fn add(self, o: Self) -> Self {
        let s = self.0 + o.0;
        Self(if s >= P { s - P } else { s })
    }
}

impl Sub for Fp {
    type Output = Self;
    #[inline(always)]
    fn sub(self, o: Self) -> Self {
        self + -o
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
        reduce(self.0 as u128 * o.0 as u128)
    }
}

/// a + b i.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct Fq {
    pub a: Fp,
    pub b: Fp,
}

/// p 2^62, a multiple of p above any product of two canonical elements.
const PP: u128 = (P as u128) << 62;

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

    /// (a + b i)(a - b i) = a^2 + b^2.
    #[inline]
    pub fn norm(self) -> Fp {
        let (a, b) = (self.a.0 as u128, self.b.0 as u128);
        reduce(a * a + b * b)
    }

    /// (a + b)(a - b) + 2ab i, the difference unreduced: 2M.
    #[inline]
    pub fn square(self) -> Self {
        let (a, b) = (self.a.0 as u128, self.b.0 as u128);
        Self {
            a: reduce((a + b) * (a + P as u128 - b)),
            b: reduce(2 * a * b),
        }
    }

    pub fn xsquare(self, n: u32) -> Self {
        (0..n).fold(self, |x, _| x.square())
    }

    /// (a - b i) / (a^2 + b^2); maps 0 to 0.
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
    /// (x + y i)^2 = a + b i is x^2 - y^2 = a, 2xy = b, and then the norm
    /// m = x^2 + y^2 is a square root of a^2 + b^2, so x^2 = (a + m)/2.
    /// With b != 0 the two signs of m give t = (a + m)/2 and (a - m)/2,
    /// whose product is -b^2/4, a non-square: exactly one is a square, x^2.
    /// Then y = b/(2x). For a + b i = w/n: t is t'/n for t' = (w.a ± m')/2
    /// and m' a root of N(w), x = sqrt(t'/n), and y = w.b/(2 n x), which is
    /// w.b s/2 for s = 1/sqrt(t' n).
    fn sqrt_over(w: Self, n: Fp) -> Option<Self> {
        if w.b.is_zero() {
            // w.a/n, or -w.a/n as -1 is a non-square, is a square in F_p
            return Some(match Fp::sqrt_ratio(w.a, n) {
                Some((x, _)) => Self::new(x, Fp::ZERO),
                None => Self::new(Fp::ZERO, Fp::sqrt_ratio(-w.a, n)?.0),
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
    /// (ac - bd) + (ad + bc) i: 4M, each half reduced once.
    #[inline(always)]
    fn mul(self, o: Self) -> Self {
        let (a, b) = (self.a.0 as u128, self.b.0 as u128);
        let (c, d) = (o.a.0 as u128, o.b.0 as u128);
        Self {
            a: reduce(a * c + PP - b * d),
            b: reduce(a * d + b * c),
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
            prop::sample::select(vec![0, 1, 2, P - 1, P - 2, 1 << 60]).prop_map(Fp),
        ]
    }

    pub fn fq() -> impl Strategy<Value = Fq> {
        (fp(), fp()).prop_map(|(a, b)| Fq::new(a, b))
    }

    fn m(a: Fp, b: Fp) -> u64 {
        (a.0 as u128 * b.0 as u128 % P as u128) as u64
    }

    proptest! {
        #[test]
        fn base_ops_match_reference(a in fp(), b in fp(), v in any::<u64>()) {
            prop_assert_eq!((a * b).0, m(a, b));
            prop_assert_eq!((a + b).0, (a.0 + b.0) % P);
            prop_assert_eq!(a - b + b, a);
            prop_assert_eq!(Fp::new(v).0, v % P);
            prop_assert_eq!(a * a.invert(), if a.is_zero() { Fp::ZERO } else { Fp::ONE });
            if let Some(s) = a.sqrt() {
                prop_assert_eq!(s.square(), a);
            } else {
                prop_assert!((-a).sqrt().is_some());
            }
        }

        #[test]
        fn mul_matches_reference(x in fq(), y in fq()) {
            let (a, b, c, d) = (x.a, x.b, y.a, y.b);
            let want = Fq::new(Fp(m(a, c)) - Fp(m(b, d)), Fp(m(a, d)) + Fp(m(b, c)));
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
    fn i_squared_is_minus_one() {
        let i = Fq::new(Fp::ZERO, Fp::ONE);
        assert_eq!(i * i, -Fq::ONE);
        assert_eq!(i.square(), -Fq::ONE);
        // and a non-square base element has its root on the i axis
        assert_eq!(
            Fq::new(-Fp::ONE, Fp::ZERO).sqrt().map(|s| s.square()),
            Some(-Fq::ONE)
        );
    }
}
