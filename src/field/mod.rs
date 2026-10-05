//! The base fields the curves sit on, and batch inversion over them.
//!
//! - `gf2_127`: GF(2^127) = F_2\[z\]/(z^127 + z^63 + 1), crrl's, for
//!   `binary127`.
//! - `gf2_109`: GF(2^109) = F_2\[z\]/(z^109 + z^5 + z^4 + z^2 + 1),
//!   with a PMULL, a PCLMULQDQ and a portable backend, for `binary109`.
//! - `gf2_122`: GF(2^122) = GF(2^61)\[u\]/(u^2 + u + 1) over
//!   F_2\[z\]/(z^61 + z^23 + z^15 + z^5 + 1), with the same three
//!   backends, for `binary122`.
//! - `fp127`: F_p, p = 2^127 - 1, portable, for `edwards127` and `weier127`.
//! - `fp107`: F_p, p = 2^107 - 1, portable, for `edwards107` and `weier107`.
//! - `fp128`: F_p, p = 2^128 - 275, for the a = -1 Edwards curves of
//!   `twisted128` (p = 5 mod 8, so -1 is a square).
//! - `fp61x2`: F_{p^2}, p = 2^61 - 1, for `edwards61x2`, `weier61x2` and
//!   `twisted61x2`.
//! - `fp64x2`: F_{p^2}, p = 2^64 - 59, its full-width counterpart, for
//!   `twisted64x2`.
//! - `goldilocks2`: F_{p^2}, p = 2^64 - 2^32 + 1, on Plonky3's Goldilocks,
//!   for `twisted_goldilocks2`.
//! - `batch`: Montgomery's batch inversion, for any of them.

pub mod batch;
pub mod fp107;
pub mod fp127;
pub mod fp128;
pub mod fp61x2;
pub mod fp64x2;
pub mod gf2_109;
pub mod gf2_122;
pub mod gf2_127;
pub mod goldilocks2;

use core::fmt::Debug;
use core::ops::{Add, AddAssign, Neg, Sub};

/// A finite field, as every curve law uses it: representations may be
/// redundant (`fp107`, `gf2_127`, `gf2_109` and `gf2_122` reduce lazily),
/// so field equality is `equals`, and `batch::Invert::inv` maps 0 to 0.
pub trait Field: batch::Invert + Debug + Add<Output = Self> {
    const ZERO: Self;
    fn is_zero(self) -> bool;
    fn equals(self, o: Self) -> bool;
    fn square(self) -> Self;
}

/// A field of odd characteristic, F_p or an extension of it, as the odd
/// curve laws and their maps use it; `==` is `equals`.
pub trait OddField: Field + Eq + Sub<Output = Self> + Neg<Output = Self> {
    /// q, the number of elements.
    const ORDER: u128;
    /// A fixed non-square. Not by itself a map's constant: each map checks
    /// its own conditions on it.
    const NON_SQUARE: Self;
    /// Some square root if `self` is a square.
    fn sqrt(self) -> Option<Self>;
    /// 0 counts as a square.
    fn is_square(self) -> bool;
    /// Some square root of n/d, if d != 0 and n/d is a square. This default
    /// inverts d; fields with a faster exponentiation override it.
    fn sqrt_ratio(n: Self, d: Self) -> Option<Self> {
        if d.is_zero() {
            return None;
        }
        (n * d.inv()).sqrt()
    }
    /// RFC 9380's sgn0: the parity of the first nonzero canonical
    /// coefficient, lowest first. For x != 0, sgn0(-x) != sgn0(x).
    fn sgn0(self) -> bool;
}

/// F_p for an odd prime p.
pub trait Prime: OddField {
    /// Reduces v mod p.
    fn new(v: u128) -> Self;
}

/// An odd field's elements as distinct integers below 2^BITS: `pack` is
/// canonical, and `unpack` accepts exactly what `pack` returns. The
/// all-ones integer of BITS bits is not a valid packing, so codecs may reserve it.
pub trait Packed: OddField {
    const BITS: u32;
    fn pack(self) -> u128;
    fn unpack(v: u128) -> Option<Self>;
    /// Some element for every v below 2^BITS, `unpack`'s where that has
    /// one: hashes cut digests to field elements by it. Not uniform: the
    /// integers past each modulus fold onto the low residues.
    fn reduce(v: u128) -> Self;
}

/// `Field` for a type whose inherent methods have the trait's names, and
/// whose `==` is field equality; `#[inline]` so generic callers see through
/// the delegation across codegen units.
macro_rules! field {
    ($f:ty) => {
        impl crate::field::Field for $f {
            const ZERO: Self = <$f>::ZERO;
            #[inline]
            fn is_zero(self) -> bool {
                <$f>::is_zero(self)
            }
            #[inline]
            fn equals(self, o: Self) -> bool {
                self == o
            }
            #[inline]
            fn square(self) -> Self {
                <$f>::square(self)
            }
        }
    };
}

macro_rules! prime {
    ($m:ident, $non_square:expr, $bits:expr) => {
        field!($m::Fp);

        impl crate::field::OddField for $m::Fp {
            const ORDER: u128 = $m::P;
            const NON_SQUARE: Self = $non_square;
            #[inline]
            fn sqrt(self) -> Option<Self> {
                $m::Fp::sqrt(self)
            }
            #[inline]
            fn is_square(self) -> bool {
                $m::Fp::is_square(self)
            }
            #[inline]
            fn sqrt_ratio(n: Self, d: Self) -> Option<Self> {
                $m::Fp::sqrt_ratio(n, d)
            }
            #[inline]
            fn sgn0(self) -> bool {
                $m::Fp::is_odd(self)
            }
        }

        impl crate::field::Prime for $m::Fp {
            #[inline]
            fn new(v: u128) -> Self {
                $m::Fp::new(v)
            }
        }

        impl crate::field::Packed for $m::Fp {
            const BITS: u32 = $bits;
            #[inline]
            fn pack(self) -> u128 {
                $m::Fp::value(self)
            }
            #[inline]
            fn unpack(v: u128) -> Option<Self> {
                (v < $m::P).then(|| $m::Fp::new(v))
            }
            #[inline]
            fn reduce(v: u128) -> Self {
                $m::Fp::new(v)
            }
        }
    };
}
// p = 3 mod 4 for the Mersenne primes, and 5 mod 8 for 2^128 - 275.
prime!(fp107, fp107::Fp::new(fp107::P - 1), 107);
prime!(fp127, fp127::Fp::new(fp127::P - 1), 127);
prime!(fp128, fp128::Fp::new(2), 128);

/// `OddField` and `Packed` for a + b i in module `$m`, given a non-square
/// and the bits of p. Square testing uses the norm in F_p; the sign comes
/// from a, or from b when a = 0. Packing uses a | b << bits without the
/// spare bits used by the arithmetic.
macro_rules! quadratic {
    ($m:ident, $non_square:expr, $bits:expr) => {
        field!($m::Fq);

        impl crate::field::OddField for $m::Fq {
            const ORDER: u128 = $m::P as u128 * $m::P as u128;
            const NON_SQUARE: Self = $non_square;
            #[inline]
            fn sqrt(self) -> Option<Self> {
                $m::Fq::sqrt(self)
            }
            fn is_square(self) -> bool {
                self.norm().sqrt().is_some()
            }
            #[inline]
            fn sqrt_ratio(n: Self, d: Self) -> Option<Self> {
                $m::Fq::sqrt_ratio(n, d)
            }
            #[inline]
            fn sgn0(self) -> bool {
                let (a, b) = (self.a.value(), self.b.value());
                (if a != 0 { a } else { b }) & 1 == 1
            }
        }

        impl crate::field::Packed for $m::Fq {
            const BITS: u32 = 2 * $bits;
            #[inline]
            fn pack(self) -> u128 {
                self.a.value() as u128 | (self.b.value() as u128) << $bits
            }
            #[inline]
            fn unpack(v: u128) -> Option<Self> {
                let (a, b) = (v as u64 & u64::MAX >> (64 - $bits), (v >> $bits) as u64);
                (v >> $bits >> $bits == 0 && a < $m::P && b < $m::P)
                    .then(|| Self::new($m::Fp::new(a), $m::Fp::new(b)))
            }
            #[inline]
            fn reduce(v: u128) -> Self {
                let a = v as u64 & u64::MAX >> (64 - $bits);
                Self::new($m::Fp::new(a), $m::Fp::new((v >> $bits) as u64))
            }
        }
    };
}
// Norms 4^2 + 1 = 17 mod 2^61 - 1, and -2 mod 2^64 - 59: non-squares.
quadratic!(
    fp61x2,
    fp61x2::Fq::new(fp61x2::Fp::new(4), fp61x2::Fp::ONE),
    61
);
quadratic!(
    fp64x2,
    fp64x2::Fq::new(fp64x2::Fp::ZERO, fp64x2::Fp::ONE),
    64
);

impl Field for goldilocks2::Fq {
    const ZERO: Self = <Self as p3_field::PrimeCharacteristicRing>::ZERO;
    #[inline]
    fn is_zero(self) -> bool {
        self == Self::ZERO
    }
    #[inline]
    fn equals(self, o: Self) -> bool {
        self == o
    }
    #[inline]
    fn square(self) -> Self {
        p3_field::PrimeCharacteristicRing::square(&self)
    }
}

impl OddField for goldilocks2::Fq {
    const ORDER: u128 = goldilocks2::P as u128 * goldilocks2::P as u128;
    /// i: its norm, -7, is a non-square, as 7 is and -1 isn't.
    const NON_SQUARE: Self = goldilocks2::new(goldilocks2::Fp::new(0), goldilocks2::Fp::new(1));
    #[inline]
    fn sqrt(self) -> Option<Self> {
        goldilocks2::sqrt(self)
    }
    fn is_square(self) -> bool {
        goldilocks2::sqrt_base(goldilocks2::norm(self)).is_some()
    }
    #[inline]
    fn sqrt_ratio(n: Self, d: Self) -> Option<Self> {
        goldilocks2::sqrt_ratio(n, d)
    }
    #[inline]
    fn sgn0(self) -> bool {
        use p3_field::PrimeField64;
        let (a, b) = goldilocks2::parts(self);
        let (a, b) = (a.as_canonical_u64(), b.as_canonical_u64());
        (if a != 0 { a } else { b }) & 1 == 1
    }
}

impl Packed for goldilocks2::Fq {
    const BITS: u32 = 128;
    #[inline]
    fn pack(self) -> u128 {
        use p3_field::PrimeField64;
        let (a, b) = goldilocks2::parts(self);
        a.as_canonical_u64() as u128 | (b.as_canonical_u64() as u128) << 64
    }
    #[inline]
    fn unpack(v: u128) -> Option<Self> {
        use goldilocks2::{Fp, P, new};
        let (a, b) = (v as u64, (v >> 64) as u64);
        (a < P && b < P).then(|| new(Fp::new(a), Fp::new(b)))
    }
    #[inline]
    fn reduce(v: u128) -> Self {
        use goldilocks2::{Fp, new};
        new(Fp::new(v as u64), Fp::new((v >> 64) as u64))
    }
}

/// GF(2^m), as the binary curve laws use it: the arithmetic, and a
/// canonical integer for each element. + is also -.
pub trait Binary: Field + AddAssign {
    /// The bits a canonical integer may have set, at most m of them.
    const MASK: u128;
    /// The element whose integer is v's bits in `MASK`.
    fn new(v: u128) -> Self;
    /// The canonical integer, in `MASK`: `new(value(x))` is x.
    fn value(self) -> u128;
    /// The square root, which every element has.
    fn sqrt(self) -> Self;
    /// Tr(self), 0 or 1.
    fn trace(self) -> u32;
    /// A root w of w^2 + w = self, given Tr(self) = 0; w + 1 is the other.
    fn solve(self) -> Self;
}

/// Unlike field!, routes is_zero and equals through the module's canonicalizing
/// is_zero/eq, since representatives are reduced lazily.
macro_rules! binary {
    ($m:ident, $mask:expr, $solve:path) => {
        impl crate::field::Field for $m::Gf {
            const ZERO: Self = $m::Gf::ZERO;
            #[inline]
            fn is_zero(self) -> bool {
                $m::is_zero(self)
            }
            #[inline]
            fn equals(self, o: Self) -> bool {
                $m::eq(self, o)
            }
            #[inline]
            fn square(self) -> Self {
                $m::Gf::square(self)
            }
        }

        impl crate::field::Binary for $m::Gf {
            const MASK: u128 = $mask;
            #[inline]
            fn new(v: u128) -> Self {
                $m::from_u128(v)
            }
            #[inline]
            fn value(self) -> u128 {
                $m::to_u128(self)
            }
            #[inline]
            fn sqrt(self) -> Self {
                $m::Gf::sqrt(self)
            }
            #[inline]
            fn trace(self) -> u32 {
                $m::Gf::trace(self)
            }
            #[inline]
            fn solve(self) -> Self {
                $solve(self)
            }
        }
    };
}
// For odd m the half-trace solves it: H(c)^2 + H(c) = c + Tr(c).
binary!(gf2_127, gf2_127::MASK127, gf2_127::halftrace8);
binary!(gf2_109, gf2_109::MASK109, gf2_109::halftrace8);
// For m = 122, Tr(1) = 0 and the half-trace is undefined: Pornin's QSolve.
binary!(gf2_122, gf2_122::MASK122, gf2_122::Gf::qsolve);

/// Half-trace window tables for a binary field of degree M, from
/// h\[i\] = H(z^i): `T[j][b]` = H(sum of z^(W j + k) over the bits k of b),
/// for W-bit windows. Entries for bits past z^(M - 1) are never read, since
/// the input is canonical.
pub(crate) const fn window_tables<const M: usize, const N: usize, const S: usize>(
    h: &[u128; M],
) -> [[u128; S]; N] {
    let w = S.trailing_zeros() as usize;
    let mut t = [[0; S]; N];
    let mut j = 0;
    while j < N {
        let mut b = 1;
        while b < S {
            // b without its lowest bit is already done
            let i = w * j + b.trailing_zeros() as usize;
            t[j][b] = t[j][b & (b - 1)] ^ if i < M { h[i] } else { 0 };
            b += 1;
        }
        j += 1;
    }
    t
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn pow<F: Field>(x: F, e: u128) -> F {
        (0..128 - e.leading_zeros()).rev().fold(F::ONE, |r, i| {
            if e >> i & 1 == 1 {
                r.square() * x
            } else {
                r.square()
            }
        })
    }

    /// Shared field and codec laws; k coefficients occupy BITS/k bits each.
    fn odd<F: Packed>(k: u32, x: F, y: F, v: u128) -> Result<(), TestCaseError> {
        prop_assert!(pow(x, F::ORDER).equals(x));
        prop_assert!(!F::NON_SQUARE.is_square());
        prop_assert_eq!(x.is_square(), x.sqrt().is_some());
        if let Some(s) = x.sqrt() {
            prop_assert_eq!(s.square(), x);
        }
        prop_assert!((x * x).is_square());
        prop_assert!(F::sqrt_ratio(x, F::ZERO).is_none());
        if !y.is_zero() {
            let r = F::sqrt_ratio(x, y);
            prop_assert_eq!(r.is_some(), (x * y.inv()).is_square());
            if let Some(s) = r {
                prop_assert_eq!(s.square() * y, x);
            }
        }
        let w = F::BITS / k;
        let c = x.pack();
        prop_assert_eq!(c >> (F::BITS - 1) >> 1, 0);
        let first = (0..k)
            .map(|i| c >> (w * i) & (u128::MAX >> (128 - w)))
            .find(|&a| a != 0);
        prop_assert_eq!(x.sgn0(), first.is_some_and(|a| a & 1 == 1));
        if !x.is_zero() {
            prop_assert_ne!((-x).sgn0(), x.sgn0());
        }
        prop_assert_eq!(F::unpack(c), Some(x));
        let ones = u128::MAX >> (128 - F::BITS);
        prop_assert!(F::unpack(ones).is_none());
        let v = v & ones;
        if let Some(z) = F::unpack(v) {
            prop_assert_eq!(z.pack(), v);
            prop_assert_eq!(F::reduce(v), z);
        }
        Ok(())
    }

    /// `new` keeps the bits of `MASK`, and `solve` finds a root on its
    /// domain Tr(c) = 0, which every c^2 + c is in.
    fn solve<F: Binary>(v: u128) -> Result<(), TestCaseError> {
        let c = F::new(v);
        prop_assert_eq!(c.value(), v & F::MASK);
        let c = c.square() + c;
        let w = c.solve();
        prop_assert!((w.square() + w).equals(c));
        Ok(())
    }

    /// Zero denominators and packing boundaries can be missed by round trips.
    /// The field has k coefficients of BITS/k bits each.
    fn boundaries<F: Packed>(k: u32, p: u128) {
        let (zero, one) = (F::ZERO, F::ONE);
        assert!(zero.inv().is_zero());
        assert!(F::sqrt_ratio(zero, zero).is_none());
        assert!(F::sqrt_ratio(one, zero).is_none());
        assert_eq!(F::sqrt_ratio(zero, one), Some(zero));
        assert_eq!(one.pack(), 1);
        let w = F::BITS / k;
        let top = (0..k).fold(0, |c, i| c | (p - 1) << (w * i));
        assert_eq!(F::unpack(top).map(F::pack), Some(top));
        for i in 0..k {
            assert!(F::unpack(p << (w * i)).is_none());
        }
        for b in F::BITS..128 {
            assert!(F::unpack(1 << b).is_none());
        }
        let ones = u128::MAX >> (128 - w);
        let folded = (0..k).fold(0, |c, i| c | (ones % p) << (w * i));
        let all = (0..k).fold(0, |c, i| c | ones << (w * i));
        assert_eq!(F::reduce(all).pack(), folded);
    }

    /// i packs as 1 << w, the largest packing exceeds q, and sgn0 reads a
    /// before b, whatever b's parity.
    fn quadratic<F: Packed>(p: u128, i: F) {
        boundaries::<F>(2, p);
        let w = F::BITS / 2;
        assert_eq!(i.pack(), 1 << w);
        assert!((p - 1) | (p - 1) << w >= F::ORDER);
        let sgn0 = |a: u128, b: u128| F::unpack(a | b << w).unwrap().sgn0();
        assert!(sgn0(0, 1) && !sgn0(0, 2) && sgn0(1, 2) && !sgn0(2, 1));
    }

    #[test]
    fn boundaries_prime() {
        boundaries::<fp107::Fp>(1, fp107::P);
        boundaries::<fp127::Fp>(1, fp127::P);
        boundaries::<fp128::Fp>(1, fp128::P);
    }

    #[test]
    fn boundaries_quadratic() {
        use fp61x2 as a;
        use fp64x2 as b;
        use goldilocks2 as g;
        quadratic(a::P as u128, a::Fq::new(a::Fp::ZERO, a::Fp::ONE));
        quadratic(b::P as u128, b::Fq::new(b::Fp::ZERO, b::Fp::ONE));
        quadratic(g::P as u128, g::new(g::Fp::new(0), g::Fp::new(1)));
    }

    /// The Goldilocks words p and p + 1 are 0 and 1, and p is a zero
    /// denominator.
    #[test]
    fn goldilocks2_words_p_and_p_plus_1() {
        use goldilocks2::{Fp, P, new};
        let (z, o) = (new(Fp::new(P), Fp::new(P)), new(Fp::new(P + 1), Fp::new(P)));
        assert!(z.is_zero() && z.pack() == 0 && !z.sgn0());
        assert!(batch::Invert::inv(z).is_zero() && OddField::sqrt_ratio(o, z).is_none());
        assert!(o.equals(<goldilocks2::Fq as batch::Invert>::ONE) && o.pack() == 1 && o.sgn0());
    }

    proptest! {
        #[test]
        fn odd_fp107(x in fp107::tests::fp(), y in fp107::tests::fp(), v in any::<u128>()) {
            odd(1, x, y, v)?;
        }

        #[test]
        fn odd_fp127(x in fp127::tests::fp(), y in fp127::tests::fp(), v in any::<u128>()) {
            odd(1, x, y, v)?;
        }

        #[test]
        fn odd_fp128(x in fp128::tests::fp(), y in fp128::tests::fp(), v in any::<u128>()) {
            odd(1, x, y, v)?;
        }

        #[test]
        fn odd_fp61x2(x in fp61x2::tests::fq(), y in fp61x2::tests::fq(), v in any::<u128>()) {
            odd(2, x, y, v)?;
        }

        #[test]
        fn odd_fp64x2(x in fp64x2::tests::fq(), y in fp64x2::tests::fq(), v in any::<u128>()) {
            odd(2, x, y, v)?;
        }

        #[test]
        fn odd_goldilocks2(x in goldilocks2::tests::fq(), y in goldilocks2::tests::fq(), v in any::<u128>()) {
            odd(2, x, y, v)?;
        }

        /// Goldilocks words may be p + a: field-equal to a, signed,
        /// packed and rooted as a.
        #[test]
        fn goldilocks2_redundant_words(a in 0..=u64::MAX - goldilocks2::P, b in any::<u64>()) {
            use goldilocks2::{Fp, P, new};
            let b = Fp::new(b);
            let (x, y) = (new(Fp::new(P + a), b), new(Fp::new(a), b));
            prop_assert!(x.equals(y));
            prop_assert_eq!(x.sgn0(), y.sgn0());
            prop_assert_eq!(x.pack(), y.pack());
            let one = <goldilocks2::Fq as batch::Invert>::ONE;
            let r = OddField::sqrt_ratio(x, x).map(|s| s.square());
            prop_assert_eq!(r, (!y.is_zero()).then_some(one));
            prop_assert_eq!(OddField::sqrt_ratio(one, x).is_some(), OddField::sqrt_ratio(one, y).is_some());
        }

        #[test]
        fn solve_gf2_127(v in any::<u128>()) {
            solve::<gf2_127::Gf>(v)?;
        }

        #[test]
        fn solve_gf2_109(v in any::<u128>()) {
            solve::<gf2_109::Gf>(v)?;
        }

        #[test]
        fn solve_gf2_122(v in any::<u128>()) {
            solve::<gf2_122::Gf>(v)?;
        }
    }
}
