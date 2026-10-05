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
//! - `batch`: Montgomery's batch inversion, for any of them.

pub mod batch;
pub mod fp107;
pub mod fp127;
pub mod fp128;
pub mod gf2_109;
pub mod gf2_122;
pub mod gf2_127;
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

    #[test]
    fn boundaries_prime() {
        boundaries::<fp107::Fp>(1, fp107::P);
        boundaries::<fp127::Fp>(1, fp127::P);
        boundaries::<fp128::Fp>(1, fp128::P);
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
