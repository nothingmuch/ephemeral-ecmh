//! The base fields the curves sit on, and batch inversion over them.
//!
//! - `gf2_127`: GF(2^127) = F_2\[z\]/(z^127 + z^63 + 1), crrl's, for
//!   `binary127`.
//! - `batch`: Montgomery's batch inversion, for any of them.

pub mod batch;
pub mod gf2_127;
use core::fmt::Debug;
use core::ops::{Add, AddAssign};

/// A finite field, as every curve law uses it: representations may be
/// redundant (`fp107`, `gf2_127`, `gf2_109` and `gf2_122` reduce lazily),
/// so field equality is `equals`, and `batch::Invert::inv` maps 0 to 0.
pub trait Field: batch::Invert + Debug + Add<Output = Self> {
    const ZERO: Self;
    fn is_zero(self) -> bool;
    fn equals(self, o: Self) -> bool;
    fn square(self) -> Self;
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

    proptest! {
        #[test]
        fn solve_gf2_127(v in any::<u128>()) {
            solve::<gf2_127::Gf>(v)?;
        }
    }
}
