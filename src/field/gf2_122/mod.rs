//! GF(2^122) = GF(2^61)\[u\]/(u^2 + u + 1): Pornin's GLS254 field shape
//! (ePrint 2023/1688, GF(2^254) = GF(2^127)\[u\]) at half the size.
//!
//! # Why this field
//!
//! With m = 61 a base element is one machine word, so a GF(2^122) element
//! is 2 words and an (X : S : Z : T) point is 8, as over GF(2^127).
//! Pornin's complete formulas carry over. A multiplication is 3 one-word
//! carry-less products (Karatsuba over the base field; crrl's GF(2^127)
//! uses 4, schoolbook), a multiplication by a constant in GF(2^61) is 2,
//! and one by u is additions only.
//!
//! 61 is prime, so the descents are to GF(2^61) with n = 2 (genus <= 2,
//! index calculus about as costly as rho) and to GF(4) and F_2; but rho is
//! about 2^60.3 (2^59.8 for a GLS curve, whose endomorphism reduces the
//! rho work factor by sqrt 2) against GF(2^127)'s 2^62.8, and the base
//! modulus is a pentanomial, reduced in two folds, where GF(2^127)'s
//! trinomial folds once. benches/compare122.rs compares the two.
//!
//! # Construction
//!
//! - Base field ([`gf2_61`]): F_2\[z\]/(f), f = z^61 + z^23 + z^15 + z^5 + 1.
//!   - There is no irreducible trinomial of degree 61 (Swan: 61 = 5 mod 8).
//!   - Elements are one u64 reduced only mod z^3 f = z^64 + rho,
//!     rho = z^26 + z^18 + z^8 + z^3, so products fold at word boundaries
//!     in two steps. Any f with k3 <= 28 folds at the same cost (4 terms);
//!     the lowest is z^61 + z^5 + z^2 + z + 1.
//!   - This f is the lowest with all middle exponents odd
//!     (sage/gf2_122.sage). That gives Tr(v) = v_0 on canonical elements
//!     (v_0 + v_61 on any word), where the lowest f's trace is 3 bits, and a
//!     4-term sqrt(z) = z^31 + z^12 + z^8 + z^3, so a square root is a bit
//!     unshuffle plus 4 shifts where the lowest f needs a multiplication.
//! - Tower: u^2 + u + 1 is irreducible over GF(2^61) because 61 is odd
//!   (a root would solve x^2 + x = 1, which needs Tr(1) = 0). An element is
//!   a0 + a1 u, stored as the two words (a0, a1), each any representative.
//! - Canonical form: [`to_u128`] = a0 | a1 << 64 with both halves reduced
//!   (bits 61..64 and 125..128 are zero); [`to_bytes`] is its 16
//!   little-endian bytes.
//!
//! # Costs
//!
//! M61 and S61 include base-field reduction; on the hardware backends,
//! each uses one carry-less product and two folds. H61 is the base half-trace
//! (8 table lookups).
//!
//! | operation | cost | how |
//! |---|---|---|
//! | mul | 3 carry-less products, 2 base reductions | Karatsuba: c0 = a0b0 + a1b1, c1 = (a0+a1)(b0+b1) + a0b0; the three 128-bit products are combined unreduced |
//! | square | 2 S61 | c0 = (a0+a1)^2, c1 = a1^2 |
//! | mul by u | 1 add | (a1, a0+a1) |
//! | sqrt | 2 base sqrt (shifts) | b1 = sqrt(a1), b0 = sqrt(a0+a1) |
//! | invert | 3 M61 + 1 S61 + 1 I61 | N = a0(a0+a1) + a1^2, 1/a = ((a0+a1) + a1 u)/N; I61 = 7 M61 + 60 S61 |
//! | trace | 1 base trace (2 bits) | Tr_122(a) = Tr_61(a1) |
//! | qsolve (z^2 + z = d) | 2 H61 | see [`Gf::qsolve`] |
//!
//! # Backends
//!
//! Picked at compile time, as in `gf2_109`:
//! - `pmull`: aarch64 with `target_feature = "aes"` (PMULL, every Apple
//!   core). An element is one NEON register.
//! - `pclmul`: x86_64 with `target_feature = "pclmulqdq"` (build with
//!   `-C target-cpu=native` or `+pclmulqdq`). One SSE register.
//! - `soft`: everything else: integer-multiply carry-less products and shift
//!   folds. Tests compare it with the selected backend and bit-serial reference.
//!
//! Carry-less multiply instructions (64 x 64 -> 128) per GF(2^122)
//! multiplication in `pmull` and `pclmul`: 3 for the products and 4 for the
//! two reductions (2 per word) = 7, against crrl's GFb127 at 4 + 0 (its
//! trinomial folds with 3 shifts). Squaring: 2 + 4, against 2 + 0.
//!
//! # What a curve over this field must check
//!
//! For E : y^2 + xy = x^3 + a x^2 + B over GF(2^122), a in {1, u}:
//! - If both a and B lie in GF(2^61) (a = 1, B = b0), E is defined over the
//!   subfield and #E(GF(2^61)) divides #E(GF(2^122)): never 2 * prime.
//!   With a = u (Tr_122(u) = 1), B in GF(2^61) is GLS's quadratic twist of
//!   the subfield curve, of order (q - 1)^2 + t^2, which can be 2r; it has
//!   the GLS endomorphism psi, which lowers its rho work factor by sqrt 2.
//! - #E = 2r with r prime, inside the Hasse bound; E\[r\] = 2E is then
//!   {P : Tr(x(P)) = Tr(a)} + O.
//! - Embedding degree greater than 2^20 (MOV / Frey-Rück). The order
//!   policy excludes anomalous curves (#E = 2^122), since #E = 2r with r
//!   odd.
//! - Weil descent: the subfields are GF(2^61), GF(4) and F_2, so the
//!   descents are of degree n = 2 (genus <= 2, index calculus about
//!   q = rho; a prime m leaves no other intermediate field), 61 and 122
//!   (genus 2^(mb - 1) for the magic number mb); check mb as sage/ghs.sage
//!   does for prime m.

use core::fmt;
use core::ops::{Add, AddAssign, Div, Mul, MulAssign};

pub mod gf2_61;
#[cfg(test)]
mod kats;
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
pub use gf2_61::Gf61;
use gf2_61::MASK61;

/// The bits of a canonical [`to_u128`]: a0 in 0..61, a1 in 64..125.
pub const MASK122: u128 = (MASK61 as u128) << 64 | MASK61 as u128;

/// Canonical form a0 | a1 << 64 (bit i of each word = coefficient of z^i).
#[inline]
pub fn to_u128(x: Gf) -> u128 {
    let (a0, a1) = x.parts();
    gf2_61::to_u64(a0) as u128 | (gf2_61::to_u64(a1) as u128) << 64
}

/// Bits outside [`MASK122`] are dropped (not reduced), so the map is onto
/// and exact.
#[inline]
pub fn from_u128(v: u128) -> Gf {
    let v = v & MASK122;
    Gf::w64le(v as u64, (v >> 64) as u64)
}

/// The canonical 16-byte encoding: [`to_u128`], little-endian.
pub fn to_bytes(x: Gf) -> [u8; 16] {
    to_u128(x).to_le_bytes()
}

/// Inverse of [`to_bytes`]; rejects bits outside [`MASK122`].
pub fn from_bytes(b: [u8; 16]) -> Option<Gf> {
    let v = u128::from_le_bytes(b);
    (v & !MASK122 == 0).then(|| from_u128(v))
}

#[inline]
pub fn is_zero(x: Gf) -> bool {
    to_u128(x) == 0
}

#[inline]
pub fn eq(a: Gf, b: Gf) -> bool {
    is_zero(a + b)
}

impl Gf {
    pub const ZERO: Self = Self::w64le(0, 0);
    pub const ONE: Self = Self::w64le(1, 0);
    pub const U: Self = Self::w64le(0, 1);

    /// a0 + a1 u.
    #[inline(always)]
    pub fn from_parts(a0: Gf61, a1: Gf61) -> Self {
        Self::w64le(a0.limb(), a1.limb())
    }

    /// (a0, a1) for a0 + a1 u.
    #[inline(always)]
    pub fn parts(self) -> (Gf61, Gf61) {
        let [a0, a1] = self.limbs();
        (Gf61::w64(a0), Gf61::w64(a1))
    }

    /// The subfield element c, embedded.
    #[inline(always)]
    pub fn from_base(c: Gf61) -> Self {
        Self::w64le(c.limb(), 0)
    }

    #[inline(always)]
    pub fn square(self) -> Self {
        backend::square(self)
    }

    pub fn xsquare(self, n: u32) -> Self {
        (0..n).fold(self, |x, _| x.square())
    }

    /// a u = a1 + (a0 + a1) u, as u^2 = u + 1.
    #[inline(always)]
    pub fn mul_u(self) -> Self {
        let [a0, a1] = self.limbs();
        Self::w64le(a1, a0 ^ a1)
    }

    /// a c for c in GF(2^61): 2 M61.
    #[inline(always)]
    pub fn mul_base(self, c: Gf61) -> Self {
        let (a0, a1) = self.parts();
        Self::from_parts(a0 * c, a1 * c)
    }

    /// The conjugate a^(2^61) = (a0 + a1) + a1 u (u -> u + 1).
    #[inline(always)]
    pub fn conjugate(self) -> Self {
        let [a0, a1] = self.limbs();
        Self::w64le(a0 ^ a1, a1)
    }

    /// N(a) = a a^(2^61) = a0^2 + a0 a1 + a1^2 = a0 (a0 + a1) + a1^2.
    #[inline(always)]
    pub fn norm(self) -> Gf61 {
        let (a0, a1) = self.parts();
        a0 * (a0 + a1) + a1.square()
    }

    /// 1/a = a^(2^61) / N(a): 3 M61 + 1 S61 + one GF(2^61) inversion. 0 -> 0.
    pub fn invert(self) -> Self {
        self.conjugate().mul_base(self.norm().invert())
    }

    /// b^2 = (b0^2 + b1^2) + b1^2 u, so b1 = sqrt(a1) and b0 = sqrt(a0 + a1).
    #[inline(always)]
    pub fn sqrt(self) -> Self {
        let (a0, a1) = self.parts();
        Self::from_parts((a0 + a1).sqrt(), a1.sqrt())
    }

    /// Tr_122(a) = Tr_61(Tr_{122/61}(a)) = Tr_61(a1), as a + a^(2^61) = a1.
    /// So Tr(u) = Tr_61(1) = 1.
    #[inline(always)]
    pub fn trace(self) -> u32 {
        self.parts().1.trace()
    }

    /// A z with z^2 + z = d + u Tr(d), for any d (Pornin's QSolve, ePrint
    /// 2023/1688 §2.1): so z^2 + z = d exactly when Tr(d) = 0, which is
    /// when the equation is solvable at all; the other root is z + 1.
    ///
    /// With z = z0 + z1 u, z^2 + z = (z0^2 + z0 + z1^2) + (z1^2 + z1) u:
    /// - z1 = H61(d1) gives z1^2 + z1 = d1 + Tr61(d1), and Tr_122(d) = Tr61(d1);
    /// - z0 needs z0^2 + z0 = d0 + z1^2, solvable iff Tr61(d0) = Tr61(z1).
    ///   Both roots z1, z1 + 1 work above and their traces differ
    ///   (Tr61(1) = 1), so take the one with Tr61(z1) = Tr61(d0): this
    ///   second condition never fails;
    /// - then z0 = H61(d0 + z1^2). With h = H61(d1) and the chosen flip bit f,
    ///   z1^2 = h + d1 + Tr61(d1) + f, so no squaring is needed.
    ///
    /// 2 H61 and a few word operations; both output coefficients are canonical.
    #[inline(always)]
    pub fn qsolve(self) -> Self {
        let (d0, d1) = self.parts();
        let h = d1.halftrace();
        let flip = Gf61::w64((h.trace() ^ d0.trace()) as u64);
        let z1 = h + flip;
        let z1sq = h + d1 + Gf61::w64(d1.trace() as u64) + flip;
        Self::from_parts((d0 + z1sq).halftrace(), z1)
    }

    /// One root z of z^2 + z = d, if Tr(d) = 0; the other is z + 1.
    pub fn solve_quadratic(self) -> Option<Self> {
        (self.trace() == 0).then(|| self.qsolve())
    }
}

impl fmt::Debug for Gf {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let v = to_u128(*self);
        write!(f, "Gf({:#018x} + {:#018x} u)", v as u64, (v >> 64) as u64)
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
    use super::gf2_61::tests::{mul_ref as mul61_ref, reduce_ref, trace_ref};
    use super::*;
    use proptest::prelude::*;

    fn split(v: u128) -> (u64, u64) {
        (reduce_ref(v as u64), reduce_ref((v >> 64) as u64))
    }

    fn join(a0: u64, a1: u64) -> u128 {
        a0 as u128 | (a1 as u128) << 64
    }

    /// Schoolbook in the tower over the bit-serial GF(2^61) reference:
    /// (a0 + a1 u)(b0 + b1 u) = a0 b0 + a1 b1 + (a0 b1 + a1 b0 + a1 b1) u.
    pub fn mul_ref(a: u128, b: u128) -> u128 {
        let ((a0, a1), (b0, b1)) = (split(a), split(b));
        let m = mul61_ref;
        let a1b1 = m(a1, b1);
        join(m(a0, b0) ^ a1b1, m(a0, b1) ^ m(a1, b0) ^ a1b1)
    }

    /// a^(2^e)
    fn frob(a: u128, e: u32) -> u128 {
        (0..e).fold(to_u128(raw_of(a)), |x, _| mul_ref(x, x))
    }

    fn raw_of(v: u128) -> Gf {
        Gf::w64le(v as u64, (v >> 64) as u64)
    }

    pub fn fe() -> impl Strategy<Value = u128> {
        any::<u128>().prop_map(|v| v & MASK122)
    }

    /// Any pair of 64-bit representatives, as products leave them.
    pub fn raw() -> impl Strategy<Value = Gf> {
        any::<u128>().prop_map(raw_of)
    }

    proptest! {
        #[test]
        fn roundtrip(a in fe()) {
            prop_assert_eq!(to_u128(from_u128(a)), a);
            prop_assert_eq!(from_bytes(to_bytes(from_u128(a))).map(to_u128), Some(a));
        }

        #[test]
        fn from_bytes_rejects_non_canonical(v in any::<u128>()) {
            prop_assert_eq!(from_bytes(v.to_le_bytes()).is_some(), v & !MASK122 == 0);
        }

        #[test]
        fn normalization_reduces_each_half(v in any::<u128>()) {
            let (a0, a1) = split(v);
            prop_assert_eq!(to_u128(raw_of(v)), join(a0, a1));
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
        fn soft_backend_matches_native(a in any::<u128>(), b in any::<u128>()) {
            let s = |v: u128| soft::Gf::w64le(v as u64, (v >> 64) as u64);
            let n = |x: soft::Gf| to_u128(Gf::w64le(x.limbs()[0], x.limbs()[1]));
            let (x, y) = (raw_of(a), raw_of(b));
            prop_assert_eq!(n(soft::mul(s(a), s(b))), to_u128(x * y));
            prop_assert_eq!(n(soft::square(s(a))), to_u128(x.square()));
            prop_assert_eq!(n(soft::add(s(a), s(b))), to_u128(x + y));
            prop_assert_eq!(n(soft::mul(s(a), s(b))), mul_ref(to_u128(x), to_u128(y)));
            let (a0, b0) = (a as u64, b as u64);
            prop_assert_eq!(reduce_ref(soft::mul61(a0, b0)), gf2_61::to_u64(Gf61::w64(a0) * Gf61::w64(b0)));
            prop_assert_eq!(reduce_ref(soft::square61(a0)), gf2_61::to_u64(Gf61::w64(a0).square()));
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
        fn norm_is_a_times_its_conjugate(a in raw()) {
            let c = frob(to_u128(a), 61);
            prop_assert_eq!(to_u128(a.conjugate()), c);
            prop_assert_eq!(mul_ref(to_u128(a), c), gf2_61::to_u64(a.norm()) as u128);
        }

        #[test]
        fn sqrt(a in raw()) {
            let s = a.sqrt();
            prop_assert!(eq(s.square(), a));
            prop_assert_eq!(to_u128(s), frob(to_u128(a), 121));
        }

        #[test]
        fn trace_matches_definition(a in raw()) {
            // Tr(v) = sum of v^(2^i), i < 122, which lies in F_2
            let v = to_u128(a);
            let t = (0..122).fold((0, v), |(t, x), _| (t ^ x, mul_ref(x, x))).0;
            prop_assert!(t <= 1);
            prop_assert_eq!(a.trace() as u128, t);
            prop_assert_eq!(a.trace() as u64, trace_ref(a.limbs()[1]));
        }

        #[test]
        fn qsolve_solves_the_quadratic(d in raw()) {
            let z = d.qsolve();
            let tr = d.trace();
            let want = if tr == 1 { d + Gf::U } else { d };
            prop_assert!(eq(z.square() + z, want));
            let zv = to_u128(z);
            prop_assert_eq!(mul_ref(zv, zv) ^ zv, to_u128(want));
            match d.solve_quadratic() {
                Some(z) => prop_assert!(tr == 0 && eq(z.square() + z, d)),
                None => prop_assert_eq!(tr, 1),
            }
        }

        #[test]
        fn mul_u_and_mul_base(a in raw(), c in any::<u64>()) {
            prop_assert!(eq(a.mul_u(), a * Gf::U));
            let c = Gf61::w64(c);
            prop_assert!(eq(a.mul_base(c), a * Gf::from_base(c)));
        }

        #[test]
        fn field_axioms(a in raw(), b in raw(), c in raw()) {
            prop_assert!(eq(a * (b + c), a * b + a * c));
            prop_assert!(eq((a * b) * c, a * (b * c)));
            prop_assert!(eq(a * b, b * a));
            prop_assert!(eq(a * Gf::ONE, a) && is_zero(a * Gf::ZERO) && is_zero(a + a));
            prop_assert!(eq((a + b).square(), a.square() + b.square()));
            prop_assert!(eq(Gf::U.square() + Gf::U, Gf::ONE));
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
    fn known_answers() {
        use super::kats::{MUL, UNARY};
        for &(a, b, c) in MUL {
            assert_eq!(to_u128(from_u128(a) * from_u128(b)), c, "{a:#x} * {b:#x}");
        }
        for &(a, sq, inv, sqrt, tr) in UNARY {
            let x = from_u128(a);
            assert_eq!(to_u128(x.square()), sq, "{a:#x}^2");
            assert_eq!(to_u128(x.invert()), inv, "1/{a:#x}");
            assert_eq!(to_u128(x.sqrt()), sqrt, "sqrt {a:#x}");
            assert_eq!(x.trace(), tr, "Tr {a:#x}");
            let z = x.qsolve();
            let rhs = if tr == 1 { x + Gf::U } else { x };
            assert!(eq(z.square() + z, rhs), "qsolve {a:#x}");
        }
    }
}
