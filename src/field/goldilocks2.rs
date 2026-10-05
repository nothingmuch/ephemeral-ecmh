//! F_{p^2}, p = 2^64 - 2^32 + 1 (Goldilocks), i^2 = 7, represented in
//! two words using Plonky3's base field and binomial extension.
//!
//! - Base field: `Fp`, `p3_goldilocks::Goldilocks`. Reduction is sparse,
//!   2^64 = 2^32 - 1, and p - 1 = 2^32 (2^32 - 1), which suits FFTs but
//!   not square roots: Tonelli-Shanks over a 2-Sylow subgroup of order
//!   2^32.
//! - Extension: `Fq`, Plonky3's binomial extension F_p\[i\]/(i^2 - 7). 7
//!   generates F_p^*, so it is a non-square and i^2 - 7 irreducible.
//! - Inversion: Plonky3's, through the norm a^2 - 7b^2.
//! - Square root: the complex method, as in `fp64x2`, but from Plonky3's
//!   Tonelli-Shanks and with an inversion for 1/x, which also serves the
//!   1/N(d) of a ratio n/d's root, through n conj(d) / N(d).
//!
//! Curves: `curve::twisted_goldilocks2`.

use p3_field::extension::BinomialExtensionField;
use p3_field::{BasedVectorSpace, Field, PrimeCharacteristicRing};

pub use p3_goldilocks::Goldilocks as Fp;
pub type Fq = BinomialExtensionField<Fp, 2>;

pub const P: u64 = 0xffff_ffff_0000_0001;

const SEVEN: Fp = Fp::new(7);
/// 1/7, avoiding a second inversion in the pure-base square-root branch.
const SEVEN_INV: Fp = Fp::new(0x2492_4924_6db6_db6e);

#[inline]
pub const fn new(a: Fp, b: Fp) -> Fq {
    Fq::new([a, b])
}

/// (a, b) of a + b i.
#[inline]
pub fn parts(x: Fq) -> (Fp, Fp) {
    let c = x.as_basis_coefficients_slice();
    (c[0], c[1])
}

/// (a + b i)(a - b i) = a^2 - 7b^2.
#[inline]
pub fn norm(x: Fq) -> Fp {
    let (a, b) = parts(x);
    a.square() - SEVEN * b.square()
}

/// 1/x; maps 0 to 0.
pub fn invert(x: Fq) -> Fq {
    x.try_inverse().unwrap_or(Fq::ZERO)
}

/// 1/x in F_p; maps 0 to 0.
pub fn invert_base(x: Fp) -> Fp {
    x.try_inverse().unwrap_or(Fp::ZERO)
}

/// Some square root of x in F_p if it is a square.
pub fn sqrt_base(x: Fp) -> Option<Fp> {
    x.try_sqrt()
}

/// Some square root if `x` is a square, equivalently if its norm is a
/// square in F_p.
pub fn sqrt(x: Fq) -> Option<Fq> {
    sqrt_over(x, Fp::ONE)
}

/// Some square root of n/d, if d != 0 and n/d is a square: of
/// n conj(d) / N(d), with the one inversion `sqrt` has.
pub fn sqrt_ratio(n: Fq, d: Fq) -> Option<Fq> {
    let nd = norm(d);
    if nd.is_zero() {
        return None;
    }
    let (a, b) = parts(d);
    sqrt_over(n * new(a, -b), nd)
}

/// Some square root of w/n, n != 0 in F_p, by the complex method.
///
/// (x + y i)^2 = a + b i is x^2 + 7y^2 = a, 2xy = b, and then the norm
/// m = x^2 - 7y^2 is a square root of a^2 - 7b^2, so x^2 = (a + m)/2. With
/// b != 0 the two signs of m give t = (a + m)/2 and (a - m)/2, whose
/// product is 7b^2/4, a non-square: exactly one is a square, x^2. Then
/// y = b/(2x). For a + b i = w/n: t is t'/n for t' = (w.a ± m')/2 and m'
/// a root of N(w), x = r/n for r a root of t' n, and y = w.b/(2r), both
/// from one inversion of 2 r n.
fn sqrt_over(w: Fq, n: Fp) -> Option<Fq> {
    let (a, b) = parts(w);
    if b.is_zero() {
        // a/n, or a/(7n) as 7 is a non-square, is a square in F_p
        let n_inv = invert_base(n);
        return Some(match (a * n).try_sqrt() {
            Some(r) => new(r * n_inv, Fp::ZERO),
            None => new(Fp::ZERO, (a * n * SEVEN_INV).try_sqrt()? * n_inv),
        });
    }
    let m = norm(w).try_sqrt()?;
    let r = ((a + m).halve() * n)
        .try_sqrt()
        .or_else(|| ((a - m).halve() * n).try_sqrt())?;
    let k = (r.double() * n).inverse();
    Some(new(r.double() * r * k, b * n * k))
}

impl crate::field::batch::Invert for Fq {
    const ONE: Self = <Fq as PrimeCharacteristicRing>::ONE;
    fn inv(self) -> Self {
        invert(self)
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use p3_field::PrimeField64;
    use proptest::prelude::*;

    pub fn fp() -> impl Strategy<Value = Fp> {
        prop_oneof![
            (0..P).prop_map(Fp::new),
            prop::sample::select(vec![0, 1, 2, 7, 1 << 32, (1 << 32) - 1, P - 1, P - 2])
                .prop_map(Fp::new),
        ]
    }

    pub fn fq() -> impl Strategy<Value = Fq> {
        (fp(), fp()).prop_map(|(a, b)| new(a, b))
    }

    #[test]
    fn modulus() {
        assert_eq!(P as u128, (1 << 64) - (1 << 32) + 1);
        assert_eq!(Fp::new(P - 1) + Fp::ONE, Fp::ZERO);
        // 7 is a non-square: 7^((p-1)/2) = -1
        assert_eq!(SEVEN.exp_u64((P - 1) / 2), -Fp::ONE);
        assert!(sqrt_base(SEVEN).is_none());
        assert_eq!(SEVEN * SEVEN_INV, Fp::ONE);
    }

    proptest! {
        #[test]
        fn base_ops_match_reference(a in fp(), b in fp()) {
            let (x, y) = (a.as_canonical_u64() as u128, b.as_canonical_u64() as u128);
            prop_assert_eq!((a * b).as_canonical_u64() as u128, x * y % P as u128);
            prop_assert_eq!((a + b).as_canonical_u64() as u128, (x + y) % P as u128);
            prop_assert_eq!(a * invert_base(a), if a.is_zero() { Fp::ZERO } else { Fp::ONE });
            match sqrt_base(a) {
                Some(s) => prop_assert_eq!(s.square(), a),
                None => prop_assert!(sqrt_base(a * SEVEN).is_some()),
            }
        }

        #[test]
        fn mul_matches_reference(x in fq(), y in fq()) {
            let ((a, b), (c, d)) = (parts(x), parts(y));
            prop_assert_eq!(x * y, new(a * c + SEVEN * b * d, a * d + b * c));
            prop_assert_eq!(x.square(), x * x);
            prop_assert_eq!(norm(x), parts(x * new(a, -b)).0);
        }

        #[test]
        fn field_laws(x in fq(), y in fq(), z in fq()) {
            prop_assert_eq!(x * (y + z), x * y + x * z);
            prop_assert_eq!(x - y + y, x);
            prop_assert_eq!(x * invert(x), if x.is_zero() { Fq::ZERO } else { Fq::ONE });
        }

        #[test]
        fn sqrt_exactly_on_squares(x in fq()) {
            let s = sqrt(x.square()).unwrap();
            prop_assert!(s == x || s == -x);
            match sqrt(x) {
                Some(s) => prop_assert_eq!(s.square(), x),
                None => prop_assert!(sqrt_base(norm(x)).is_none()),
            }
        }

        /// `sqrt_ratio` where n conj(d) is in F_p, here a N(d): a square,
        /// as every element of F_p is in F_{p^2}.
        #[test]
        fn sqrt_ratio_of_base_elements(a in fp(), d in fq()) {
            let n = new(a, Fp::ZERO) * d;
            match sqrt_ratio(n, d) {
                Some(s) => prop_assert_eq!(s.square() * d, n),
                None => prop_assert!(d.is_zero()),
            }
        }
    }

    #[test]
    fn i_squared_is_seven() {
        let i = new(Fp::ZERO, Fp::ONE);
        let seven = new(SEVEN, Fp::ZERO);
        assert_eq!(i * i, seven);
        assert_eq!(sqrt(seven).map(|s| s.square()), Some(seven));
    }
}
