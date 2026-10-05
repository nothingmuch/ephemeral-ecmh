//! Simplified SWU for y² = x³ − 3x + b (RFC 9380, §6.6.2 and Appendix F.2).

use super::curve::cubic_irreducible;
use super::{Affine, Curve};
use crate::curve::{encoding, h2c::Map};
use crate::field::{OddField, Packed, Prime, fp127};
use core::fmt::Debug;

/// The ratio-root operation used by SSWU, including its nonsquare branch.
pub trait SswuField: OddField {
    type Ratio: Copy + Debug;

    /// Cache the field-specific constant for a nonsquare Z.
    fn prepare_ratio(z: Self) -> Self::Ratio;

    /// For d != 0, return (true, sqrt(n/d)), or (false, sqrt(Z*n/d)).
    fn ratio_root(n: Self, d: Self, cache: &Self::Ratio) -> (bool, Self);
}

macro_rules! three_mod_four {
    ($f:ty) => {
        impl SswuField for $f {
            type Ratio = Self;

            fn prepare_ratio(z: Self) -> Self {
                (-z).sqrt().expect("-Z is square when q = 3 mod 4")
            }

            #[inline]
            fn ratio_root(n: Self, d: Self, sqrt_neg_z: &Self) -> (bool, Self) {
                // r²d = n * (nd)^((q-1)/2), hence either n or -n.
                let r = n * (n * d).pow_p34();
                let square = r.square() * d == n;
                (square, if square { r } else { r * *sqrt_neg_z })
            }
        }
    };
}

three_mod_four!(fp127::Fp);

/// A deterministic map with per-curve setup amortized across inputs.
///
/// A single map does not have uniform output. Hash expansion and combining
/// independent mapped points belong to the hash-to-curve construction.
#[derive(Clone, Copy, Debug)]
pub struct Sswu<F: SswuField> {
    curve: Curve<F>,
    z: F,
    ratio: F::Ratio,
}

impl<F: SswuField> Sswu<F> {
    /// Validate the curve parameters and all four RFC 9380 conditions on Z.
    /// This validates the map, not the curve's group-order certificate.
    pub fn new(curve: Curve<F>, z: F) -> Option<Self> {
        let two = F::ONE + F::ONE;
        let a = -(two + F::ONE);
        if a.is_zero() || curve.b.is_zero() || curve.b == two || curve.b == -two {
            return None;
        }
        if z.is_square() || z == -F::ONE {
            return None;
        }
        let exceptional_x = curve.b * (z * a).inv();
        let exceptional_y2 = exceptional_x * (exceptional_x.square() + a) + curve.b;
        if !exceptional_y2.is_square() || !cubic_irreducible(a, curve.b - z) {
            return None;
        }
        Some(Self {
            curve,
            z,
            ratio: F::prepare_ratio(z),
        })
    }

    /// Map a field element; the output sign follows RFC 9380's sgn0(u).
    /// One ratio-root exponentiation and one inversion for the prime fields.
    pub fn map_to_curve(&self, u: F) -> Affine<F> {
        let zu2 = self.z * u.square();
        let denominator = zu2.square() + zu2;
        let numerator = self.curve.b * (denominator + F::ONE);
        let divisor = if denominator.is_zero() {
            self.z
        } else {
            -denominator
        };
        let divisor = -(divisor + divisor + divisor);
        // Evaluate g(numerator/divisor) as a ratio with denominator divisor³.
        let d2 = divisor.square();
        let d3 = d2 * divisor;
        let n = (numerator.square() - (d2 + d2 + d2)) * numerator + self.curve.b * d3;
        let (square, y) = F::ratio_root(n, d3, &self.ratio);
        let (x, y) = if square {
            (numerator, y)
        } else {
            (zu2 * numerator, zu2 * u * y)
        };
        Affine {
            x: x * divisor.inv(),
            y: if u.sgn0() == y.sgn0() { y } else { -y },
        }
    }
}

impl<F: SswuField> Sswu<F> {
    /// RFC 9380 Appendix H.2's search: Z = g + k, then -(g + k), for
    /// k = 0, 1, ..., count - 1, with g a generator of F over its prime
    /// field. Return the first valid setup, or None if this finite search
    /// fails.
    pub fn search_from(curve: Curve<F>, g: F, count: u32) -> Option<Self> {
        core::iter::successors(Some(g), |&z| Some(z + F::ONE))
            .take(count as usize)
            .flat_map(|z| [z, -z])
            .find_map(|z| Self::new(curve, z))
    }
}

impl<F: SswuField + Prime> Sswu<F> {
    /// Search Z in the order 1, -1, 2, -2, ..., max_abs_z, -max_abs_z.
    pub fn search(curve: Curve<F>, max_abs_z: u32) -> Option<Self> {
        Self::search_from(curve, F::ONE, max_abs_z)
    }
}

impl<F: SswuField + Packed> Map<Affine<F>> for Sswu<F> {
    /// Reduce the low BITS bits of a digest half, then apply the field map.
    fn map(&self, c: u128) -> Affine<F> {
        self.map_to_curve(encoding::x::<F>(c))
    }
}

#[cfg(test)]
mod tests;
