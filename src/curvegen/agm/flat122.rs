//! Dense `binary122` curves: GF(2^122) from the tower `field::gf2_122`,
//! GF(2^61)\[u\]/(u^2 + u + 1), to the flat polynomial basis `Gf122`,
//! F_2\[z\]/(f) with f = z^122 + z^6 + z^2 + z + 1, which `Zq` lifts. A GLS
//! curve needs no such map (`Counter<Gf61>::order`); a dense-family B is
//! not restricted to GF(2^61), so it is counted in the full 122-bit field.
//!
//! G is a root of f in the tower, so z -> G is a ring map F_2\[z\]/(f) ->
//! GF(2^122), whose matrix has the columns G^i, i < 122. `ToFlat::new`
//! checks that f(G) = 0 and that the columns are independent; then G has
//! a minimal polynomial of degree 122, which is f, so f is irreducible and
//! the map is a field isomorphism. Its inverse, by elimination, is the
//! change of basis. Any root would serve, since an isomorphism preserves
//! #E; G is the least of the 122 as a `to_u128` integer (PARI's
//! `polrootsmod`, in a degree-122 field that the tower embeds in).
//!
//! f is the lowest irreducible pentanomial of degree 122; there is no
//! irreducible trinomial of that degree.

use super::{Counter, Gf122, checked, twist_order};
use crate::curve::binary122;
use crate::field::gf2_122::{Gf, eq, from_u128, to_u128};

/// The least root of f in the tower, as `to_u128`.
const G: u128 = 0x1c_2130_d4a4_0400_05c7_f442_4f3e_1f90;

/// An elimination basis of the map z -> G: entry k, if any, pairs a tower
/// vector whose highest set bit is k with the flat vector it is the image
/// of.
type Basis = [Option<(u128, u128)>; 128];

/// Clears the set bits of the tower vector t from the top while the basis
/// has their pivots, carrying the flat vector f along: t + G(f) is
/// invariant. Stops at the first pivot missing.
fn reduce(basis: &Basis, mut t: u128, mut f: u128) -> (u128, u128) {
    while t != 0 {
        let Some((bt, bf)) = basis[127 - t.leading_zeros() as usize] else {
            break;
        };
        (t, f) = (t ^ bt, f ^ bf);
    }
    (t, f)
}

/// The change of basis from the tower to `Gf122`, built once.
pub struct ToFlat(Basis);

impl Default for ToFlat {
    fn default() -> Self {
        Self::new()
    }
}

impl ToFlat {
    /// Panics if G is not a root of f, or its powers are dependent.
    pub fn new() -> Self {
        let g = from_u128(G);
        let mut basis = [None; 128];
        let mut power = Gf::ONE;
        for i in 0..122 {
            let (t, f) = reduce(&basis, to_u128(power), 1 << i);
            assert!(t != 0, "G^{i} depends on the lower powers of G");
            basis[127 - t.leading_zeros() as usize] = Some((t, f));
            power *= g;
        }
        let g2 = g * g;
        let g6 = g2 * g2 * g2;
        assert!(eq(power, g6 + g2 + g + Gf::ONE), "G is not a root of f");
        Self(basis)
    }

    /// The flat coordinates of v: bit i the coefficient of z^i.
    pub fn apply(&self, v: Gf) -> u128 {
        let (t, f) = reduce(&self.0, to_u128(v), 0);
        assert_eq!(t, 0, "the basis spans GF(2^122)");
        f
    }
}

/// Dense `binary122` curves, counted in the flat basis.
pub struct Counter122 {
    counter: Counter<Gf122>,
    to_flat: ToFlat,
}

impl Default for Counter122 {
    fn default() -> Self {
        Self::new()
    }
}

impl Counter122 {
    pub fn new() -> Self {
        Self {
            counter: Counter::new(),
            to_flat: ToFlat::new(),
        }
    }

    /// #E for E : y^2 + xy = x^3 + u x^2 + B over GF(2^122). m is even, so
    /// Tr(1) = 0, but Tr(u) = 1: E is the quadratic twist of the a2 = 0
    /// curve, and #E = q + 1 + t as for odd m and a = 1. Panics if the
    /// count fails `checked`'s sanity check.
    pub fn order(&self, c: &binary122::Dense) -> u128 {
        let t = self.counter.trace(self.to_flat.apply(c.big_b));
        checked(c, twist_order(122, t, 1 << 62))
    }
}

#[cfg(test)]
mod tests {
    use super::super::{Zq, counter61, counter122};
    use super::*;
    use crate::curve::binary122::tests::{dense, gls};
    use crate::field::gf2_122::MASK122;
    use proptest::prelude::*;

    fn flat(v: u128) -> u128 {
        counter122().to_flat.apply(from_u128(v))
    }

    #[test]
    fn to_flat_fixes_the_prime_field() {
        assert_eq!(flat(0), 0);
        assert_eq!(flat(1), 1);
        assert_eq!(flat(G), 2);
    }

    proptest! {
        /// Additive by construction; multiplicative, so a ring map.
        #[test]
        fn to_flat_is_a_ring_map(a in any::<u128>(), b in any::<u128>()) {
            let (a, b) = (a & MASK122, b & MASK122);
            let ab = to_u128(from_u128(a) * from_u128(b));
            let z = |v| Zq::<Gf122>::from_bits(flat(v));
            prop_assert_eq!(flat(a ^ b), flat(a) ^ flat(b));
            prop_assert_eq!(flat(ab), (z(a) * z(b)).bits());
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(16))]

        #[test]
        fn order_kills_points((c, ps) in dense::curve_and_points(2)) {
            let n = counter122().order(&c);
            prop_assert_eq!(n % 4, 2);
            for p in ps {
                let o = |k| c.mul(&c.from_affine(&p), k).is_identity();
                prop_assert!(o(n) && o(n / 2));
            }
        }

        /// A GLS curve is a dense curve too: the flat count of its B and
        /// the subfield identity agree.
        #[test]
        fn gls_by_either_counter(c in gls::curve()) {
            let dense = binary122::Dense::new(Gf::from_base(c.big_b));
            prop_assert_eq!(counter122().order(&dense), counter61().order(&c));
        }
    }
}
