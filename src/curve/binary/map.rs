//! Pornin's deterministic map to E\[r\] (crrl `gls254::map_to_curve`),
//! for every family whose modulus supplies its constants ([`Pornin`]):
//! 1 I, shared by a batch, and 2 `solve`s, with no retry, against
//! try-and-increment's expected two tries.
//!
//! For a constant k with Tr(k) = 0 and an input c, take m1 = c,
//! m2 = c + k and m3 = c + c^2/k. As 1/c + 1/(c + k) = k/(c (c + k)) =
//! 1/m3, the e_i = b/m_i (b = sqrt(B)) satisfy e1 + e2 + e3 = 0 for every
//! c, k and B, so at least one e_i has trace 0, on random curves too. If
//! moreover
//!
//! - Tr(c) = Tr(a), and
//! - Tr(c^2/k) = Tr(c/sqrt(k)) = 0,
//!
//! then Tr(m_i) = Tr(a) = 1 for all three, so none is zero. For the first
//! m with Tr(b/m) = 0, d = sqrt(m) has Tr(d) = Tr(a) and Tr(b/d^2) = 0:
//! exactly the conditions under which Pornin's decoding of a w with
//! w^2 + w = d + a lands in E\[r\] (ePrint 2022/1325, §4.3; `wcodec`).
//! There x = d f for f^2 + f = b/d^2, the root of x^2 + d x = b with
//! Tr(x) = 0 (x = b/x̄, the extended coordinate), and s = x w^2. The
//! roots w and w + 1 are P and -P, and the input's bit `Model::SIGN`
//! picks one.
//!
//! The two trace conditions are affine conditions on c's bits, so a
//! family forces them by overwriting two bits that each enter one
//! condition. With k = z^2 the second condition is Tr(c/z) = 0, and the
//! division is a shift. This leaves m - 2 bits of c and the sign: the map
//! is not uniform on its own, and two maps added (`hash_to_curve_map2`)
//! are crrl's hash.
//!
//! The map's decoding is `wcodec`'s with x = b/x̄ kept instead of x̄, so
//! it lands on the other models' addends without a lift: the λ-affine
//! (x̄, λ) = (x + d, w^2 + 1 + a) and the unscaled (u, v) = (x/b, w^2 u),
//! each in a few field operations past the one shared inversion
//! ([`MapRoot`]), where try-and-increment pays an inversion to lift its
//! point.

use super::{Constant, Curve, Model};
use crate::field::Binary;
use crate::field::batch::Invert;

/// What the map takes from a family's modulus: k, c with both trace
/// conditions forced, c^2/k, and a bit that tells w from w + 1.
pub trait Pornin: Model {
    /// k, with Tr(k) = 0.
    fn k() -> Self::F;
    /// c: v's bits in `F::MASK`, two of them overwritten so that
    /// Tr(c) = Tr(a) and Tr(c^2/k) = 0.
    fn force(v: u128) -> Self::F;
    /// c^2/k.
    fn sq_over_k(c: Self::F) -> Self::F;
    /// A linear form that is 1 at 1, so it differs on w and w + 1.
    fn root_bit(w: Self::F) -> u32;
}

/// The input's part of the map, before the one inversion.
#[derive(Clone, Copy, Debug)]
pub struct MapState<F> {
    sign: u32,
    /// m1, m2, m3.
    pub(crate) m: [F; 3],
    mm: F,
    /// m1 m2 m3, the one inversion.
    pub den: F,
}

/// Where the map lands, before coordinates: the root x of X^2 + d X = b
/// with Tr(x) = 0 (x = b/x̄, the extended coordinate), the chosen w and
/// d = w^2 + w + a, `wcodec`'s d.
#[derive(Clone, Copy, Debug)]
pub struct MapRoot<F> {
    pub x: F,
    pub w: F,
    pub d: F,
}

impl<M: Pornin> Curve<M> {
    pub fn map_prepare(c: u128) -> MapState<M::F> {
        let sign = (c >> M::SIGN) as u32 & 1;
        let c = M::force(c);
        let (m1, m2, m3) = (c, c + M::k(), c + M::sq_over_k(c));
        let mm = m1 * m2;
        MapState {
            sign,
            m: [m1, m2, m3],
            mm,
            den: mm * m3,
        }
    }

    pub fn map_root(&self, st: &MapState<M::F>, inv_den: M::F) -> MapRoot<M::F> {
        let MapState {
            sign,
            m: [m1, m2, m3],
            mm,
            ..
        } = *st;
        let ii = self.b.mul(inv_den);
        let e3 = ii * mm;
        let jj = ii * m3;
        let (e1, e2) = (jj * m2, jj * m1);
        // At least one e_i has trace 0 since they sum to 0.
        let (m, e) = if e1.trace() == 0 {
            (m1, e1)
        } else if e2.trace() == 0 {
            (m2, e2)
        } else {
            (m3, e3)
        };
        // d = sqrt(m), w^2 + w = d + a (Tr(d) = Tr(a)), e = b/d^2
        let d = Binary::sqrt(m);
        let mut w = (d + M::A).solve();
        if M::root_bit(w) != sign {
            w += M::F::ONE;
        }
        let mut x = d * e.solve();
        if x.trace() == 1 {
            x += d;
        }
        MapRoot { x, w, d }
    }
}

/// The map's properties, for each family's `map` tests module, over
/// curves from `$curve`.
#[cfg(test)]
macro_rules! map_suite {
    ($m:ty, $curve:expr) => {
        use crate::curve::binary::{Constant, Model};
        use crate::field::Field;
        use crate::field::batch::Invert;
        use proptest::prelude::*;

        type M = $m;
        type Curve = crate::curve::binary::Curve<M>;

        proptest! {
            #[test]
            fn candidates_have_trace_a_and_cancel(c in $curve, v in any::<u128>()) {
                // Tr(m_i) = Tr(a) for all three, and b/m1 + b/m2 + b/m3 = 0
                let [m1, m2, m3] = Curve::map_prepare(v).m;
                for m in [m1, m2, m3] {
                    prop_assert_eq!(m.trace(), M::A.trace());
                }
                let e = |m: <M as Model>::F| c.b.mul(m.inv());
                prop_assert!((e(m1) + e(m2) + e(m3)).is_zero());
            }
        }
    };
}
#[cfg(test)]
pub(crate) use map_suite;
