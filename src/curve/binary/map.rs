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

use super::{Constant, Curve, Model, Point};
use crate::curve::h2c::Map;
use crate::field::batch::{Invert, invert as batch_invert};
use crate::field::{Binary, Field};
use crate::hash::{Salted, halves};

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

impl<M: Pornin> Map<Point<M>> for Curve<M> {
    /// The map of c, whose bit `Model::SIGN` selects between +P and -P.
    fn map(&self, c: u128) -> Point<M> {
        let st = Self::map_prepare(c);
        self.map_finish(&st, st.den.inv())
    }
}

impl<M: Pornin> Curve<M> {
    /// Pornin's map on many inputs with one shared inversion.
    pub fn map_to_curve_batch(&self, cs: &[u128]) -> Vec<Point<M>> {
        self.map_batch(cs, |st, inv| self.map_finish(st, inv))
    }

    /// `finish` of each input's state and the inverse of its denominator,
    /// the inversions shared.
    pub(super) fn map_batch<P>(
        &self,
        cs: &[u128],
        finish: impl Fn(&MapState<M::F>, M::F) -> P,
    ) -> Vec<P> {
        let ms: Vec<MapState<M::F>> = cs.iter().map(|&c| Self::map_prepare(c)).collect();
        let mut inv: Vec<M::F> = ms.iter().map(|m| m.den).collect();
        batch_invert(&mut inv);
        ms.iter().zip(inv).map(|(m, i)| finish(m, i)).collect()
    }

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

    pub fn map_finish(&self, st: &MapState<M::F>, inv_den: M::F) -> Point<M> {
        let MapRoot { x, w, .. } = self.map_root(st, inv_den);
        let s = x * w.square();
        Point {
            x,
            s: self.beta.mul(s),
            z: self.beta.gf(),
            t: self.beta.mul(x),
        }
    }

    pub fn hash_to_curve_map1_batch(&self, h: &Salted, msgs: &[&[u8]]) -> Vec<Point<M>> {
        let cs: Vec<u128> = msgs.iter().map(|m| halves(&h.digest(m, 0))[0]).collect();
        self.map_to_curve_batch(&cs)
    }

    /// crrl-style map(h0) + map(h1), conjectured indifferentiable from a
    /// random oracle. The sum of two encodings is indifferentiable when the
    /// map is well distributed (Farashahi, Fouque, Shparlinski, Tibouchi and
    /// Voloch, 2013); that property has not been established for these
    /// instantiations of the map.
    pub fn hash_to_curve_map2(&self, h: &Salted, msg: &[u8]) -> Point<M> {
        let [c0, c1] = halves(&h.digest(msg, 0));
        self.add(&self.map(c0), &self.map(c1))
    }
}

/// The map's affine outputs against (input, x, y) known answers.
#[cfg(test)]
pub(crate) fn check_map_vectors<M: Pornin>(c: &Curve<M>, vectors: &[(u128, u128, u128)]) {
    for &(v, x, y) in vectors {
        let a = c.to_affine(&c.map(v));
        assert_eq!((a.x.value(), a.y.value()), (x, y), "input {v:#x}");
    }
}

/// The map's properties, for each family's `map` tests module, over
/// curves from `$curve`.
#[cfg(test)]
macro_rules! map_suite {
    ($m:ty, $curve:expr) => {
        use crate::curve::binary::{Constant, Model};
        use crate::curve::h2c::{Map, map1};
        use crate::field::batch::Invert;
        use crate::field::Field;
        use crate::hash::Salted;
        use proptest::prelude::*;

        type M = $m;
        type Curve = crate::curve::binary::Curve<M>;

        proptest! {
            #[test]
            fn map_lands_in_the_group_and_sign_negates(c in $curve, v in any::<u128>()) {
                let p = c.map(v);
                let a = c.to_affine(&p);
                prop_assert!(!a.is_identity() && c.is_on_curve(&a) && a.x.trace() == M::A.trace());
                prop_assert_eq!(c.to_affine(&c.map(v ^ 1 << M::SIGN)), a.neg());
            }

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

            #[test]
            fn map_output_decodes_to_itself(c in $curve, v in any::<u128>()) {
                let a = c.to_affine(&c.map(v));
                prop_assert_eq!(c.decode(a.encode()), Some(a));
            }

            #[test]
            fn map_hashes_land_in_group(c in $curve, salt in any::<[u8; 32]>(), msg in any::<[u8; 8]>()) {
                let h = Salted::new(b"test", &salt);
                for q in [map1(&c, &h, &msg), c.hash_to_curve_map2(&h, &msg)] {
                    let q = c.to_affine(&q);
                    prop_assert!(c.is_on_curve(&q) && (q.is_identity() || q.x.trace() == M::A.trace()));
                }
            }

            #[test]
            fn batched_map_matches_single(c in $curve, salt in any::<[u8; 32]>(), msgs in prop::collection::vec(prop::collection::vec(any::<u8>(), 0..8), 0..20)) {
                let h = Salted::new(b"test", &salt);
                let refs: Vec<&[u8]> = msgs.iter().map(|m| m.as_slice()).collect();
                let batch = c.hash_to_curve_map1_batch(&h, &refs);
                prop_assert_eq!(batch.len(), refs.len());
                for (m, p) in refs.iter().zip(batch) {
                    prop_assert_eq!(c.to_affine(&p), c.to_affine(&map1(&c, &h, m)));
                }
            }
        }
    };
}
#[cfg(test)]
pub(crate) use map_suite;
