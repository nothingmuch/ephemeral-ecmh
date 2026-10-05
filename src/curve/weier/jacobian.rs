//! Jacobian accumulators (X : Y : Z), x = X/Z^2, y = Y/Z^3. The addend is
//! the `Affine` point itself, by the mixed addition madd-2007-bl,
//! 7M + 4S; doubling, dbl-2001-b for a = -3, is 3M + 5S. Neither formula
//! involves b.
//!
//! Against RCB's 11M + 2 m_b this trades completeness for branches: O as
//! either input and x_P = x_Q (P = Q or P = -Q) are tested for explicitly.
//! Branches are permitted because the inputs are public. A RIBLT's
//! additions reach these cases when a cell's sum is O, equals the addend,
//! or is its negative, as when an item is subtracted from a cell whose sum
//! is that item's point.

use super::Affine;
use crate::field::OddField;
use crate::field::batch::invert as batch_invert;
use crate::group::{Accumulate, Decode, Encode, Group, HashToCurve, Negate, SumBatch};
use crate::hash::Salted;

/// A `weier::OddCurve` under the Jacobian `Group` impl: the same points,
/// hashes, encodings and digests.
#[derive(Clone, Copy, Debug)]
pub struct Curve<F>(pub super::OddCurve<F>);

/// Jacobian (X : Y : Z), x = X/Z^2, y = Y/Z^3; Z = 0 is O.
#[derive(Clone, Copy, Debug)]
pub struct Point<F> {
    pub x: F,
    pub y: F,
    pub z: F,
}

impl<F: OddField> From<Affine<F>> for Point<F> {
    fn from(a: Affine<F>) -> Self {
        if a.is_identity() {
            return Self::IDENTITY;
        }
        Self {
            x: a.x,
            y: a.y,
            z: F::ONE,
        }
    }
}

impl<F: OddField> Point<F> {
    /// (1 : 1 : 0); `is_identity` takes any Z = 0 as O.
    pub const IDENTITY: Self = Self {
        x: F::ONE,
        y: F::ONE,
        z: F::ZERO,
    };

    pub fn is_identity(&self) -> bool {
        self.z.is_zero()
    }

    pub fn neg(&self) -> Self {
        Self {
            y: -self.y,
            ..*self
        }
    }

    /// dbl-2001-b, a = -3: 3M + 5S. Y = 0 would be 2-torsion, which
    /// odd-order curves lack; Z = 0 stays O.
    pub fn double(&self) -> Self {
        let delta = self.z.square();
        let gamma = self.y.square();
        let beta = self.x * gamma;
        let t = (self.x - delta) * (self.x + delta);
        let alpha = t + t + t;
        let beta4 = beta + beta;
        let beta4 = beta4 + beta4;
        let x = alpha.square() - (beta4 + beta4);
        let z = (self.y + self.z).square() - gamma - delta;
        let g2 = gamma.square();
        let g2 = g2 + g2;
        let g2 = g2 + g2;
        Self {
            x,
            y: alpha * (beta4 - x) - (g2 + g2),
            z,
        }
    }

    /// madd-2007-bl, 7M + 4S, with branches for its exceptions.
    pub fn add_affine(&self, q: &Affine<F>) -> Self {
        if q.is_identity() {
            return *self;
        }
        if self.is_identity() {
            return (*q).into();
        }
        let z1z1 = self.z.square();
        let u2 = q.x * z1z1;
        let s2 = q.y * (self.z * z1z1);
        let h = u2 - self.x;
        let r = s2 - self.y;
        // x_P = x_Q: P = Q, where the formulas collapse to (0 : 0 : 0), or
        // P = -Q
        if h.is_zero() {
            return if r.is_zero() {
                self.double()
            } else {
                Self::IDENTITY
            };
        }
        let r = r + r;
        let hh = h.square();
        let i = hh + hh;
        let i = i + i;
        let j = h * i;
        let v = self.x * i;
        let x = r.square() - j - (v + v);
        let yj = self.y * j;
        Self {
            x,
            y: r * (v - x) - (yj + yj),
            z: (self.z + h).square() - z1z1 - hh,
        }
    }

    /// X1 Z2^2 = X2 Z1^2 and Y1 Z2^3 = Y2 Z1^3: 6M + 2S, no inversion.
    pub fn equals(&self, o: &Self) -> bool {
        match (self.is_identity(), o.is_identity()) {
            (true, true) => true,
            (false, false) => {
                let (zz1, zz2) = (self.z.square(), o.z.square());
                self.x * zz2 == o.x * zz1 && self.y * (zz2 * o.z) == o.y * (zz1 * self.z)
            }
            _ => false,
        }
    }

    /// One inversion (of Z) and 3M + 1S.
    pub fn to_affine(&self) -> Affine<F> {
        self.normalize_with(self.z.inv())
    }

    fn normalize_with(&self, u: F) -> Affine<F> {
        if self.is_identity() {
            return Affine::IDENTITY;
        }
        let uu = u.square();
        Affine {
            x: self.x * uu,
            y: self.y * (uu * u),
        }
    }
}

/// `Point::to_affine` on many points with one shared inversion.
pub fn normalize_batch<F: OddField>(ps: &[Point<F>]) -> Vec<Affine<F>> {
    let mut u: Vec<F> = ps
        .iter()
        .map(|p| if p.is_identity() { F::ONE } else { p.z })
        .collect();
    batch_invert(&mut u);
    ps.iter().zip(u).map(|(p, u)| p.normalize_with(u)).collect()
}

impl<F: OddField> Group for Curve<F> {
    type Affine = Affine<F>;
    type Point = Point<F>;
    fn identity(&self) -> Point<F> {
        Point::IDENTITY
    }
    fn is_identity(&self, p: &Point<F>) -> bool {
        p.is_identity()
    }
    fn to_affine(&self, p: &Point<F>) -> Affine<F> {
        p.to_affine()
    }
    fn to_affine_batch(&self, ps: &[Point<F>]) -> Vec<Affine<F>> {
        normalize_batch(ps)
    }
}

impl<F: OddField> Accumulate for Curve<F> {
    type Addend = Affine<F>;
    fn prepare(&self, a: &Affine<F>) -> Affine<F> {
        *a
    }
    fn add(&self, p: &Point<F>, a: &Affine<F>) -> Point<F> {
        p.add_affine(a)
    }
    fn add_affine(&self, p: &Point<F>, a: &Affine<F>) -> Point<F> {
        p.add_affine(a)
    }
}

impl<F: OddField> Negate for Curve<F> {
    fn neg(&self, a: &Affine<F>) -> Affine<F> {
        a.neg()
    }
    fn neg_addend(&self, a: &Affine<F>) -> Affine<F> {
        a.neg()
    }
    /// X = x Z^2 and Y = y Z^3: 3M + 1S, and 1M + 1S when the x test fails.
    fn equals_addend(&self, p: &Point<F>, a: &Affine<F>) -> bool {
        if a.is_identity() || p.is_identity() {
            return a.is_identity() && p.is_identity();
        }
        let zz = p.z.square();
        p.x == a.x * zz && p.y == a.y * (zz * p.z)
    }
}

impl<F: OddField> SumBatch for Curve<F> {
    fn sum_batch(&self, a: &[Affine<F>]) -> Affine<F> {
        self.0.sum_batch(a)
    }
}

impl<F: OddField> HashToCurve for Curve<F>
where
    super::OddCurve<F>: HashToCurve<Affine = Affine<F>>,
{
    fn hash(&self, h: &Salted, msg: &[u8]) -> Affine<F> {
        self.0.hash(h, msg)
    }
}

impl<F: OddField> Encode for Curve<F>
where
    super::OddCurve<F>: Encode<Affine = Affine<F>>,
{
    type Encoding = <super::OddCurve<F> as Encode>::Encoding;
    fn encode(&self, a: &Affine<F>) -> Self::Encoding {
        self.0.encode(a)
    }
}

impl<F: OddField> Decode for Curve<F>
where
    super::OddCurve<F>: Decode<Affine = Affine<F>>,
{
    fn decode(&self, e: &Self::Encoding) -> Option<Affine<F>> {
        self.0.decode(e)
    }
}

#[cfg(test)]
mod tests {
    macro_rules! suite {
        ($family:ident, $m:ident, $f:ident, $elem:ident) => {
            mod $family {
                use super::super::normalize_batch;
                use crate::curve::$family::Affine;
                use crate::curve::$family::tests::{curve_and_points, odd_curve_and_points};
                use crate::field::$m::tests::$elem as elem;
                use proptest::prelude::*;

                type Point = super::super::Point<crate::field::$m::$f>;

        proptest! {
            #[test]
            fn from_affine_roundtrip((c, ps) in curve_and_points(1)) {
                let p = ps[0];
                prop_assert_eq!(Point::from(p).to_affine(), p);
                prop_assert_eq!(Point::from(p).neg().to_affine(), p.neg());
                prop_assert!(Point::from(Affine::IDENTITY).is_identity());
                prop_assert_eq!(Point::IDENTITY.to_affine(), Affine::IDENTITY);
                prop_assert!(Point::IDENTITY.double().is_identity());
                prop_assert_eq!(Point::from(p).double().to_affine(), c.add(&p, &p));
            }

            #[test]
            fn mixed_matches_affine((c, ps) in curve_and_points(3)) {
                let (p, q, s) = (ps[0], ps[1], ps[2]);
                let o = Affine::IDENTITY;
                let jp = Point::from(p);
                prop_assert_eq!(jp.add_affine(&q).to_affine(), c.add(&p, &q));
                prop_assert_eq!(Point::IDENTITY.add_affine(&q).to_affine(), q);
                prop_assert_eq!(jp.add_affine(&o).to_affine(), p);
                prop_assert!(Point::IDENTITY.add_affine(&o).is_identity());
                prop_assert_eq!(jp.add_affine(&p).to_affine(), c.add(&p, &p));
                prop_assert!(jp.add_affine(&p.neg()).is_identity());
                // the exceptional cases again with Z != 1
                let pq = jp.add_affine(&q);
                let spq = c.add(&p, &q);
                prop_assert_eq!(pq.add_affine(&spq).to_affine(), c.add(&spq, &spq));
                prop_assert!(pq.add_affine(&spq.neg()).is_identity());
                prop_assert_eq!(pq.add_affine(&s).to_affine(), c.add(&spq, &s));
                prop_assert_eq!(pq.double().to_affine(), c.add(&spq, &spq));
            }

            #[test]
            fn chains_match_affine((c, ps) in curve_and_points(8), ops in prop::collection::vec((0usize..8, 0u8..4), 0..48)) {
                // Ops 2 and 3 add or subtract the running sum, exercising
                // doubling and cancellation after arbitrary preceding additions.
                let (mut acc, mut want) = (Point::IDENTITY, Affine::IDENTITY);
                for (i, op) in ops {
                    let a = match op {
                        0 => ps[i],
                        1 => ps[i].neg(),
                        2 => want,
                        _ => want.neg(),
                    };
                    acc = acc.add_affine(&a);
                    want = c.add(&want, &a);
                    prop_assert_eq!(acc.to_affine(), want);
                    prop_assert_eq!(acc.is_identity(), want.is_identity());
                }
            }

            #[test]
            fn equality((_c, ps) in curve_and_points(2)) {
                let (p, q) = (ps[0], ps[1]);
                let (jp, jq) = (Point::from(p), Point::from(q));
                // same point, different representatives
                prop_assert!(jp.add_affine(&q).equals(&jq.add_affine(&p)));
                prop_assert!(!jp.equals(&jq) && !jp.equals(&jp.neg()) && !jp.equals(&Point::IDENTITY));
                prop_assert!(Point::IDENTITY.equals(&jp.add_affine(&p.neg())));
            }

            #[test]
            fn equals_addend_matches_default((c, ps) in odd_curve_and_points(2), l in elem().prop_filter("zero scale", |l| !l.is_zero())) {
                use crate::group::{Accumulate, Group, Negate};
                let g = super::super::Curve(crate::curve::weier::OddCurve::new(c).unwrap());
                let (p, q) = (ps[0], ps[1]);
                let pq = c.add(&p, &q);
                let pts = [p, q, pq, p.neg(), Affine::IDENTITY];
                let ll = l.square();
                let scale = |e: Point| Point { x: e.x * ll, y: e.y * ll * l, z: e.z * l };
                let (jp, jq) = (Point::from(p), Point::from(q));
                // representatives by different routes and scales, O among them
                let accs = [
                    (jp, p),
                    (scale(jp), p),
                    (jp.add_affine(&q), pq),
                    (scale(jq.add_affine(&p)), pq),
                    (scale(jp.neg()), p.neg()),
                    (Point::IDENTITY, Affine::IDENTITY),
                    (scale(Point::IDENTITY), Affine::IDENTITY),
                    (scale(jp.add_affine(&q).add_affine(&pq.neg())), Affine::IDENTITY),
                ];
                for (e, x) in accs {
                    for y in &pts {
                        let default = Group::is_identity(&g, &Accumulate::add(&g, &e, &Negate::neg_addend(&g, y)));
                        prop_assert_eq!(Negate::equals_addend(&g, &e, y), default);
                        prop_assert_eq!(default, x == *y);
                    }
                }
            }

            #[test]
            fn normalize_batch_matches((_c, ps) in curve_and_points(6)) {
                let mut pts: Vec<Point> = ps.windows(2).map(|w| Point::from(w[0]).add_affine(&w[1])).collect();
                pts.push(Point::IDENTITY);
                let want: Vec<Affine> = pts.iter().map(|p| p.to_affine()).collect();
                prop_assert_eq!(normalize_batch(&pts), want);
            }
        }            }
        };
    }

    suite!(weier127, fp127, Fp, fp);
    suite!(weier107, fp107, Fp, fp);
    suite!(weier61x2, fp61x2, Fq, fq);
}
