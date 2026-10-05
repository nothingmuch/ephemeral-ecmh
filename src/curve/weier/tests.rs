//! Shared property tests and curve/point strategies for `weier` instances.
//! `suite!(fp127, F, fp)` takes the field module, type and element strategy.

macro_rules! suite {
    ($m:ident, $f:ident, $elem:ident) => {
        use crate::curve::encoding::{Signed, mask};
        use crate::field::$m::tests::$elem as elem;
        use crate::field::Packed;
        use crate::hash::Salted;
        use proptest::prelude::*;

        type F = crate::field::$m::$f;
        type Curve = crate::curve::weier::Curve<F>;
        type OddCurve = crate::curve::weier::OddCurve<F>;
        type Affine = crate::curve::weier::Affine<F>;
        type Point = crate::curve::weier::Point<F>;

        /// Curves of odd order, where every formula is complete.
        pub fn odd_curve() -> impl Strategy<Value = OddCurve> {
            elem().prop_filter_map("even order", |b| Curve::new(b).and_then(OddCurve::new))
        }

        /// Curves of any order: the chord law, codec, hash and sieve hold on
        /// every one, and certificate verification computes on all of them.
        pub fn curve() -> impl Strategy<Value = Curve> {
            elem().prop_filter_map("singular", Curve::new)
        }

        pub fn point(c: Curve) -> impl Strategy<Value = Affine> {
            any::<u128>().prop_filter_map("not on curve", move |v| {
                c.decode(F::to_bytes(v & mask::<F>())).filter(|p| !p.is_identity())
            })
        }

        pub fn curve_and_points(n: usize) -> impl Strategy<Value = (Curve, Vec<Affine>)> {
            curve().prop_flat_map(move |c| (Just(c), prop::collection::vec(point(c), n)))
        }

        /// Points on the raw form of an odd-order curve, for RCB's formulas.
        pub fn odd_curve_and_points(n: usize) -> impl Strategy<Value = (Curve, Vec<Affine>)> {
            odd_curve().prop_flat_map(move |c| (Just(c.curve()), prop::collection::vec(point(c.curve()), n)))
        }

        proptest! {
            #[test]
            fn odd_curve_refuses_a_root(x0 in elem()) {
                // (x0, 0) is a point of order 2
                let b = (F::ONE + F::ONE + F::ONE) * x0 - x0 * x0.square();
                let c = Curve { b };
                prop_assert!(OddCurve::new(c).is_none());
            }

            #[test]
            fn x_on_curve_matches_decode(c in curve(), x in elem()) {
                prop_assert_eq!(c.x_on_curve(x), c.decode(F::to_bytes(x.pack())).is_some());
            }

            #[test]
            fn encode_roundtrip((c, ps) in curve_and_points(1)) {
                let p = ps[0];
                prop_assert!(c.is_on_curve(&p));
                prop_assert_eq!(c.decode(p.encode()), Some(p));
                prop_assert_eq!(c.decode(Affine::IDENTITY.encode()), Some(Affine::IDENTITY));
            }

            #[test]
            fn decode_accepts_only_canonical((c, v) in (curve(), any::<<F as Signed>::Bytes>())) {
                if let Some(p) = c.decode(v) {
                    prop_assert!(c.is_on_curve(&p));
                    prop_assert_eq!(p.encode(), v);
                }
            }

            #[test]
            fn group_laws((c, ps) in curve_and_points(3)) {
                let (p, q, s) = (ps[0], ps[1], ps[2]);
                let o = Affine::IDENTITY;
                prop_assert!(c.is_on_curve(&c.add(&p, &q)));
                prop_assert_eq!(c.add(&p, &q), c.add(&q, &p));
                prop_assert_eq!(c.add(&c.add(&p, &q), &s), c.add(&p, &c.add(&q, &s)));
                prop_assert_eq!(c.add(&p, &o), p);
                prop_assert_eq!(c.add(&p, &p.neg()), o);
                prop_assert_eq!(c.add(&c.add(&p, &p), &p.neg()), p);
            }

            #[test]
            fn batch_matches_single((c, ps) in curve_and_points(24), mask in any::<u32>()) {
                let a: Vec<Affine> = ps[..12].to_vec();
                let b: Vec<Affine> = ps[12..].iter().enumerate().map(|(i, q)| match mask >> (2 * i) & 3 {
                    0 => *q,
                    1 => a[i],
                    2 => a[i].neg(),
                    _ => Affine::IDENTITY,
                }).collect();
                let mut out = vec![Affine::IDENTITY; 12];
                c.add_batch(&a, &b, &mut out);
                for i in 0..12 {
                    prop_assert_eq!(out[i], c.add(&a[i], &b[i]));
                }
            }

            #[test]
            fn sum_batch_matches_fold((c, ps) in curve_and_points(17), k in 0usize..17) {
                let seq = ps[..k].iter().fold(Affine::IDENTITY, |acc, p| c.add(&acc, p));
                prop_assert_eq!(c.sum_batch(&ps[..k]), seq);
            }

            #[test]
            fn projective_matches_affine((c, ps) in odd_curve_and_points(3)) {
                let (p, q, s) = (ps[0], ps[1], ps[2]);
                let (ep, eq_) = (c.from_affine(&p), c.from_affine(&q));
                prop_assert_eq!(c.to_affine(&c.from_affine(&p)), p);
                prop_assert_eq!(c.to_affine(&ep.neg()), p.neg());
                prop_assert_eq!(c.to_affine(&c.add_ext(&ep, &eq_)), c.add(&p, &q));
                prop_assert_eq!(c.to_affine(&c.add_ext(&ep, &ep)), c.add(&p, &p));
                prop_assert_eq!(c.to_affine(&c.add_ext(&ep, &ep.neg())), Affine::IDENTITY);
                prop_assert_eq!(c.to_affine(&c.add_ext(&Point::IDENTITY, &ep)), p);
                let pqs = c.add_affine(&c.add_ext(&ep, &eq_), &s);
                prop_assert_eq!(c.to_affine(&pqs), c.add(&c.add(&p, &q), &s));
                prop_assert!(c.add_ext(&ep, &eq_).equals(&c.add_affine(&eq_, &p)));
                prop_assert!(!ep.equals(&eq_) && !ep.equals(&ep.neg()) && !ep.equals(&Point::IDENTITY));
                prop_assert!(c.add_ext(&ep, &ep.neg()).is_identity() && !ep.is_identity());
                let degenerate = Point { x: F::ZERO, y: F::ZERO, z: F::ZERO };
                prop_assert!(!degenerate.is_identity());
            }

            #[test]
            fn equals_addend_matches_default((c, ps) in odd_curve_and_points(2), l in elem().prop_filter("zero scale", |l| !l.is_zero())) {
                use crate::group::{Accumulate, Group, Negate};
                let g = OddCurve::new(c).unwrap();
                let (p, q) = (ps[0], ps[1]);
                let pq = c.add(&p, &q);
                let pts = [p, q, pq, p.neg(), Affine::IDENTITY];
                let scale = |e: Point| Point { x: e.x * l, y: e.y * l, z: e.z * l };
                let (ep, eq_) = (c.from_affine(&p), c.from_affine(&q));
                // representatives by different routes and scales, O among them
                let accs = [
                    (ep, p),
                    (scale(ep), p),
                    (c.add_ext(&ep, &eq_), pq),
                    (scale(c.add_affine(&eq_, &p)), pq),
                    (scale(ep.neg()), p.neg()),
                    (Point::IDENTITY, Affine::IDENTITY),
                    (scale(Point::IDENTITY), Affine::IDENTITY),
                    (scale(c.add_affine(&ep, &p.neg())), Affine::IDENTITY),
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
            fn mixed_matches_projective((c, ps) in odd_curve_and_points(3)) {
                // the addend at Z = 1 through Alg. 4, against Alg. 5, including O
                // on either side, doubling and cancellation, with Z1 != 1
                let (p, q, s) = (ps[0], ps[1], ps[2]);
                let acc = c.add_ext(&c.from_affine(&p), &c.from_affine(&q));
                let pq = c.add(&p, &q);
                for (a, b) in [(acc, s), (acc, pq), (acc, pq.neg()), (acc, Affine::IDENTITY), (Point::IDENTITY, s), (Point::IDENTITY, Affine::IDENTITY)] {
                    let want = c.add_ext(&a, &c.from_affine(&b));
                    let got = c.add_affine(&a, &b);
                    prop_assert!(got.equals(&want) && got.is_identity() == want.is_identity());
                    prop_assert_eq!(c.to_affine(&got), c.to_affine(&want));
                }
            }

            #[test]
            fn normalize_batch_matches((c, ps) in odd_curve_and_points(6)) {
                let mut ext: Vec<Point> = ps.windows(2).map(|w| c.add_affine(&c.from_affine(&w[0]), &w[1])).collect();
                ext.push(Point::IDENTITY);
                let want: Vec<Affine> = ext.iter().map(|p| c.to_affine(p)).collect();
                prop_assert_eq!(c.normalize_batch(&ext), want);
            }

            #[test]
            fn scalar_mul_is_repeated_addition((c, ps) in odd_curve_and_points(1), k in 0u128..40) {
                let p = ps[0];
                let want = (0..k).fold(Affine::IDENTITY, |acc, _| c.add(&acc, &p));
                prop_assert_eq!(c.to_affine(&c.mul(&c.from_affine(&p), k)), want);
            }

            #[test]
            fn hash_lands_on_curve(c in curve(), salt in any::<[u8; 32]>(), msg in any::<[u8; 8]>()) {
                let p = c.hash_to_curve(&Salted::new(b"test", &salt), &msg);
                prop_assert!(!p.is_identity() && c.is_on_curve(&p));
            }
        }
    };
}

pub(crate) use suite;
