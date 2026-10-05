//! Shared property tests and curve/point strategies for `twisted` instances.
//! `suite!(fp128, Fp, fp)` takes the field module, type and element strategy.

macro_rules! suite {
    ($m:ident, $f:ident, $elem:ident) => {
        use crate::curve::twisted::BYTES;
        use crate::field::$m::tests::$elem as elem;
        use crate::field::batch::Invert;
        use crate::field::{Field, OddField, Packed};
        use crate::hash::Salted;
        use proptest::prelude::*;

        type F = crate::field::$m::$f;
        type Curve = crate::curve::twisted::Curve<F>;
        type Affine = crate::curve::twisted::Affine<F>;
        type Point = crate::curve::twisted::Point<F>;

        pub fn curve() -> impl Strategy<Value = Curve> {
            elem().prop_filter_map("Edwards model incomplete", Curve::new)
        }

        pub fn point(c: Curve) -> impl Strategy<Value = Affine> {
            elem().prop_filter_map("not on curve", move |u| {
                c.decode_u128(u.pack()).filter(|p| !p.is_identity())
            })
        }

        pub fn curve_and_points(n: usize) -> impl Strategy<Value = (Curve, Vec<Affine>)> {
            curve().prop_flat_map(move |c| (Just(c), prop::collection::vec(point(c), n)))
        }

        /// The affine addition law, with its two inversions: the reference for
        /// the extended and cached formulas.
        fn add_ref(c: &Curve, p: &Affine, q: &Affine) -> Affine {
            let e = c.d * p.u * q.u * p.v * q.v;
            Affine {
                u: (p.u * q.v + p.v * q.u) * Invert::inv(F::ONE + e),
                v: (p.v * q.v + p.u * q.u) * Invert::inv(F::ONE - e),
            }
        }

        /// (0, -1), of order 2.
        fn t() -> Affine {
            Affine {
                u: <F as Field>::ZERO,
                v: -F::ONE,
            }
        }

        /// (i, 0) and (-i, 0), the points of order 4.
        fn order_four() -> [Affine; 2] {
            let i = OddField::sqrt(-F::ONE).unwrap();
            [i, -i].map(|u| Affine { u, v: <F as Field>::ZERO })
        }

        fn ext_eq(c: &Curve, p: &Point, q: &Affine) -> bool {
            c.to_affine(p).equals(q)
        }

        proptest! {
            #[test]
            fn decoded_points_are_on_curve((c, ps) in curve_and_points(1)) {
                prop_assert!(c.is_on_curve(&ps[0]));
                prop_assert!(c.is_on_curve(&t()));
                for q in order_four() {
                    prop_assert!(c.is_on_curve(&q));
                }
            }

            #[test]
            fn decode_is_exact(c in curve(), x in elem()) {
                // accepted iff (1 + u^2)/(1 - d u^2) is a square, but for ±i with sgn0 1
                let ratio = (F::ONE + x.square()) * Invert::inv(F::ONE - c.d * x.square());
                let odd_i = ratio.is_zero() && x.sgn0();
                prop_assert_eq!(c.decode_u128(x.pack()).is_some(), OddField::is_square(ratio) && !odd_i);
            }

            #[test]
            fn encode_roundtrip((c, ps) in curve_and_points(1)) {
                let p = ps[0];
                let q = c.decode(p.encode()).unwrap();
                prop_assert!(q.equals(&p));
                // p and p + T = (-u, -v) are one class, with one encoding
                let pt = Affine { u: -p.u, v: -p.v };
                prop_assert_eq!(pt.encode(), p.encode());
                // -p is another class but for G's points of order 2, (±i, 0)
                prop_assert_eq!(p.neg().equals(&p), p.v.is_zero());
                prop_assert_eq!(p.neg().encode() == p.encode(), p.v.is_zero());
            }

            #[test]
            fn decode_accepts_only_canonical(c in curve(), e in any::<[u8; BYTES]>()) {
                if let Some(p) = c.decode(e) {
                    prop_assert!(c.is_on_curve(&p));
                    prop_assert_eq!(p.encode(), e);
                }
            }

            #[test]
            fn torsion_in_the_quotient(c in curve()) {
                let o = Affine::IDENTITY;
                prop_assert_eq!(o.encode(), [0; BYTES]);
                prop_assert_eq!(t().encode(), [0; BYTES]);
                prop_assert!(t().equals(&o) && t().is_identity());
                prop_assert!(c.decode([0; BYTES]).unwrap().is_identity());
                // (i, 0) + T = (-i, 0): one class, of order 2 in G
                let [i, j] = order_four();
                prop_assert!(i.equals(&add_ref(&c, &j, &t())));
                prop_assert_eq!(i.encode(), j.encode());
                prop_assert!(!i.is_identity());
                prop_assert!(add_ref(&c, &i, &i).equals(&t()));
                let e = i.encode();
                prop_assert!(c.decode(e).unwrap().equals(&i));
                // the odd one of ±i is the non-canonical encoding of that class
                let odd = if i.u.sgn0() { i.u } else { j.u };
                prop_assert!(c.decode(odd.pack().to_le_bytes()).is_none());
            }

            #[test]
            fn formulas_match_affine_law((c, ps) in curve_and_points(3)) {
                let (p, q, s) = (ps[0], ps[1], ps[2]);
                let (ep, eq_) = (c.from_affine(&p), c.from_affine(&q));
                let pq = add_ref(&c, &p, &q);
                prop_assert!(c.is_on_curve(&pq));
                prop_assert!(ext_eq(&c, &c.add_ext(&ep, &eq_), &pq));
                prop_assert!(ext_eq(&c, &c.add_cached(&ep, &c.to_cached(&q)), &pq));
                // unified: the same formulas double
                prop_assert!(ext_eq(&c, &c.add_ext(&ep, &ep), &add_ref(&c, &p, &p)));
                prop_assert!(ext_eq(&c, &c.add_cached(&ep, &c.to_cached(&p)), &add_ref(&c, &p, &p)));
                let pqs = c.add_cached(&c.add_ext(&ep, &eq_), &c.to_cached(&s));
                prop_assert!(ext_eq(&c, &pqs, &add_ref(&c, &pq, &s)));
                prop_assert_eq!(pqs.t * pqs.z, pqs.x * pqs.y);
            }

            #[test]
            fn group_laws((c, ps) in curve_and_points(3)) {
                let [p, q, s] = [ps[0], ps[1], ps[2]].map(|a| c.from_affine(&a));
                let o = Point::IDENTITY;
                prop_assert!(c.add_ext(&p, &q).equals(&c.add_ext(&q, &p)));
                prop_assert!(c.add_ext(&c.add_ext(&p, &q), &s).equals(&c.add_ext(&p, &c.add_ext(&q, &s))));
                prop_assert!(c.add_ext(&p, &o).equals(&p));
                prop_assert!(c.add_ext(&p, &p.neg()).is_identity());
                prop_assert!(c.add_cached(&p, &c.to_cached(&ps[0].neg())).is_identity());
                prop_assert!(c.add_cached(&p, &c.to_cached(&ps[0]).neg()).is_identity());
                prop_assert!(c.add_ext(&p, &c.from_affine(&t())).equals(&p));
                prop_assert_eq!(p.equals(&q), ps[0].encode() == ps[1].encode());
                prop_assert_eq!(p.equals(&p.neg()), ps[0].v.is_zero());
            }

            #[test]
            fn equals_addend_matches_default((c, ps) in curve_and_points(2), l in elem().prop_filter("zero scale", |l| !l.is_zero())) {
                use crate::group::{Accumulate, Group, Negate};
                let (p, q) = (ps[0], ps[1]);
                let pq = add_ref(&c, &p, &q);
                let pt = Affine { u: -p.u, v: -p.v };
                let [i, j] = order_four();
                let pts = [p, q, pq, p.neg(), pt, t(), Affine::IDENTITY, i, j];
                let scale = |e: Point| Point { x: e.x * l, y: e.y * l, z: e.z * l, t: e.t * l };
                let (ep, eq_) = (c.from_affine(&p), c.from_affine(&q));
                // representatives by different routes and scales, and classes
                // of G met through either point: p and p + T, O and T, i and -i
                let accs = [
                    (ep, p),
                    (scale(ep), p),
                    (scale(c.add_cached(&ep, &c.to_cached(&t()))), p),
                    (c.from_affine(&pt), p),
                    (c.add_ext(&ep, &eq_), pq),
                    (scale(c.add_cached(&eq_, &c.to_cached(&p))), pq),
                    (c.from_affine(&t()), Affine::IDENTITY),
                    (Point::IDENTITY, Affine::IDENTITY),
                    (scale(c.add_ext(&ep, &ep.neg())), Affine::IDENTITY),
                    (c.from_affine(&i), i),
                    (scale(c.add_cached(&c.from_affine(&j), &c.to_cached(&t()))), i),
                ];
                for (e, x) in accs {
                    for y in &pts {
                        let a = c.to_cached(y);
                        let default = Group::is_identity(&c, &Accumulate::add(&c, &e, &Negate::neg_addend(&c, &a)));
                        prop_assert_eq!(Negate::equals_addend(&c, &e, &a), default);
                        prop_assert_eq!(e.equals(&scale(c.from_affine(y))), default);
                        prop_assert_eq!(default, x.equals(y));
                    }
                }
            }

            #[test]
            fn normalize_batch_matches((c, ps) in curve_and_points(6)) {
                let mut ext: Vec<Point> = ps.windows(2)
                    .map(|w| c.add_ext(&c.from_affine(&w[0]), &c.from_affine(&w[1]))).collect();
                ext.push(Point::IDENTITY);
                ext.push(c.from_affine(&t()));
                let got = c.normalize_batch(&ext);
                for (p, a) in ext.iter().zip(&got) {
                    let want = c.to_affine(p);
                    prop_assert!(a.u == want.u && a.v == want.v);
                }
            }

            #[test]
            fn sum_batch_matches_fold((c, ps) in curve_and_points(17), k in 0usize..17) {
                let seq = ps[..k].iter().fold(Affine::IDENTITY, |acc, p| add_ref(&c, &acc, p));
                prop_assert!(c.sum_batch(&ps[..k]).equals(&seq));
            }

            #[test]
            fn scalar_mul_is_repeated_addition((c, ps) in curve_and_points(1), k in 0u128..40) {
                let p = ps[0];
                let want = (0..k).fold(Affine::IDENTITY, |acc, _| add_ref(&c, &acc, &p));
                prop_assert!(ext_eq(&c, &c.mul(&c.from_affine(&p), k), &want));
            }

            #[test]
            fn elligator2_lands_on_curve(c in curve(), v in any::<u128>()) {
                use crate::curve::h2c::{Elligator2, Map, Montgomery};
                use crate::curve::twisted::MontgomeryModel;
                let mm = MontgomeryModel::new(c);
                // (0, 0), the Montgomery model's 2-torsion, is T
                let o = mm.rational_map(<F as Field>::ZERO, <F as Field>::ZERO);
                prop_assert!(o.u.is_zero() && o.v == -F::ONE);
                let e = Elligator2::new(mm).unwrap();
                let m = e.map(v);
                prop_assert!(c.is_on_curve(&m));
                // r and -r share x; sgn0(y) = sgn0(r) tells them apart: (u, v)
                // and (-u, v)
                let r: F = crate::curve::encoding::x(v);
                if !r.is_zero() {
                    let n = e.map((-r).pack());
                    prop_assert!(n.u == -m.u && n.v == m.v);
                }
            }

            #[test]
            fn hashes_land_on_curve(c in curve(), salt in any::<[u8; 32]>(), msg in any::<[u8; 8]>()) {
                let p = c.hash_to_curve(&Salted::new(b"test", &salt), &msg);
                prop_assert!(c.is_on_curve(&p));
                prop_assert_eq!(c.decode(p.encode()).map(|q| q.u == p.u && q.v == p.v), Some(true));
            }
        }
    };
}

pub(crate) use suite;
