//! Shared property tests and curve/point strategies for `edwards` instances.
//! `suite!(fp127, F, fp)` takes the field module, type and element strategy.

macro_rules! suite {
    ($m:ident, $f:ident, $elem:ident) => {
        use crate::curve::h2c::Elligator2;
        use crate::curve::encoding::{Signed, mask};
        use crate::curve::h2c::Map;
        use crate::field::$m::tests::$elem as elem;
        use crate::field::{OddField, Packed};
        use crate::hash::Salted;
        use proptest::prelude::*;

        type F = crate::field::$m::$f;
        type Curve = crate::curve::edwards::Curve<F>;
        type Affine = crate::curve::edwards::Affine<F>;
        type Point = crate::curve::edwards::Point<F>;

        pub fn curve() -> impl Strategy<Value = Curve> {
            elem().prop_filter_map("Edwards model incomplete", Curve::new)
        }

        pub fn point(c: Curve) -> impl Strategy<Value = Affine> {
            any::<u128>().prop_filter_map("not on curve", move |v| {
                c.decode(F::to_bytes(v & mask::<F>())).filter(|p| !p.is_identity())
            })
        }

        pub fn curve_and_points(n: usize) -> impl Strategy<Value = (Curve, Vec<Affine>)> {
            curve().prop_flat_map(move |c| (Just(c), prop::collection::vec(point(c), n)))
        }

        proptest! {
            #[test]
            fn models_agree(c in curve()) {
                // Montgomery B y^2 = x^3 + A x^2 + x scaled by 1/B, A = 2(1+d)/(1-d)
                let (e, two) = (F::ONE - c.d, F::ONE + F::ONE);
                let a = two * (F::ONE + c.d) * e.invert();
                prop_assert_eq!(c.a2 * c.bm, a);
                prop_assert_eq!(c.a4 * c.bm.square(), F::ONE);
                // x^2 + a2 x + a4 has discriminant d: (0, 0) is the only 2-torsion
                prop_assert_eq!(c.a2.square() - two * two * c.a4, c.d);
            }

            #[test]
            fn edwards_points_satisfy_edwards_equation((c, ps) in curve_and_points(1)) {
                let e = c.from_affine(&ps[0]);
                let zi = e.z.invert();
                let (u, v) = (e.x * zi, e.y * zi);
                prop_assert_eq!(u.square() + v.square(), F::ONE + c.d * u.square() * v.square());
                prop_assert_eq!(e.t * zi, u * v);
            }

            #[test]
            fn edwards_native_hash(c in curve(), salt in any::<[u8; 32]>(), msg in prop::collection::vec(any::<u8>(), 0..8)) {
                let q = c.hash_to_edwards(&Salted::new(b"test", &salt), &msg);
                let (u2, v2) = (q.u.square(), q.v.square());
                prop_assert_eq!(u2 + v2, F::ONE + c.d * u2 * v2);
                prop_assert_eq!(q.duv, c.d * q.u * q.v);
                let m = c.to_affine(&c.add_cached(&Point::IDENTITY, &q));
                prop_assert!(c.is_on_curve(&m));
                prop_assert_eq!(c.to_cached_batch(&[m])[0].v, q.v);
            }

            #[test]
            fn coordinate_checks_match_decoding(c in curve(), x in elem()) {
                let v = x.pack();
                prop_assert_eq!(c.x_on_curve(x), c.decode(F::to_bytes(v)).is_some());
                prop_assert_eq!(c.u_on_curve(x), c.edwards_from_u(v).is_some());
            }

            #[test]
            fn edwards_from_u_is_exact(c in curve(), v in any::<u128>()) {
                // accepted iff (1 - u^2)/(1 - d u^2) is a square, with u the
                // low BITS bits of v, a mask written independently of encoding's
                let u = F::unpack(v & ((1 << F::BITS) - 1));
                prop_assume!(u.is_some());
                let u = u.unwrap();
                let ratio = (F::ONE - u.square()) * (F::ONE - c.d * u.square()).invert();
                prop_assert_eq!(c.edwards_from_u(v).is_some(), ratio.is_square());
            }

            #[test]
            fn encode_roundtrip((c, ps) in curve_and_points(1)) {
                let p = ps[0];
                prop_assert!(c.is_on_curve(&p));
                prop_assert_eq!(c.decode(p.encode()), Some(p));
                prop_assert_eq!(c.decode(Affine::IDENTITY.encode()), Some(Affine::IDENTITY));
                let t = Affine { x: F::ZERO, y: F::ZERO };
                prop_assert_eq!(c.decode(t.encode()), Some(t));
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
                let t = Affine { x: F::ZERO, y: F::ZERO };
                prop_assert!(c.is_on_curve(&c.add(&p, &q)));
                prop_assert_eq!(c.add(&p, &q), c.add(&q, &p));
                prop_assert_eq!(c.add(&c.add(&p, &q), &s), c.add(&p, &c.add(&q, &s)));
                prop_assert_eq!(c.add(&p, &o), p);
                prop_assert_eq!(c.add(&p, &p.neg()), o);
                prop_assert_eq!(c.add(&c.add(&p, &p), &p.neg()), p);
                prop_assert_eq!(c.add(&t, &t), o);
                prop_assert_eq!(c.add(&c.add(&p, &t), &t), p);
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
            fn edwards_matches_montgomery((c, ps) in curve_and_points(3)) {
                let (p, q, s) = (ps[0], ps[1], ps[2]);
                let t = Affine { x: F::ZERO, y: F::ZERO };
                let (ep, eq_) = (c.from_affine(&p), c.from_affine(&q));
                for a in [p, t, Affine::IDENTITY] {
                    prop_assert_eq!(c.to_affine(&c.from_affine(&a)), a);
                }
                prop_assert_eq!(c.to_affine(&ep.neg()), p.neg());
                prop_assert_eq!(c.to_affine(&c.add_ext(&ep, &eq_)), c.add(&p, &q));
                prop_assert_eq!(c.to_affine(&c.add_ext(&ep, &ep)), c.add(&p, &p));
                prop_assert_eq!(c.to_affine(&c.add_ext(&ep, &ep.neg())), Affine::IDENTITY);
                let cached = c.to_cached_batch(&[q, s, p.neg()]);
                prop_assert_eq!(c.to_affine(&c.add_cached(&ep, &cached[0])), c.add(&p, &q));
                let pqs = c.add_cached(&c.add_cached(&ep, &cached[0]), &cached[1]);
                prop_assert_eq!(c.to_affine(&pqs), c.add(&c.add(&p, &q), &s));
                prop_assert!(c.add_cached(&ep, &cached[2]).equals(&Point::IDENTITY));
                prop_assert!(c.add_ext(&ep, &eq_).equals(&c.add_cached(&eq_, &c.to_cached_batch(&[p])[0])));
                prop_assert!(!ep.equals(&eq_) && !ep.equals(&ep.neg()));
            }

            #[test]
            fn equals_addend_matches_default((c, ps) in curve_and_points(2), l in elem().prop_filter("zero scale", |l| !l.is_zero())) {
                use crate::group::{Accumulate, Group, Negate};
                let (p, q) = (ps[0], ps[1]);
                let t = Affine { x: F::ZERO, y: F::ZERO };
                let pq = c.add(&p, &q);
                let pts = [p, q, pq, p.neg(), t, c.add(&p, &t), Affine::IDENTITY];
                let adds = c.to_cached_batch(&pts);
                let scale = |e: Point| Point { x: e.x * l, y: e.y * l, z: e.z * l, t: e.t * l };
                let (ep, eq_) = (c.from_affine(&p), c.from_affine(&q));
                // representatives of one point by different routes and scales
                let accs = [
                    (ep, p),
                    (scale(ep), p),
                    (c.add_ext(&ep, &eq_), pq),
                    (scale(c.add_cached(&eq_, &adds[0])), pq),
                    (c.from_affine(&t), t),
                    (scale(c.add_ext(&ep, &c.from_affine(&t))), c.add(&p, &t)),
                    (Point::IDENTITY, Affine::IDENTITY),
                    (scale(c.add_ext(&ep, &ep.neg())), Affine::IDENTITY),
                ];
                for (e, x) in accs {
                    for (a, y) in adds.iter().zip(&pts) {
                        let default = Group::is_identity(&c, &Accumulate::add(&c, &e, &Negate::neg_addend(&c, a)));
                        prop_assert_eq!(Negate::equals_addend(&c, &e, a), default);
                        prop_assert_eq!(default, x == *y);
                    }
                }
            }

            #[test]
            fn normalize_batch_matches((c, ps) in curve_and_points(6)) {
                let mut ext: Vec<Point> = ps.windows(2)
                    .map(|w| c.add_ext(&c.from_affine(&w[0]), &c.from_affine(&w[1]))).collect();
                ext.push(Point::IDENTITY);
                ext.push(c.from_affine(&Affine { x: F::ZERO, y: F::ZERO }));
                let want: Vec<Affine> = ext.iter().map(|p| c.to_affine(p)).collect();
                prop_assert_eq!(c.normalize_batch(&ext), want);
            }

            #[test]
            fn scalar_mul_is_repeated_addition((c, ps) in curve_and_points(1), k in 0u128..40) {
                let p = ps[0];
                let want = (0..k).fold(Affine::IDENTITY, |acc, _| c.add(&acc, &p));
                prop_assert_eq!(c.to_affine(&c.mul(&c.from_affine(&p), k)), want);
            }

            #[test]
            fn hashes_land_on_curve(c in curve(), salt in any::<[u8; 32]>(), msg in any::<[u8; 8]>(), v in any::<u128>()) {
                let h = Salted::new(b"test", &salt);
                let p = c.hash_to_curve(&h, &msg);
                prop_assert!(!p.is_identity() && c.is_on_curve(&p));
                // The strategy favours d = -1, where Montgomery a2 vanishes
                // and Elligator2 is unavailable; see `elligator2_needs_a2`.
                if let Some(e) = Elligator2::new(c) {
                    let m = e.map(v);
                    prop_assert!(c.is_on_curve(&m));
                    // r and -r share x; sgn0(y) = sgn0(r) tells them apart
                    let r: F = crate::curve::encoding::x(v);
                    prop_assert!(m.y.is_zero() || m.y.sgn0() == r.sgn0());
                    if !r.is_zero() {
                        prop_assert_eq!(e.map((-r).pack()), m.neg());
                    }
                }
            }
        }

        #[test]
        fn elligator2_needs_a2() {
            // d = -1 is complete iff -1 is a non-square (q = 3 mod 4), with a2 = 0.
            match Curve::new(-F::ONE) {
                Some(c) => assert!(c.a2.is_zero() && Elligator2::new(c).is_none()),
                None => assert!(OddField::is_square(-F::ONE)),
            }
        }
    };
}

pub(crate) use suite;
