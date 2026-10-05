//! The properties every instance of `binary` must have, as one suite for
//! each instance's `tests` module (`suite!(gf2_127, M127)`), with the strategies
//! other modules' tests draw curves and points from. Curves are B uniform
//! in the field, unless the family draws B from elsewhere
//! (`suite!(gf2_122, M122Gls, <strategy>)`).

macro_rules! suite {
    ($gf:ident, $m:ty) => {
        crate::curve::binary::tests::suite!(
            $gf,
            $m,
            fe().prop_filter("B != 0", |&v| v != 0)
                .prop_map(|v| Curve::new(Gf::new(v)))
        );
    };
    ($gf:ident, $m:ty, $curve:expr) => {
        use crate::curve::binary::{add_batch, candidate, sum_batch};
        use crate::ecmh::{Ecmh, digest_batch};
        use crate::field::$gf::tests::fe;
        use crate::field::batch::Invert;
        use crate::curve::binary::{Constant, Model};
        use crate::field::{Binary, Field};
        use crate::hash::Salted;
        use proptest::prelude::*;

        type M = $m;
        type Gf = crate::field::$gf::Gf;
        type Curve = crate::curve::binary::Curve<M>;
        type Affine = crate::curve::binary::Affine<M>;
        type Point = crate::curve::binary::Point<M>;

        pub fn curve() -> impl Strategy<Value = Curve> {
            $curve
        }

        /// A uniformly random x with Tr(x) = 1 lands in E\[r\] about half the time.
        pub fn point(c: Curve) -> impl Strategy<Value = Affine> {
            any::<u128>().prop_filter_map("not on curve", move |v| {
                let (x, sign) = candidate::<M>(v);
                c.decode_with_inverse(x, x.inv(), sign)
            })
        }

        pub fn curve_and_points(n: usize) -> impl Strategy<Value = (Curve, Vec<Affine>)> {
            curve().prop_flat_map(move |c| (Just(c), prop::collection::vec(point(c), n)))
        }

        proptest! {
            #[test]
            fn empty_digest_encodes_the_identity(c in curve()) {
                prop_assert_eq!(Ecmh::new(c, &[0; 32]).digest(), M::to_bytes(0));
                prop_assert_eq!(digest_batch(c, &[0; 32], &[]), M::to_bytes(0));
            }

            #[test]
            fn decoded_points_are_in_the_group((c, ps) in curve_and_points(1)) {
                let p = ps[0];
                prop_assert!(c.is_on_curve(&p) && p.x.trace() == 1);
            }

            #[test]
            fn bytes_keep_the_bits_that_fit(v in any::<u128>()) {
                let n = 8 * M::to_bytes(0).as_ref().len() as u32;
                prop_assert!(M::SIGN < n && Gf::MASK & 1 << M::SIGN == 0);
                prop_assert_eq!(M::from_bytes(M::to_bytes(v)), v & u128::MAX >> (128 - n));
            }

            #[test]
            fn encode_roundtrip((c, ps) in curve_and_points(1)) {
                let p = ps[0];
                prop_assert_eq!(c.decode(p.encode()), Some(p));
                // -P flips Tr(y), the encoding's sign bit
                let n = M::from_bytes(p.encode()) ^ 1 << M::SIGN;
                prop_assert_eq!(p.neg().encode(), M::to_bytes(n));
                prop_assert_eq!(c.decode(Affine::IDENTITY.encode()), Some(Affine::IDENTITY));
            }

            #[test]
            fn decode_accepts_only_canonical((c, v) in (curve(), any::<<M as Model>::Bytes>())) {
                if let Some(p) = c.decode(v) {
                    prop_assert!(c.is_on_curve(&p));
                    prop_assert_eq!(p.encode(), v);
                }
            }

            #[test]
            fn candidate_has_trace_one(v in any::<u128>()) {
                let (x, sign) = candidate::<M>(v);
                let y = Gf::new(v);
                prop_assert_eq!(x.trace(), 1);
                prop_assert!(Field::equals(x, y) || Field::equals(x, y + M::A));
                prop_assert_eq!(sign as u128, v >> M::SIGN & 1);
            }

            #[test]
            fn group_laws((c, ps) in curve_and_points(3)) {
                let (p, q, s) = (ps[0], ps[1], ps[2]);
                let o = Affine::IDENTITY;
                prop_assert!(c.is_on_curve(&p.add(&q)) && c.is_on_curve(&p.add(&p)));
                prop_assert_eq!(p.add(&q), q.add(&p));
                prop_assert_eq!(p.add(&q).add(&s), p.add(&q.add(&s)));
                prop_assert_eq!(p.add(&o), p);
                prop_assert_eq!(o.add(&p), p);
                prop_assert_eq!(p.add(&p.neg()), o);
                prop_assert_eq!(p.add(&p).add(&p.neg()), p);
                prop_assert!(p.add(&p).x.trace() == 1);
                let pq = p.add(&q);
                prop_assert!(pq.is_identity() || pq.x.trace() == 1);
            }

            #[test]
            fn batch_matches_single((_c, ps) in curve_and_points(24), mask in any::<u32>()) {
                // The mask mixes identity addends, doublings and cancellations
                // into ordinary additions to exercise the batch exceptions.
                let a: Vec<Affine> = ps[..12].to_vec();
                let b: Vec<Affine> = ps[12..].iter().enumerate().map(|(i, q)| match mask >> (2 * i) & 3 {
                    0 => *q,
                    1 => a[i],
                    2 => a[i].neg(),
                    _ => Affine::IDENTITY,
                }).collect();
                let mut out = vec![Affine::IDENTITY; 12];
                add_batch(&a, &b, &mut out);
                for i in 0..12 {
                    prop_assert_eq!(out[i], a[i].add(&b[i]));
                }
            }

            #[test]
            fn sum_batch_matches_fold((_c, ps) in curve_and_points(17), k in 0usize..17) {
                let seq = ps[..k].iter().fold(Affine::IDENTITY, |acc, p| acc.add(p));
                prop_assert_eq!(sum_batch(&ps[..k]), seq);
            }

            #[test]
            fn extended_roundtrip((c, ps) in curve_and_points(1)) {
                let p = ps[0];
                prop_assert_eq!(c.to_affine(&c.from_affine(&p)), p);
                prop_assert_eq!(c.to_affine(&c.neutral()), Affine::IDENTITY);
                prop_assert_eq!(c.to_affine(&c.from_affine(&p).neg()), p.neg());
            }

            #[test]
            fn extended_matches_affine((c, ps) in curve_and_points(3)) {
                let (p, q, s) = (ps[0], ps[1], ps[2]);
                let (ep, eq_) = (c.from_affine(&p), c.from_affine(&q));
                prop_assert_eq!(c.to_affine(&c.add(&ep, &eq_)), p.add(&q));
                prop_assert_eq!(c.to_affine(&c.add_affine(&ep, &q)), p.add(&q));
                // Exceptional inputs of the affine law are ordinary here:
                // P = Q, P = -Q and O, each way round.
                prop_assert_eq!(c.to_affine(&c.add(&ep, &ep)), p.add(&p));
                prop_assert_eq!(c.to_affine(&c.add_affine(&ep, &p)), p.add(&p));
                prop_assert!(c.add(&ep, &ep.neg()).is_identity());
                prop_assert_eq!(c.to_affine(&c.add_affine(&ep, &p.neg())), Affine::IDENTITY);
                prop_assert_eq!(c.to_affine(&c.add(&c.neutral(), &ep)), p);
                prop_assert_eq!(c.to_affine(&c.add(&ep, &c.neutral())), p);
                prop_assert!(c.add(&c.neutral(), &c.neutral()).is_identity());
                prop_assert_eq!(c.to_affine(&c.add_affine(&c.neutral(), &p)), p);
                prop_assert_eq!(c.to_affine(&c.add_affine(&ep, &Affine::IDENTITY)), p);
                let pqs = c.add_affine(&c.add(&ep, &eq_), &s);
                prop_assert_eq!(c.to_affine(&pqs), p.add(&q).add(&s));
            }

            #[test]
            fn associative((c, ps) in curve_and_points(3)) {
                let [p, q, s] = [ps[0], ps[1], ps[2]].map(|p| c.from_affine(&p));
                let l = c.add(&c.add(&p, &q), &s);
                let r = c.add(&p, &c.add(&q, &s));
                prop_assert!(l.equals(&r));
                prop_assert!(c.add(&p, &q).equals(&c.add(&q, &p)));
            }

            #[test]
            fn chains_match_affine((c, ps) in curve_and_points(6), ops in prop::collection::vec((0usize..6, 0u8..4), 0..32)) {
                // Ops 2 and 3 add or subtract the running sum, exercising
                // doubling and cancellation after arbitrary preceding additions.
                let (mut acc, mut want) = (c.neutral(), Affine::IDENTITY);
                for (i, op) in ops {
                    let a = match op {
                        0 => ps[i],
                        1 => ps[i].neg(),
                        2 => want,
                        _ => want.neg(),
                    };
                    acc = c.add(&acc, &c.from_affine(&a));
                    want = want.add(&a);
                    prop_assert_eq!(c.to_affine(&acc), want);
                    prop_assert_eq!(acc.is_identity(), want.is_identity());
                }
            }

            #[test]
            fn extended_equality((c, ps) in curve_and_points(2)) {
                let (p, q) = (ps[0], ps[1]);
                let (ep, eq_) = (c.from_affine(&p), c.from_affine(&q));
                // Same point, different projective representatives.
                prop_assert!(c.add(&ep, &eq_).equals(&c.add_affine(&eq_, &p)));
                prop_assert!(!ep.equals(&eq_) && !ep.equals(&ep.neg()) && !ep.equals(&c.neutral()));
                prop_assert!(c.neutral().equals(&c.add_affine(&ep, &p.neg())));
            }

            #[test]
            fn normalize_batch_matches((c, ps) in curve_and_points(6)) {
                let mut ext: Vec<Point> = ps.windows(2).map(|w| c.add_affine(&c.from_affine(&w[0]), &w[1])).collect();
                ext.push(c.neutral());
                let want: Vec<Affine> = ext.iter().map(|p| c.to_affine(p)).collect();
                prop_assert_eq!(c.normalize_batch(&ext), want);
            }

            #[test]
            fn double_matches_add((c, ps) in curve_and_points(2)) {
                let (p, q) = (ps[0], ps[1]);
                let n = c.neutral();
                prop_assert!(c.double(&n).is_identity() && c.double(&n).equals(&n));
                let pq = c.add(&c.from_affine(&p), &c.from_affine(&q));
                for x in [c.from_affine(&p), pq, c.double(&pq)] {
                    let d = c.double(&x);
                    prop_assert!(d.equals(&c.add(&x, &x)));
                    prop_assert_eq!(c.to_affine(&d), c.to_affine(&x).add(&c.to_affine(&x)));
                    prop_assert!(Field::equals(d.t, d.x * d.z));
                }
            }

            #[test]
            fn scalar_mul_is_repeated_addition((c, ps) in curve_and_points(1), k in 0u128..40) {
                let p = ps[0];
                let want = (0..k).fold(Affine::IDENTITY, |acc, _| acc.add(&p));
                prop_assert_eq!(c.to_affine(&c.mul(&c.from_affine(&p), k)), want);
            }

            #[test]
            fn hash_to_curve_lands_in_group(c in curve(), salt in any::<[u8; 32]>(), msg in any::<[u8; 8]>()) {
                let h = Salted::new(b"test", &salt);
                let p = c.hash_to_curve(&h, &msg);
                prop_assert!(!p.is_identity() && c.is_on_curve(&p) && p.x.trace() == 1);
                prop_assert_eq!(p, c.hash_to_curve(&h, &msg));
            }

            #[test]
            fn x_on_curve_matches_decode(c in curve(), v in fe()) {
                // 0 is O, so leave it out
                prop_assume!(v != 0);
                prop_assert_eq!(c.x_on_curve(Gf::new(v)), c.decode(M::to_bytes(v)).is_some());
            }

            #[test]
            fn batched_hashing_matches_single(c in curve(), salt in any::<[u8; 32]>(), msgs in prop::collection::vec(prop::collection::vec(any::<u8>(), 0..8), 0..20)) {
                let h = Salted::new(b"test", &salt);
                let refs: Vec<&[u8]> = msgs.iter().map(|m| m.as_slice()).collect();
                let single: Vec<Affine> = refs.iter().map(|m| c.hash_to_curve(&h, m)).collect();
                prop_assert_eq!(c.hash_to_curve_batch(&h, &refs), single);
            }

            #[test]
            fn ecmh_streaming_matches_batch(c in curve(), salt in any::<[u8; 32]>(), xs in prop::collection::vec(prop::collection::vec(any::<u8>(), 0..8), 0..12)) {
                let refs: Vec<&[u8]> = xs.iter().map(|m| m.as_slice()).collect();
                let mut e = Ecmh::new(c, &salt);
                refs.iter().for_each(|x| e.insert(x));
                prop_assert_eq!(e.digest(), digest_batch(c, &salt, &refs));
                refs.iter().rev().for_each(|x| e.remove(x));
                prop_assert_eq!(e.digest(), M::to_bytes(0));
            }
        }

        /// `binary::lambda` against the affine law.
        mod lambda {
            use super::*;
            use crate::curve::binary::lambda::{lift_batch, normalize_batch};

            type LAffine = crate::curve::binary::lambda::Affine<M>;
            type LPoint = crate::curve::binary::lambda::Point<M>;

            fn lam(p: &Affine) -> LAffine {
                LAffine::lift(p)
            }

            proptest! {
                #[test]
                fn lift_roundtrip((c, ps) in curve_and_points(1)) {
                    let p = ps[0];
                    prop_assert_eq!(LPoint::from(lam(&p)).to_affine(), p);
                    prop_assert_eq!(LPoint::from(lam(&p).neg()).to_affine(), p.neg());
                    prop_assert!(lam(&Affine::IDENTITY).is_identity());
                    prop_assert_eq!(LPoint::IDENTITY.to_affine(), Affine::IDENTITY);
                    prop_assert!(LAffine::IDENTITY.double().is_identity());
                    // on the λ form of the curve equation
                    let (x, l) = (p.x, lam(&p).l);
                    prop_assert!(Field::equals((l.square() + l + M::A) * x.square(), x.square().square() + c.big_b.gf()));
                }

                #[test]
                fn lift_batch_matches_single((_c, mut ps) in curve_and_points(6)) {
                    ps.push(Affine::IDENTITY);
                    let want: Vec<LPoint> = ps.iter().map(|p| lam(p).into()).collect();
                    let got: Vec<LPoint> = lift_batch(&ps).into_iter().map(LPoint::from).collect();
                    for (w, g) in want.iter().zip(&got) {
                        prop_assert!(w.equals(g));
                    }
                }

                #[test]
                fn mixed_matches_affine((_c, ps) in curve_and_points(3)) {
                    let (p, q, s) = (ps[0], ps[1], ps[2]);
                    let o = Affine::IDENTITY;
                    let lp = LPoint::from(lam(&p));
                    prop_assert_eq!(lp.add_affine(&lam(&q)).to_affine(), p.add(&q));
                    prop_assert_eq!(LPoint::IDENTITY.add_affine(&lam(&q)).to_affine(), q);
                    prop_assert_eq!(lp.add_affine(&lam(&o)).to_affine(), p);
                    prop_assert!(LPoint::IDENTITY.add_affine(&lam(&o)).is_identity());
                    prop_assert_eq!(lp.add_affine(&lam(&p)).to_affine(), p.add(&p));
                    prop_assert!(lp.add_affine(&lam(&p).neg()).is_identity());
                    prop_assert_eq!(lam(&p).double().to_affine(), p.add(&p));
                    // the exceptional cases again with Z != 1
                    let pq = lp.add_affine(&lam(&q));
                    let spq = p.add(&q);
                    prop_assert_eq!(pq.add_affine(&lam(&spq)).to_affine(), spq.add(&spq));
                    prop_assert!(pq.add_affine(&lam(&spq).neg()).is_identity());
                    prop_assert_eq!(pq.add_affine(&lam(&s)).to_affine(), spq.add(&s));
                }

                #[test]
                fn chains_match_affine((_c, ps) in curve_and_points(8), ops in prop::collection::vec((0usize..8, 0u8..4), 0..48)) {
                    // Ops 2 and 3 add or subtract the running sum, exercising
                    // doubling and cancellation after arbitrary preceding additions.
                    let (mut acc, mut want) = (LPoint::IDENTITY, Affine::IDENTITY);
                    for (i, op) in ops {
                        let a = match op {
                            0 => ps[i],
                            1 => ps[i].neg(),
                            2 => want,
                            _ => want.neg(),
                        };
                        acc = acc.add_affine(&lam(&a));
                        want = want.add(&a);
                        prop_assert_eq!(acc.to_affine(), want);
                        prop_assert_eq!(acc.is_identity(), want.is_identity());
                    }
                }

                #[test]
                fn equality((_c, ps) in curve_and_points(2)) {
                    let (p, q) = (ps[0], ps[1]);
                    let (lp, lq) = (LPoint::from(lam(&p)), LPoint::from(lam(&q)));
                    // same point, different representatives
                    prop_assert!(lp.add_affine(&lam(&q)).equals(&lq.add_affine(&lam(&p))));
                    // p and q may coincide; p and -p never do, in odd order
                    prop_assert_eq!(lp.equals(&lq), p == q);
                    prop_assert!(!lp.equals(&lam(&p).neg().into()) && !lp.equals(&LPoint::IDENTITY));
                    prop_assert!(LPoint::IDENTITY.equals(&lp.add_affine(&lam(&p).neg())));
                }

                #[test]
                fn equals_affine_is_affine_equality((_c, ps) in curve_and_points(3)) {
                    // accumulators at Z = 1 and Z != 1 against each addend, O included
                    let (p, q, r) = (ps[0], ps[1], ps[2]);
                    let pq = LPoint::from(lam(&p)).add_affine(&lam(&q));
                    let cancelled = pq.add_affine(&lam(&q).neg());
                    let o = Affine::IDENTITY;
                    for acc in [LPoint::from(lam(&p)), pq, cancelled, LPoint::IDENTITY] {
                        for a in [p, q, p.add(&q), r, p.neg(), o] {
                            prop_assert_eq!(acc.equals_affine(&lam(&a)), acc.to_affine() == a);
                        }
                    }
                }

                #[test]
                fn normalize_batch_matches((_c, ps) in curve_and_points(6)) {
                    let mut pts: Vec<LPoint> = ps.windows(2).map(|w| LPoint::from(lam(&w[0])).add_affine(&lam(&w[1]))).collect();
                    pts.push(LPoint::IDENTITY);
                    let want: Vec<Affine> = pts.iter().map(|p| p.to_affine()).collect();
                    prop_assert_eq!(normalize_batch(&pts), want);
                }

                #[test]
                fn hash_to_lambda_is_the_lift((c, _ps) in curve_and_points(0), salt in any::<[u8; 32]>(), msgs in prop::collection::vec(prop::collection::vec(any::<u8>(), 0..8), 0..12)) {
                    let h = Salted::new(b"test", &salt);
                    let refs: Vec<&[u8]> = msgs.iter().map(|m| m.as_slice()).collect();
                    let want: Vec<LAffine> = refs.iter().map(|m| lam(&c.hash_to_curve(&h, m))).collect();
                    let one: Vec<LAffine> = refs.iter().map(|m| c.hash_to_lambda(&h, m)).collect();
                    prop_assert_eq!(&one, &want);
                    prop_assert_eq!(c.hash_to_lambda_batch(&h, &refs), want);
                }
            }
        }

        /// `binary::unscaled`, with points checked through the (x, y) law.
        mod unscaled {
            use super::*;
            use crate::group::{Accumulate, Encode, Group, SumBatch};

            type UCurve = crate::curve::binary::unscaled::Curve<M>;
            type UAffine = crate::curve::binary::unscaled::Affine<M>;
            type UPoint = crate::curve::binary::unscaled::Point<M>;

            /// The (x, y) point an addend stands for: x = 1/u and
            /// y = (v + 1) x^2 + (1 + a) x.
            fn xy(a: &UAffine) -> Affine {
                if a.is_identity() {
                    return Affine::IDENTITY;
                }
                let x = a.u.inv();
                Affine { x, y: (a.v + Gf::ONE) * x.square() + (Gf::ONE + M::A) * x }
            }

            fn at(u: &UCurve, p: &UPoint) -> Affine {
                xy(&u.to_affine(p))
            }

            proptest! {
                #[test]
                fn lift_roundtrip((_c, ps) in curve_and_points(1)) {
                    let p = ps[0];
                    prop_assert_eq!(xy(&UAffine::lift(&p)), p);
                    prop_assert_eq!(xy(&UAffine::lift(&p).neg()), p.neg());
                    prop_assert_eq!(UAffine::lift(&Affine::IDENTITY), UAffine::IDENTITY);
                }

                #[test]
                fn add_matches_affine((c, ps) in curve_and_points(3)) {
                    let u = UCurve::new(c);
                    let (p, q, s) = (ps[0], ps[1], ps[2]);
                    let (lp, lq) = (UAffine::lift(&p), UAffine::lift(&q));
                    let o = UAffine::IDENTITY;
                    let n = u.neutral();
                    let pp = n.add(&lp);
                    prop_assert_eq!(at(&u, &pp), p);
                    prop_assert_eq!(at(&u, &pp.add(&lq)), p.add(&q));
                    // N, P = Q and P = -Q
                    prop_assert!(n.add(&o).is_identity());
                    prop_assert_eq!(at(&u, &pp.add(&o)), p);
                    prop_assert_eq!(at(&u, &pp.add(&lp)), p.add(&p));
                    prop_assert!(pp.add(&lp.neg()).is_identity());
                    // the same again with Z != 1
                    let pq = pp.add(&lq);
                    let spq = UAffine::lift(&p.add(&q));
                    prop_assert_eq!(at(&u, &pq.add(&spq)), p.add(&q).add(&p.add(&q)));
                    prop_assert!(pq.add(&spq.neg()).is_identity());
                    prop_assert_eq!(at(&u, &pq.add(&o)), p.add(&q));
                    prop_assert_eq!(at(&u, &pq.add(&UAffine::lift(&s))), p.add(&q).add(&s));
                    prop_assert!(!pq.z.is_zero() && !pq.add(&spq.neg()).z.is_zero());
                }

                #[test]
                fn chains_match_affine((c, ps) in curve_and_points(8), ops in prop::collection::vec((0usize..8, 0u8..4), 0..48)) {
                    // Ops 2 and 3 add or subtract the running sum, exercising
                    // doubling and cancellation after arbitrary preceding additions.
                    let u = UCurve::new(c);
                    let (mut acc, mut want) = (u.neutral(), Affine::IDENTITY);
                    for (i, op) in ops {
                        let a = match op {
                            0 => ps[i],
                            1 => ps[i].neg(),
                            2 => want,
                            _ => want.neg(),
                        };
                        acc = u.add(&acc, &UAffine::lift(&a));
                        want = want.add(&a);
                        prop_assert_eq!(at(&u, &acc), want);
                        prop_assert_eq!(acc.is_identity(), want.is_identity());
                        prop_assert!(acc.equals_affine(&UAffine::lift(&want)));
                    }
                }

                #[test]
                fn equality((c, ps) in curve_and_points(2)) {
                    let u = UCurve::new(c);
                    let (lp, lq) = (UAffine::lift(&ps[0]), UAffine::lift(&ps[1]));
                    let o = UAffine::IDENTITY;
                    let (n, pp, qq) = (u.neutral(), u.neutral().add(&lp), u.neutral().add(&lq));
                    let (pq, qp) = (pp.add(&lq), qq.add(&lp));
                    // same point, different representatives
                    prop_assert!(pq.equals(&qp) && pq.equals_affine(&UAffine::lift(&ps[0].add(&ps[1]))));
                    // p and q may coincide; p and -p never do, in odd order
                    let same = ps[0] == ps[1];
                    prop_assert_eq!(pp.equals(&qq), same);
                    prop_assert_eq!(pp.equals_affine(&lq), same);
                    prop_assert!(pp.equals(&n.add(&lp)) && pp.equals_affine(&lp));
                    prop_assert!(!pp.equals(&n) && !pp.equals(&n.add(&lp.neg())));
                    prop_assert!(!pp.equals_affine(&lp.neg()) && !pp.equals_affine(&o));
                    prop_assert!(n.equals(&pp.add(&lp.neg())) && n.equals_affine(&o) && !n.equals_affine(&lp));
                    prop_assert!(pp.add(&lp.neg()).equals_affine(&o));
                }

                #[test]
                fn to_affine_batch_matches((c, ps) in curve_and_points(6)) {
                    let u = UCurve::new(c);
                    let mut pts: Vec<UPoint> = ps.windows(2).map(|w| u.neutral().add(&UAffine::lift(&w[0])).add(&UAffine::lift(&w[1]))).collect();
                    pts.push(u.neutral());
                    let want: Vec<UAffine> = pts.iter().map(|p| u.to_affine(p)).collect();
                    prop_assert_eq!(u.to_affine_batch(&pts), want);
                    let lifted: Vec<UAffine> = ps.iter().map(UAffine::lift).collect();
                    let sum = ps.iter().fold(Affine::IDENTITY, |s, p| s.add(p));
                    prop_assert_eq!(xy(&u.sum_batch(&lifted)), sum);
                }

                #[test]
                fn codec_is_wcodecs((c, ps) in curve_and_points(3)) {
                    let u = UCurve::new(c);
                    let w = &u.w;
                    let lam = crate::curve::binary::lambda::Affine::<M>::lift;
                    let a = UAffine::lift(&ps[0]);
                    let e = u.encode(&a);
                    prop_assert_eq!(e, w.encode(&lam(&ps[0])));
                    prop_assert_eq!(u.decode(e), Some(a));
                    prop_assert_eq!(u.encode(&a.neg()), w.encode(&lam(&ps[0].neg())));
                    prop_assert_eq!(u.encode(&UAffine::IDENTITY), M::to_bytes(0));
                    prop_assert_eq!(u.decode(M::to_bytes(0)), Some(UAffine::IDENTITY));
                    // from accumulators, one at a time or batched, with O among them
                    let pts = [u.neutral().add(&a), u.neutral().add(&a).add(&UAffine::lift(&ps[1])), u.neutral(), u.neutral().add(&a).add(&a.neg())];
                    let want: Vec<_> = pts.iter().map(|p| u.encode(&u.to_affine(p))).collect();
                    let direct: Vec<_> = pts.iter().map(|p| u.encode_point(p)).collect();
                    prop_assert_eq!(direct[0], e);
                    prop_assert_eq!(direct[1], w.encode(&lam(&ps[0].add(&ps[1]))));
                    prop_assert_eq!(&direct, &want);
                    prop_assert_eq!(Encode::encode_batch(&u, &pts), want);
                }

                #[test]
                fn decode_accepts_only_canonical((c, v) in (curve(), any::<u128>())) {
                    let u = UCurve::new(c);
                    let e = M::to_bytes(v & (Gf::MASK | 1 << M::SIGN));
                    let got = u.decode(e);
                    prop_assert_eq!(got.map(|a| xy(&a)), u.w.decode(e).map(|a| crate::curve::binary::lambda::Point::from(a).to_affine()));
                    if let Some(a) = got {
                        prop_assert_eq!(u.encode(&a), e);
                    }
                }

                #[test]
                fn decode_batch_matches_single((c, ps) in curve_and_points(4), junk in prop::collection::vec(any::<u128>(), 0..4)) {
                    let u = UCurve::new(c);
                    let mut es: Vec<_> = ps.iter().map(|p| u.encode(&UAffine::lift(p))).collect();
                    es.push(M::to_bytes(0));
                    es.extend(junk.into_iter().map(M::to_bytes));
                    let want: Vec<_> = es.iter().map(|&e| u.decode(e)).collect();
                    prop_assert_eq!(u.decode_batch(&es), want);
                }

                #[test]
                fn hash_is_wcodecs((c, msgs) in (curve(), prop::collection::vec(any::<[u8; 4]>(), 0..12))) {
                    let u = UCurve::new(c);
                    let h = Salted::new(b"wcodec test", &[0; 32]);
                    let refs: Vec<&[u8]> = msgs.iter().map(|m| &m[..]).collect();
                    let to_xy = |a: &crate::curve::binary::lambda::Affine<M>| crate::curve::binary::lambda::Point::from(*a).to_affine();
                    let want: Vec<_> = refs.iter().map(|m| to_xy(&u.w.hash(&h, m))).collect();
                    let one: Vec<_> = refs.iter().map(|m| xy(&crate::group::HashToCurve::hash(&u, &h, m))).collect();
                    let batch: Vec<_> = crate::group::HashToCurve::hash_batch(&u, &h, &refs).iter().map(xy).collect();
                    prop_assert_eq!(&one, &want);
                    prop_assert_eq!(&batch, &want);
                }
            }
        }

        /// `binary::wcodec`, with points checked through the (x, y) law.
        mod wcodec {
            use super::*;
            use crate::curve::binary::wcodec::to_lambda_batch;

            type WCurve = crate::curve::binary::wcodec::Curve<M>;
            type WAffine = crate::curve::binary::lambda::Affine<M>;
            type WPoint = crate::curve::binary::lambda::Point<M>;

            /// The (x, y) point a λ-affine one stands for.
            fn xy(a: &WAffine) -> Affine {
                WPoint::from(*a).to_affine()
            }

            proptest! {
                #[test]
                fn roundtrip((c, ps) in curve_and_points(1)) {
                    let w = WCurve::new(c);
                    let a = WAffine::lift(&ps[0]);
                    let e = w.encode(&a);
                    prop_assert_eq!(w.decode(e).map(|b| xy(&b)), Some(ps[0]));
                    // -P is w + 1
                    let mut n = e;
                    n[0] ^= 1;
                    prop_assert_eq!(w.encode(&a.neg()), n);
                    prop_assert_eq!(w.encode(&WAffine::IDENTITY), M::to_bytes(0));
                    prop_assert!(w.decode(M::to_bytes(0)).unwrap().is_identity());
                }

                #[test]
                fn decode_accepts_only_canonical((c, v) in (curve(), any::<u128>())) {
                    let w = WCurve::new(c);
                    // the bits a point encoding may have, of which the sign
                    // bit is outside the mask and nothing encodes
                    let v = v & (Gf::MASK | 1 << M::SIGN);
                    let e = M::to_bytes(v);
                    if let Some(a) = w.decode(e) {
                        let p = xy(&a);
                        prop_assert!(c.is_on_curve(&p) && (p.is_identity() || p.x.trace() == 1));
                        prop_assert_eq!(w.encode(&a), e);
                        prop_assert!(v & !Gf::MASK == 0);
                    }
                }

                #[test]
                fn hash_lands_in_the_group((c, msg) in (curve(), any::<[u8; 8]>())) {
                    let w = WCurve::new(c);
                    let h = Salted::new(b"wcodec test", &[0; 32]);
                    let a = w.hash(&h, &msg);
                    let p = xy(&a);
                    prop_assert!(!a.is_identity() && c.is_on_curve(&p) && p.x.trace() == 1);
                    prop_assert_eq!(w.decode(w.encode(&a)).map(|b| xy(&b)), Some(p));
                }

                #[test]
                fn hash_batch_matches_single((c, msgs) in (curve(), prop::collection::vec(any::<[u8; 4]>(), 0..12))) {
                    let w = WCurve::new(c);
                    let h = Salted::new(b"wcodec test", &[0; 32]);
                    let refs: Vec<&[u8]> = msgs.iter().map(|m| &m[..]).collect();
                    let got: Vec<_> = w.hash_batch(&h, &refs).iter().map(xy).collect();
                    let want: Vec<_> = refs.iter().map(|m| xy(&w.hash(&h, m))).collect();
                    prop_assert_eq!(got, want);
                }

                #[test]
                fn to_lambda_matches((c, ps) in curve_and_points(6)) {
                    let w = WCurve::new(c);
                    let mut pts: Vec<WPoint> = ps
                        .windows(2)
                        .map(|p| WPoint::from(WAffine::lift(&p[0])).add_affine(&WAffine::lift(&p[1])))
                        .collect();
                    pts.push(WPoint::IDENTITY);
                    let want: Vec<_> = pts.iter().map(|p| p.to_affine()).collect();
                    let single: Vec<_> = pts.iter().map(|p| xy(&p.to_lambda())).collect();
                    let batch: Vec<_> = to_lambda_batch(&pts).iter().map(xy).collect();
                    prop_assert_eq!(&single, &want);
                    prop_assert_eq!(&batch, &want);
                    let lam: Vec<WAffine> = ps.iter().map(WAffine::lift).collect();
                    let sum = ps.iter().fold(Affine::IDENTITY, |s, p| s.add(p));
                    prop_assert_eq!(xy(&crate::group::SumBatch::sum_batch(&w, &lam)), sum);
                }
            }

            proptest! {
                #[test]
                fn w0_decodes_to_nothing(c in curve()) {
                    // Tr(b/d^2) = 1, with w0 = 0 iff Tr(b/a^2) = 1, as d(0) = a
                    let w = WCurve::new(c);
                    let d = w.w0.square() + w.w0 + M::A;
                    prop_assert_eq!((c.b.gf() * d.square().inv()).trace(), 1);
                }
            }
        }
    };
}
pub(crate) use suite;
