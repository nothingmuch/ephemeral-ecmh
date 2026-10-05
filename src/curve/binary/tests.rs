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
        use crate::curve::binary::Model;
        use crate::curve::binary::sum_batch;
        use crate::field::Binary;
        use crate::field::$gf::tests::fe;
        use proptest::prelude::*;

        type M = $m;
        type Gf = crate::field::$gf::Gf;
        type Curve = crate::curve::binary::Curve<M>;
        type Affine = crate::curve::binary::Affine<M>;
        pub fn curve() -> impl Strategy<Value = Curve> {
            $curve
        }

        /// A uniformly random x with Tr(x) = 1 lands in E\[r\] about half the time.
        pub fn point(c: Curve) -> impl Strategy<Value = Affine> {
            any::<u128>().prop_filter_map("not on curve", move |v| c.decode(M::to_bytes(v)))
        }

        pub fn curve_and_points(n: usize) -> impl Strategy<Value = (Curve, Vec<Affine>)> {
            curve().prop_flat_map(move |c| (Just(c), prop::collection::vec(point(c), n)))
        }

        proptest! {
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
            fn sum_batch_matches_fold((_c, ps) in curve_and_points(17), k in 0usize..17) {
                let seq = ps[..k].iter().fold(Affine::IDENTITY, |acc, p| acc.add(p));
                prop_assert_eq!(sum_batch(&ps[..k]), seq);
            }
        }
    };
}
pub(crate) use suite;
