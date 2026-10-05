//! Numerical vectors from an independent integer implementation of RFC 9380
//! §6.6.2, cross-checked against Appendix F.2. The curve parameters reproduce
//! the first certified Weier107/127 curves in the Sage certificate fixtures.

use super::*;
use crate::field::{fp107, fp127};

#[test]
fn reference_vectors_107() {
    type Fp = fp107::Fp;
    let curve = Curve::new(Fp::new(0x6beb73f4f8778caf14d57fb6175)).unwrap();
    let map = Sswu::new(curve, Fp::new(0x7fffffffffffffffffffffffff8)).unwrap();
    let vectors = [
        (
            0x0,
            0x72da733c6d5b06b300fdf9e4197,
            0x762beca6971bf1938c0b4f16260,
        ),
        (
            0x7ef1b30e4343858f5454fb40b4a,
            0x72da733c6d5b06b300fdf9e4197,
            0x762beca6971bf1938c0b4f16260,
        ),
        (
            0x10e4cf1bcbc7a70abab04bf4b5,
            0x72da733c6d5b06b300fdf9e4197,
            0x9d4135968e40e6c73f4b0e9d9f,
        ),
        (
            0x1,
            0x21c839dbba61b002dc717ee2b6b,
            0x734b13a82ec324ced99f483da89,
        ),
        (
            0x2,
            0x22613614e796161a2a545021583,
            0x5199050c77f390ffcc015658d34,
        ),
        (
            0x3,
            0x16a7571ebb230da7dae3fa5d8e5,
            0x1bd4f47e80f14d645dd9516521b,
        ),
        (
            0x4,
            0x28265bc8630e89374b78068d0e0,
            0x5d87e6abc1263585a66a291b4a6,
        ),
        (
            0x5,
            0x2297139fd947712b421da8a08e9,
            0x766332aeb47f3bd8284897fb19b,
        ),
        (
            0x11,
            0x1710854cea9a33aaaf4713214e9,
            0x4493a412692be93d19568e78fe1,
        ),
        (
            0x7fffffffffffffffffffffffffe,
            0x21c839dbba61b002dc717ee2b6b,
            0xcb4ec57d13cdb312660b7c2576,
        ),
        (
            0x7fffffffffffffffffffffffffd,
            0x22613614e796161a2a545021583,
            0x2e66faf3880c6f0033fea9a72cb,
        ),
    ];
    for (u, x, y) in vectors {
        let p = map.map_to_curve(Fp::new(u));
        assert_eq!(
            p,
            Affine {
                x: Fp::new(x),
                y: Fp::new(y)
            }
        );
        assert!(curve.is_on_curve(&p));
    }
}

#[test]
fn rejects_invalid_z_107() {
    type Fp = fp107::Fp;
    let curve = Curve::new(Fp::new(0x6beb73f4f8778caf14d57fb6175)).unwrap();
    for z in [0, 1, fp107::P - 1] {
        assert!(Sswu::new(curve, Fp::new(z)).is_none());
    }
}

#[test]
fn reference_vectors_127() {
    type Fp = fp127::Fp;
    let curve = Curve::new(Fp::new(0x7fa0081246bc73f7e3c59c4fbf247061)).unwrap();
    let map = Sswu::new(curve, Fp::new(0x3)).unwrap();
    let vectors = [
        (
            0x0,
            0x71d1c636dbb22c00e6b127855c8a2c66,
            0x1e111e37b28ff74ee58f23c6357582b6,
        ),
        (
            0x6910301f653125bf71a44404739af006,
            0x71d1c636dbb22c00e6b127855c8a2c66,
            0x1e111e37b28ff74ee58f23c6357582b6,
        ),
        (
            0x16efcfe09aceda408e5bbbfb8c650ff9,
            0x71d1c636dbb22c00e6b127855c8a2c66,
            0x61eee1c84d7008b11a70dc39ca8a7d49,
        ),
        (
            0x1,
            0x4e163bcdb5fcf0fd12403f8e933eefb1,
            0x55ba3e64e31e4c8544ed6b3d26500715,
        ),
        (
            0x2,
            0x292855015e663bdf59d9de479871deab,
            0x3571fd457009a8c4afd87391a1d0b75e,
        ),
        (
            0x3,
            0x7992f367c083d824a0165f2400a6dc66,
            0x6548199c6dab659dd43852299c2b0bd7,
        ),
        (
            0x4,
            0xa8b2bdd24dbafbcdfdbfc08adeaf8f5,
            0x5f7c248afaa08e3a17d40101f2e925e2,
        ),
        (
            0x5,
            0x2b65fb3944c8b6860cc545989ccce39b,
            0x76edc22d10bad5bb8879ed55385f41bf,
        ),
        (
            0x11,
            0x38b932a8274222baa30e5c3b1d0d2835,
            0x75cf1709f1dcfce30137aa1b0e1490f7,
        ),
        (
            0x7ffffffffffffffffffffffffffffffe,
            0x4e163bcdb5fcf0fd12403f8e933eefb1,
            0x2a45c19b1ce1b37abb1294c2d9aff8ea,
        ),
        (
            0x7ffffffffffffffffffffffffffffffd,
            0x292855015e663bdf59d9de479871deab,
            0x4a8e02ba8ff6573b50278c6e5e2f48a1,
        ),
    ];
    for (u, x, y) in vectors {
        let p = map.map_to_curve(Fp::new(u));
        assert_eq!(
            p,
            Affine {
                x: Fp::new(x),
                y: Fp::new(y)
            }
        );
        assert!(curve.is_on_curve(&p));
    }
}

#[test]
fn rejects_invalid_z_127() {
    type Fp = fp127::Fp;
    let curve = Curve::new(Fp::new(0x7fa0081246bc73f7e3c59c4fbf247061)).unwrap();
    for z in [0, 1, fp127::P - 1] {
        assert!(Sswu::new(curve, Fp::new(z)).is_none());
    }
}

/// The affine procedure in RFC 9380 §6.6.2, independent of the ratio-root
/// and homogeneous evaluation used by the implementation.
fn reference<F: OddField>(curve: &Curve<F>, z: F, u: F) -> Affine<F> {
    let a = -(F::ONE + F::ONE + F::ONE);
    let zu2 = z * u.square();
    let denominator = zu2.square() + zu2;
    let x1 = if denominator.is_zero() {
        curve.b * (z * a).inv()
    } else {
        -curve.b * a.inv() * (F::ONE + denominator.inv())
    };
    let (x, y) = match curve.rhs(x1).sqrt() {
        Some(y) => (x1, y),
        None => {
            let x2 = zu2 * x1;
            (x2, curve.rhs(x2).sqrt().unwrap())
        }
    };
    Affine {
        x,
        y: if y.sgn0() == u.sgn0() { y } else { -y },
    }
}

macro_rules! suite {
    ($module:ident, $fp:ident, $b:expr, $z:expr, $reducible:expr, $bad_exception:expr) => {
        mod $module {
            use super::*;
            use crate::field::batch::Invert;
            use crate::field::$fp::{Fp, P, tests::fp};
            use proptest::prelude::*;

            #[test]
            fn setup_requires_each_rfc_condition() {
                let curve = Curve::new(Fp::new($b)).unwrap();
                let a = -Fp::new(3);
                let z = Fp::new($reducible);
                assert!(!z.is_square() && z != -Fp::ONE);
                assert!(curve.rhs(curve.b * (z * a).inv()).is_square());
                assert!(!cubic_irreducible(a, curve.b - z));
                assert!(Sswu::new(curve, z).is_none());

                let z = Fp::new($bad_exception);
                assert!(!z.is_square() && z != -Fp::ONE);
                assert!(cubic_irreducible(a, curve.b - z));
                assert!(!curve.rhs(curve.b * (z * a).inv()).is_square());
                assert!(Sswu::new(curve, z).is_none());

                for b in [Fp::ZERO, Fp::new(2), -Fp::new(2)] {
                    assert!(Sswu::new(Curve { b }, Fp::new($z)).is_none());
                }
                assert!(!cubic_irreducible(Fp::ZERO, Fp::ZERO));
                assert!(!cubic_irreducible(Fp::ZERO, -Fp::ONE));
            }

            #[test]
            fn ratio_zero_is_square() {
                let cache = Fp::prepare_ratio(Fp::new($z));
                assert_eq!(Fp::ratio_root(Fp::ZERO, Fp::ONE, &cache), (true, Fp::ZERO));
            }

            #[test]
            fn search_respects_bound_and_rfc_order() {
                let curve = Curve::new(Fp::new($b)).unwrap();
                // The independent vector generator searched 1, -1, 2, -2, ... .
                let z = Fp::new($z);
                let bound = z.value().min(P - z.value()) as u32;
                assert!(Sswu::search(curve, 0).is_none());
                assert!(Sswu::search(curve, bound - 1).is_none());
                assert_eq!(Sswu::search(curve, bound).unwrap().z, z);
                assert_eq!(Sswu::search(curve, 256).unwrap().z, z);
                for b in [Fp::ZERO, Fp::new(2), -Fp::new(2)] {
                    assert!(Sswu::search(Curve { b }, 8).is_none());
                }
            }

            #[test]
            fn digest_adapter_reduces_low_field_bits() {
                use crate::curve::h2c::Map;

                let curve = Curve::new(Fp::new($b)).unwrap();
                let z = Fp::new($z);
                let map = Sswu::new(curve, z).unwrap();
                let bits = 128 - P.leading_zeros();
                for c in [
                    0,
                    1,
                    P - 1,
                    P,
                    P + 1,
                    P + 2,
                    (P - 1) | (1 << bits),
                    u128::MAX - 1,
                    u128::MAX,
                ] {
                    // Both supported moduli are 2^bits - 1.
                    let u = Fp::new((c & P) % P);
                    assert_eq!(map.map(c), reference(&curve, z, u));
                }
            }

            #[test]
            fn map1_uses_counter_zero_low_half() {
                use crate::curve::h2c::map1;
                use crate::hash::Salted;

                let curve = Curve::new(Fp::new($b)).unwrap();
                let z = Fp::new($z);
                let map = Sswu::new(curve, z).unwrap();
                let salt = Salted::new(b"sswu-test", &[0x5a; 32]);
                for msg in [b"".as_slice(), b"a", b"ECMH point identifier"] {
                    let digest = salt.digest(msg, 0);
                    let c = u128::from_le_bytes(digest[..16].try_into().unwrap());
                    let u = Fp::new((c & P) % P);
                    assert_eq!(map1(&map, &salt, msg), reference(&curve, z, u));
                }
            }

            proptest! {
                #[test]
                fn matches_affine_reference(b in fp(), u in fp()) {
                    let Some(curve) = Curve::new(b) else { return Ok(()); };
                    let candidate = (1..=64).flat_map(|v| [Fp::new(v), -Fp::new(v)])
                        .find_map(|z| Sswu::new(curve, z).map(|map| (z, map)));
                    prop_assume!(candidate.is_some());
                    let (z, map) = candidate.unwrap();
                    let p = map.map_to_curve(u);
                    prop_assert_eq!(p, reference(&curve, z, u));
                    prop_assert!(curve.is_on_curve(&p));
                    if !u.is_zero() {
                        prop_assert_eq!(map.map_to_curve(-u), p.neg());
                    }
                }

                #[test]
                fn ratio_root_handles_both_square_classes(x in fp(), d in fp()) {
                    prop_assume!(!d.is_zero());
                    let z = Fp::new($z);
                    let cache = Fp::prepare_ratio(z);
                    let n = d * x.square();
                    let (square, r) = Fp::ratio_root(n, d, &cache);
                    prop_assert!(square);
                    prop_assert_eq!(r.square() * d, n);
                    if !x.is_zero() {
                        let n = z * n;
                        let (square, r) = Fp::ratio_root(n, d, &cache);
                        prop_assert!(!square);
                        prop_assert_eq!(r.square() * d, z * n);
                    }
                }

                #[test]
                fn reducible_cubics_are_rejected(a in fp(), root in fp()) {
                    let b = -root * (root.square() + a);
                    prop_assert!(!cubic_irreducible(a, b));
                }
            }
        }
    };
}

suite!(small, fp107, 0x6beb73f4f8778caf14d57fb6175, P - 7, 3, P - 8);
suite!(
    large,
    fp127,
    0x7fa0081246bc73f7e3c59c4fbf247061,
    3,
    P - 4,
    P - 2
);

mod fp61x2 {
    use super::*;
    use crate::field::fp61x2::{Fp, Fq, tests::fq};
    use proptest::prelude::*;

    /// RFC 9380's g for GF(p^2) = F_p(i).
    const I: Fq = Fq::new(Fp::ZERO, Fp::ONE);

    fn map(b: Fq) -> Option<(Fq, Sswu<Fq>)> {
        let curve = Curve::new(b)?;
        let map = Sswu::search_from(curve, I, 64)?;
        Some((map.z, map))
    }

    #[test]
    fn ratio_zero_is_square() {
        assert_eq!(
            Fq::ratio_root(Fq::ZERO, Fq::ONE, &Fq::NON_SQUARE),
            (true, Fq::ZERO)
        );
    }

    /// The exceptional denominator Z^2 u^4 + Z u^2 = 0 at u = 0 only: -1
    /// is a square in GF(p^2), so u^2 = -1/Z has no root for a non-square Z.
    #[test]
    fn exceptional_input_matches_affine_reference() {
        let b = Fq::new(Fp::new(5), Fp::new(7));
        let (z, map) = map(b).expect("fixed curve admits SSWU setup");
        let curve = Curve::new(b).unwrap();
        assert!(!(-z.invert()).is_square());
        let p = map.map_to_curve(Fq::ZERO);
        assert_eq!(p, reference(&curve, z, Fq::ZERO));
        assert!(curve.is_on_curve(&p));
    }

    proptest! {
        #[test]
        fn matches_affine_reference(b in fq(), u in fq()) {
            let candidate = map(b);
            prop_assume!(candidate.is_some());
            let (z, map) = candidate.unwrap();
            let curve = Curve::new(b).unwrap();
            let p = map.map_to_curve(u);
            prop_assert_eq!(p, reference(&curve, z, u));
            prop_assert!(curve.is_on_curve(&p));
            if !u.is_zero() {
                prop_assert_eq!(map.map_to_curve(-u), p.neg());
            }
        }

        #[test]
        fn ratio_root_handles_both_square_classes(x in fq(), d in fq()) {
            prop_assume!(!d.is_zero());
            let z = Fq::NON_SQUARE;
            let n = d * x.square();
            let (square, r) = Fq::ratio_root(n, d, &z);
            prop_assert!(square);
            prop_assert_eq!(r.square() * d, n);
            if !x.is_zero() {
                let n = z * n;
                let (square, r) = Fq::ratio_root(n, d, &z);
                prop_assert!(!square);
                prop_assert_eq!(r.square() * d, z * n);
            }
        }
    }
}
