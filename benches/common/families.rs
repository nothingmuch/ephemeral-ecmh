//! Typed benchmark fixtures and their representation inventory.
//!
//! Construction verifies the first KAT certificate under its selector's
//! probable-prime policy. Every fixture is verified in Rust. Workload salts
//! remain the caller's choice. Only metadata consumed by the suites lives here.

use ephemeral_ecmh::curve::{self, binary, weier};
use ephemeral_ecmh::curvegen::{select, select_fp2, select107, select109, select122, select128};
use ephemeral_ecmh::field::OddField;
use ephemeral_ecmh::group::{Decode, HashToCurve, Negate};

#[path = "../../tests/common/kats.rs"]
pub(crate) mod kats;
#[path = "../../tests/common/kats107.rs"]
pub(crate) mod kats107;
#[path = "../../tests/common/kats109.rs"]
pub(crate) mod kats109;
#[path = "../../tests/common/kats122.rs"]
pub(crate) mod kats122;
#[path = "../../tests/common/kats128.rs"]
pub(crate) mod kats128;
#[path = "../../tests/common/kats_fp2.rs"]
pub(crate) mod kats_fp2;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Representation {
    Native,
    Lambda,
    LambdaW,
    UnscaledW,
    Jacobian,
}

#[derive(Clone, Copy, Debug)]
pub struct CurveInfo {
    pub cofactor: u32,
    pub automorphisms: u32,
    // The benchmark id's spelling, which need not name the mathematical model.
    pub id_model: &'static str,
    pub id_field: &'static str,
}

#[derive(Clone, Copy, Debug)]
pub struct Family {
    pub curve: CurveInfo,
    pub representation: Representation,
}

impl Family {
    pub fn id(self) -> String {
        let suffix = match self.representation {
            Representation::Native => "",
            Representation::Lambda => "-lambda",
            Representation::LambdaW => "-w",
            Representation::UnscaledW => "-u",
            Representation::Jacobian => "-jacobian",
        };
        format!("{}{suffix}.{}", self.curve.id_model, self.curve.id_field)
    }
}

#[derive(Clone, Copy)]
pub struct Fixture<G> {
    pub group: G,
    pub seed: [u8; 32],
    pub index: u32,
    pub r: u128,
    pub family: Family,
}

impl<G> Fixture<G> {
    fn represented<H>(self, group: H, representation: Representation) -> Fixture<H> {
        Fixture {
            group,
            seed: self.seed,
            index: self.index,
            r: self.r,
            family: Family {
                representation,
                ..self.family
            },
        }
    }
}

/// A statically dispatched consumer; the timed operations keep their concrete types.
pub trait Visitor {
    fn visit<G: HashToCurve + Negate + Decode>(&mut self, fixture: Fixture<G>);
}

macro_rules! variants {
    (single, $v:ident, $f:expr) => {
        $v.visit($f);
    };
    (binary, $v:ident, $f:expr) => {{
        let f = $f;
        $v.visit(f);
        $v.visit(f.represented(binary::lambda::Curve(f.group), Representation::Lambda));
        $v.visit(f.represented(binary::wcodec::Curve::new(f.group), Representation::LambdaW));
        $v.visit(f.represented(
            binary::unscaled::Curve::new(f.group),
            Representation::UnscaledW,
        ));
    }};
    (weier, $v:ident, $f:expr) => {{
        let f = $f;
        $v.visit(f);
        $v.visit(f.represented(weier::jacobian::Curve(f.group), Representation::Jacobian));
    }};
}

/// Verification returns the raw curve. The group is its odd-order form,
/// which `OddCurve::new` checks rather than infers from the certificate,
/// whose order argument is conditional on r's primality. Every Weierstrass
/// policy has cofactor 1, so for the fixtures #E = r is an odd prime and
/// the `expect` is a consistency check, not a path.
fn odd<const N: usize, F: OddField>(
    verify: fn(&[u8; 32], &select::Certificate<N>) -> Result<weier::Curve<F>, select::Error>,
) -> impl Fn(&[u8; 32], &select::Certificate<N>) -> Result<weier::OddCurve<F>, select::Error> {
    move |seed, cert| {
        verify(seed, cert).map(|c| weier::OddCurve::new(c).expect("fixture curve of even order"))
    }
}

// The same declaration creates direct constructors for the comparison suites
// and the complete inventory for group operations and RIBLT workloads.
macro_rules! registry {
    ($( $name:ident: $ty:ty, $kat:path, $verify:expr,
        ($id_model:literal, $id_field:literal, $cofactor:literal, $automorphisms:literal), $variants:ident; )*) => {
        $(pub fn $name() -> Fixture<$ty> {
            let k = &$kat[0];
            let cert = super::certificate(k.index, k.r, k.rejections);
            Fixture {
                group: ($verify)(&k.seed, &cert).unwrap(),
                seed: k.seed,
                index: k.index,
                r: k.r,
                family: Family {
                    curve: CurveInfo {
                        cofactor: $cofactor, automorphisms: $automorphisms,
                        id_model: $id_model, id_field: $id_field,
                    },
                    representation: Representation::Native,
                },
            }
        })*

        // Generated here, where the constructor/KAT pairs are declared once;
        // it runs through tests/bench_families.rs, which includes this module.
        #[cfg(test)]
        #[test]
        fn base_fixtures_reproduce_kat_hashes_and_sums() {
            use ephemeral_ecmh::group::{Accumulate, Encode, Group};
            use ephemeral_ecmh::hash::Salted;
            $(
                let f = $name();
                let k = &$kat[0];
                let salt = Salted::new(b"kat", &f.seed);
                let mut sum = f.group.identity();
                for (i, expected) in k.hashes.iter().enumerate() {
                    let p = f.group.hash(&salt, &[i as u8]);
                    let e = f.group.encode(&p);
                    assert_eq!(e.as_ref(), &expected.to_le_bytes()[..e.as_ref().len()], stringify!($name));
                    sum = f.group.add_affine(&sum, &p);
                }
                let e = f.group.encode_point(&sum);
                assert_eq!(e.as_ref(), &k.sum.to_le_bytes()[..e.as_ref().len()], stringify!($name));
            )*
        }

        /// Verified once per suite, then reused for every measurement stage.
        pub struct Fixtures { $( pub $name: Fixture<$ty>, )* }
        impl Fixtures {
            pub fn new() -> Self { Self { $( $name: $name(), )* } }
            pub fn visit(&self, v: &mut impl Visitor) {
                $(variants!($variants, v, self.$name);)*
            }
        }
    };
}

registry! {
    binary127: curve::binary127::Curve, kats::GF2_127_CERTS, select::verify_gf2_127,
        ("binary", "127", 2, 2), binary;
    edwards127: curve::edwards127::Curve, kats::FP127_CERTS, select::verify_fp127,
        ("edwards", "127", 4, 2), single;
    twisted128: curve::twisted128::Curve, kats128::FP128_CERTS, select128::verify_fp128,
        ("twisted", "128", 4, 2), single;
    weier127: curve::weier127::OddCurve, kats::WEIER127_CERTS, odd(select::verify_weier127),
        ("weier", "127", 1, 2), weier;
    binary109: curve::binary109::Curve, kats109::GF2_109_CERTS, select109::verify,
        ("binary", "109", 2, 2), binary;
    binary122: curve::binary122::Dense, kats122::GF2_122_CERTS, select122::verify_dense,
        ("binary", "122", 2, 2), binary;
    binary122_gls: curve::binary122::Gls, kats122::GF2_122_GLS_CERTS, select122::verify_gls,
        ("binary", "122-gls", 2, 4), binary;
    edwards107: curve::edwards107::Curve, kats107::FP107_CERTS, select107::verify_fp107,
        ("edwards", "107", 4, 2), single;
    weier107: curve::weier107::OddCurve, kats107::WEIER107_CERTS, odd(select107::verify_weier107),
        ("weier", "107", 1, 2), weier;
    edwards61x2: curve::edwards61x2::Curve, kats_fp2::EDWARDS61X2_CERTS, select_fp2::verify_edwards61x2,
        ("edwards", "61x2", 4, 2), single;
    weier61x2: curve::weier61x2::OddCurve, kats_fp2::WEIER61X2_CERTS, odd(select_fp2::verify_weier61x2),
        ("weier", "61x2", 1, 2), weier;
    twisted61x2: curve::twisted61x2::Curve, kats_fp2::TWISTED61X2_CERTS, select_fp2::verify_twisted61x2,
        ("twisted", "61x2", 4, 2), single;
    twisted64x2: curve::twisted64x2::Curve, kats_fp2::TWISTED64X2_CERTS, select_fp2::verify_twisted64x2,
        ("twisted", "64x2", 4, 2), single;
    twisted_goldilocks2: curve::twisted_goldilocks2::Curve, kats_fp2::TWISTED_GOLDILOCKS2_CERTS, select_fp2::verify_twisted_goldilocks2,
        ("twisted", "goldilocks2", 4, 2), single;
}
