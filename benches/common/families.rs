//! Typed benchmark fixtures and their representation inventory.
//!
//! Construction verifies the first KAT certificate under its selector's
//! probable-prime policy. Every fixture is verified in Rust. Workload salts
//! remain the caller's choice. Only metadata consumed by the suites lives here.

use ephemeral_ecmh::curve::{self, binary};
use ephemeral_ecmh::curvegen::select;
use ephemeral_ecmh::group::{Decode, HashToCurve, Negate};

#[path = "../../tests/common/kats.rs"]
pub(crate) mod kats;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Representation {
    Native,
    Lambda,
    LambdaW,
    UnscaledW,
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
            use ephemeral_ecmh::group::{Encode, Group};
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
}
