//! The benchmark fixture registry (benches/common/families.rs): its inventory
//! of representations, and the test it generates for its own constructors,
//! which this include brings into a test binary.
#[path = "../benches/common/mod.rs"]
mod common;

use common::families::{Fixture, Fixtures, Visitor};
use ephemeral_ecmh::group::{Decode, HashToCurve, Negate};
use std::collections::BTreeSet;

#[test]
fn registry_retains_the_recorded_group_inventory() {
    struct Inventory(BTreeSet<String>);
    impl Visitor for Inventory {
        fn visit<G: HashToCurve + Negate + Decode>(&mut self, fixture: Fixture<G>) {
            assert!(self.0.insert(fixture.family.id()), "duplicate family");
        }
    }
    // A literal rather than a file under report/tests/fixtures, which the
    // crate source of the checks leaves out (nix/package.nix).
    let recorded = [
        "binary-lambda.109",
        "binary-lambda.122",
        "binary-lambda.122-gls",
        "binary-lambda.127",
        "binary-u.109",
        "binary-u.122",
        "binary-u.122-gls",
        "binary-u.127",
        "binary-w.109",
        "binary-w.122",
        "binary-w.122-gls",
        "binary-w.127",
        "binary.109",
        "binary.122",
        "binary.122-gls",
        "binary.127",
        "edwards.107",
        "edwards.127",
        "edwards.61x2",
        "twisted.128",
        "twisted.61x2",
        "twisted.64x2",
        "twisted.goldilocks2",
        "weier-jacobian.107",
        "weier-jacobian.127",
        "weier-jacobian.61x2",
        "weier.107",
        "weier.127",
        "weier.61x2",
    ]
    .into_iter()
    .map(str::to_owned)
    .collect::<BTreeSet<_>>();
    assert!(!recorded.is_empty());
    let mut inventory = Inventory(BTreeSet::new());
    Fixtures::new().visit(&mut inventory);
    assert_eq!(inventory.0, recorded);
}

/// The ristretto255 wrapper of the RIBLT comparison rows reconciles 50
/// differences in 100 cells, under both purity-test orders.
#[test]
fn ristretto255_wrapper_peels_a_difference() {
    use common::ristretto::Ristretto255;
    use ephemeral_ecmh::riblt::{Peel, Riblt};
    let r = Riblt::new(Ristretto255, &[7; 32]);
    let (c, a, b) = (
        common::items(b"c", 100),
        common::items(b"a", 30),
        common::items(b"b", 20),
    );
    let mut cells = r.cells(100);
    r.encode(&mut cells, &[c.clone(), a.clone()].concat(), 1);
    r.encode(&mut cells, &[c, b.clone()].concat(), -1);
    for prefilter in [false, true] {
        let how = Peel {
            prefilter,
            batch: !prefilter,
        };
        let (mut x, mut y) = r.peel(&mut cells.clone(), how).expect("stalled");
        x.sort();
        y.sort();
        let (mut a, mut b) = (a.clone(), b.clone());
        a.sort();
        b.sort();
        assert_eq!((x, y), (a, b));
    }
}
