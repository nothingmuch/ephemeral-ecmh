//! Quadratic selectors, codecs, signed sums, and provers against the
//! independent affine reference in `sage/kat_fp2.sage`.

#[path = "common/kats_fp2.rs"]
mod kats;

use ephemeral_ecmh::curve::{edwards61x2, twisted, weier61x2};
use ephemeral_ecmh::curvegen::criteria::OrderError;
use ephemeral_ecmh::curvegen::criteria::{self, Count};
use ephemeral_ecmh::curvegen::prove::{self, Factor, Family, NoSieve, Sieve, SmallL};
use ephemeral_ecmh::curvegen::select::{Certificate, Error};
use ephemeral_ecmh::curvegen::select_fp2 as select;
use ephemeral_ecmh::curvegen::sieve;
use ephemeral_ecmh::field::Packed;
use ephemeral_ecmh::group::{Decode, Encode, Group, HashToCurve, Negate, SumBatch};
use ephemeral_ecmh::hash::Salted;
use kats::CertKat;
use std::collections::HashMap;

trait Parameter {
    fn parameter(&self) -> u128;
}
impl<F: Packed> Parameter for twisted::Curve<F> {
    fn parameter(&self) -> u128 {
        self.d.pack()
    }
}
impl Parameter for edwards61x2::Curve {
    fn parameter(&self) -> u128 {
        self.d.pack()
    }
}
impl Parameter for weier61x2::Curve {
    fn parameter(&self) -> u128 {
        self.b.pack()
    }
}

/// What the checks need of a candidate beyond its parameter: the encodings
/// of tampered witnesses, on candidates of any order, and the group of an
/// accepted one, which a Weierstrass curve is only at odd order.
trait Candidate: Parameter + Copy {
    type Group: HashToCurve + Negate + SumBatch + Decode<Encoding = [u8; 16]>;
    fn group(self) -> Self::Group;
    fn identity_encoding(&self) -> [u8; 16];
    fn negate_encoding(&self, e: &[u8; 16]) -> Option<[u8; 16]>;
}
macro_rules! group_candidate {
    ($($curve:ty),*) => {$(
        impl Candidate for $curve {
            type Group = Self;
            fn group(self) -> Self {
                self
            }
            fn identity_encoding(&self) -> [u8; 16] {
                Encode::encode(self, &Group::to_affine(self, &Group::identity(self)))
            }
            fn negate_encoding(&self, e: &[u8; 16]) -> Option<[u8; 16]> {
                Some(Encode::encode(self, &Negate::neg(self, &Decode::decode(self, e)?)))
            }
        }
    )*};
}
group_candidate!(
    twisted::Curve<ephemeral_ecmh::field::fp61x2::Fq>,
    twisted::Curve<ephemeral_ecmh::field::fp64x2::Fq>,
    twisted::Curve<ephemeral_ecmh::field::goldilocks2::Fq>,
    edwards61x2::Curve
);
impl Candidate for weier61x2::Curve {
    type Group = weier61x2::OddCurve;
    fn group(self) -> weier61x2::OddCurve {
        weier61x2::OddCurve::new(self).expect("fixture curve of even order")
    }
    fn identity_encoding(&self) -> [u8; 16] {
        weier61x2::Affine::IDENTITY.encode()
    }
    fn negate_encoding(&self, e: &[u8; 16]) -> Option<[u8; 16]> {
        Some(self.decode(*e)?.neg().encode())
    }
}

struct Table {
    orders: HashMap<u128, u128>,
    visited: Vec<u128>,
}
impl Table {
    fn new(k: &CertKat) -> Self {
        Self {
            orders: k.orders.iter().copied().collect(),
            visited: vec![],
        }
    }
}
impl<C: Parameter> Count<C> for Table {
    fn order(&mut self, c: &C) -> u128 {
        let parameter = c.parameter();
        self.visited.push(parameter);
        self.orders[&parameter]
    }
    fn order_early_abort(&mut self, c: &C, tors: u32) -> Option<u128> {
        let n = self.order(c);
        (![2u128, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]
            .into_iter()
            .any(|l| n.is_multiple_of(l) && !u128::from(tors).is_multiple_of(l)))
        .then_some(n)
    }
}
struct Factors(HashMap<u128, u128>);
impl Factor for Factors {
    fn smallest_prime_factor(&mut self, n: u128) -> Option<u128> {
        Some(self.0[&n])
    }
}

fn certificate(k: &CertKat) -> Certificate {
    Certificate {
        index: k.index,
        r: k.r,
        rejections: k
            .rejections
            .iter()
            .map(|&(l, p)| (l, p.to_le_bytes()))
            .collect(),
    }
}

fn vectors<G>(g: G, k: &CertKat)
where
    G: HashToCurve + Negate + SumBatch + Decode<Encoding = [u8; 16]>,
{
    let h = Salted::new(b"kat", &k.seed);
    let ps: Vec<_> = (0..8).map(|i| g.hash(&h, &[i])).collect();
    for (p, expected) in ps.iter().zip(k.hashes) {
        let bytes = expected.to_le_bytes();
        assert_eq!(g.encode(p), bytes);
        assert_eq!(g.encode(&g.decode(&bytes).unwrap()), bytes);
    }
    assert_eq!(g.encode(&g.sum_batch(&ps)), k.sum.to_le_bytes());
    let (mut sum, mut signed, mut prepared) = (g.identity(), g.identity(), g.identity());
    for (i, p) in ps.iter().enumerate() {
        sum = g.add_affine(&sum, p);
        let a = g.prepare(p);
        if i % 2 == 0 {
            signed = g.add_affine(&signed, p);
            prepared = g.add(&prepared, &a);
        } else {
            signed = g.add_affine(&signed, &g.neg(p));
            prepared = g.add(&prepared, &g.neg_addend(&a));
        }
    }
    assert_eq!(g.encode(&g.to_affine(&sum)), k.sum.to_le_bytes());
    for p in [signed, prepared] {
        assert_eq!(g.encode(&g.to_affine(&p)), k.signed_sum.to_le_bytes());
    }
}

/// The odd l the sieve checks are run with.
const L_MAX: u32 = 11;

type Torsion<C> = fn(&C, u32) -> Option<[u8; 16]>;

/// The sieve against the KAT's orders: a point of order l (or 8, or 2,
/// the family's `first` rule) exactly when l | #E, on the first
/// candidates; and over the whole KAT, the sieve's rejections are the
/// KAT's labels, every one up to L_MAX, and verify, and the selection is
/// the KAT's, counting the rest only.
fn check_sieve<F: Family>(
    k: &CertKat,
    first: u32,
    torsion: Torsion<F::Curve>,
    verify: fn(&[u8; 32], &Certificate) -> Result<F::Curve, Error>,
) where
    F::Curve: Parameter,
    SmallL: Sieve<F::Curve>,
{
    let orders: HashMap<u128, u128> = k.orders.iter().copied().collect();
    for c in (0..=k.index).filter_map(|j| F::candidate(&k.seed, j)) {
        assert_eq!(F::reject_without_count(&c), F::quick_reject(&c).is_some());
    }
    let live = (0..k.index).filter_map(|j| F::candidate(&k.seed, j));
    for c in live.take(32) {
        let n = orders[&c.parameter()];
        for l in [first, 3, 5, 7, 11] {
            assert_eq!(
                torsion(&c, l).is_some(),
                n.is_multiple_of(l.into()),
                "l = {l}"
            );
        }
    }
    let mut counts = Table::new(k);
    let factors = &mut Factors(k.factors.iter().copied().collect());
    let (sieved, _) = prove::prove::<F>(&k.seed, &mut counts, factors, &mut SmallL(L_MAX));
    assert_eq!((sieved.index, sieved.r), (k.index, k.r));
    let ls = |c: &[(u128, [u8; 16])]| c.iter().map(|e| e.0).collect::<Vec<_>>();
    let kat: Vec<u128> = k.rejections.iter().map(|e| e.0).collect();
    assert_eq!(ls(&sieved.rejections), kat);
    let large = kat
        .iter()
        .filter(|&&l| l != first.into() && l > L_MAX.into())
        .count();
    assert_eq!(counts.visited.len(), large + 1);
    verify(&k.seed, &sieved).unwrap();
    let found = criteria::find::<F>(&k.seed, &mut Table::new(k), &mut SmallL(L_MAX));
    assert_eq!((found.index, found.r), (k.index, k.r));
}

fn check<F: Family>(
    fixtures: &[CertKat],
    verify: fn(&[u8; 32], &Certificate) -> Result<F::Curve, Error>,
    first: u32,
    torsion: Torsion<F::Curve>,
) where
    F::Curve: Candidate,
    SmallL: Sieve<F::Curve>,
{
    for k in fixtures {
        let cert = certificate(k);
        let g = verify(&k.seed, &cert).unwrap();
        vectors(g.group(), k);
        let mut counts = Table::new(k);
        let (generated, curve) = prove::prove::<F>(
            &k.seed,
            &mut counts,
            &mut Factors(k.factors.iter().copied().collect()),
            &mut NoSieve,
        );
        assert_eq!((generated.index, generated.r), (k.index, k.r));
        // Both implementations use the specified deterministic witness hash;
        // their independent group arithmetic must give the same encodings.
        assert_eq!(generated.rejections, cert.rejections);
        assert_eq!(curve.parameter(), g.parameter());
        assert_eq!(
            counts.visited,
            k.orders.iter().map(|row| row.0).collect::<Vec<_>>()
        );
        verify(&k.seed, &generated).unwrap();
        let found = criteria::find::<F>(&k.seed, &mut Table::new(k), &mut NoSieve);
        assert_eq!(
            (found.index, found.r, found.curve.parameter()),
            (k.index, k.r, g.parameter())
        );

        // Identity never witnesses rejection, in either the full or quotient group.
        let mut changed = cert.clone();
        let rejected = (0..k.index)
            .find_map(|j| F::candidate(&k.seed, j).map(|c| (j, c)))
            .unwrap();
        changed.rejections[0].1 = rejected.1.identity_encoding();
        assert!(
            matches!(verify(&k.seed, &changed), Err(Error::BadRejection(j)) if j == rejected.0)
        );
        let mut changed = cert.clone();
        changed.rejections.pop();
        assert!(matches!(
            verify(&k.seed, &changed),
            Err(Error::RejectionCount)
        ));
        let mut changed = cert.clone();
        changed.rejections.push((3, [0; 16]));
        assert!(matches!(
            verify(&k.seed, &changed),
            Err(Error::RejectionCount)
        ));
        let mut changed = cert.clone();
        changed.r = 15;
        assert!(matches!(
            verify(&k.seed, &changed),
            Err(Error::Order(OrderError::RNotPrime))
        ));

        // Witnesses are not unique: replacing every P by -P is valid.
        let mut negated = cert.clone();
        let mut entries = negated.rejections.iter_mut();
        for j in 0..k.index {
            if let Some(c) = F::candidate(&k.seed, j) {
                let (_, bytes) = entries.next().unwrap();
                *bytes = c.negate_encoding(bytes).unwrap();
            }
        }
        verify(&k.seed, &negated).unwrap();
        assert!(std::panic::catch_unwind(|| F::verdict(0)).is_err());

        check_sieve::<F>(k, first, torsion, verify);
    }
}

macro_rules! family_test {
    ($test:ident, $family:ident, $fixtures:ident, $verify:ident, $first:expr, $torsion:expr) => {
        #[test]
        fn $test() {
            check::<prove::$family>(kats::$fixtures, select::$verify, $first, $torsion);
        }
    };
}
family_test!(
    twisted61x2,
    Twisted61x2,
    TWISTED61X2_CERTS,
    verify_twisted61x2,
    8,
    sieve::twisted_torsion
);
family_test!(
    twisted64x2,
    Twisted64x2,
    TWISTED64X2_CERTS,
    verify_twisted64x2,
    8,
    sieve::twisted_torsion
);
family_test!(
    twisted_goldilocks2,
    TwistedGoldilocks2,
    TWISTED_GOLDILOCKS2_CERTS,
    verify_twisted_goldilocks2,
    8,
    sieve::twisted_torsion
);
family_test!(
    edwards61x2,
    Edwards61x2,
    EDWARDS61X2_CERTS,
    verify_edwards61x2,
    8,
    sieve::edwards_torsion
);
family_test!(
    weier61x2,
    Weier61x2,
    WEIER61X2_CERTS,
    verify_weier61x2,
    2,
    sieve::weier_torsion
);
