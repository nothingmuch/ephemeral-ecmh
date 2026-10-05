//! The hash constructions a RIBLT row can run under, and the adapter that
//! runs a group under one of them.
//!
//! Each curve's `HashToCurve` is its default construction, which the
//! group suite measures as `h2c`. The comparison suites measure the
//! alternatives the library provides: Pornin's map on the binary curves,
//! straight to the point or to a λ-affine or unscaled addend, and
//! Elligator 2 and SSWU on the prime-field curves. `Constructions` lists,
//! per group, every construction that yields the group's addend, under
//! the tokens the report's chooser reads and writes:
//!
//! - `ti`: the group's own hash, then `prepare`;
//! - `ti-addend`: try-and-increment straight to the addend
//!   (`hash_to_lambda`, `hash_to_edwards`);
//! - `pornin`: Pornin's map to the extended point;
//! - `pornin-addend`: Pornin's map to the λ-affine or unscaled addend;
//! - `elligator2`, `sswu`: one map to the affine point, then `prepare`.
//!
//! The w codec's accumulator is λ-projective with λ-affine addends, the
//! type `map_to_lambda` produces, so its `pornin-addend` is the λ
//! family's. Batched variants share one inversion where the library has
//! one; Elligator 2 and SSWU have none, so their batch maps each item.
//!
//! `RIBLT_PLAN` names a TSV of `<family>\t<token>\t<scope>\t<proj>\t<mapproj>`
//! lines, the family being the RIBLT id's function component
//! (`binary-u.127`), the scope `full`, `buffer` or `spot`: how much of each
//! swept dimension the family measures under that construction, and the
//! projection fields (`PROJECTIONS`) its ID-keyed rows use for the checksum
//! and for the mapping, or `-` in both for none. Several lines per family
//! are allowed. With a plan set, a family it omits is spot-checked under
//! every construction it has, with no ID-keyed rows; unset, every family
//! runs every construction in full and under every projection, the
//! mapping's in the checksum's field, for exploration and for `--list`.

use std::collections::HashMap;

use ephemeral_ecmh::curve::binary::{self, Pornin};
use ephemeral_ecmh::curve::encoding::Signed;
use ephemeral_ecmh::curve::h2c::{Elligator2, Map, map1};
use ephemeral_ecmh::curve::twisted::MontgomeryModel;
use ephemeral_ecmh::curve::weier::{Sswu, SswuField};
use ephemeral_ecmh::curve::{edwards, twisted, weier};
use ephemeral_ecmh::field::{OddField, Packed, fp61x2, fp107, fp127};
use ephemeral_ecmh::group::{Accumulate, Group, HashToCurve, Negate};
use ephemeral_ecmh::hash::{Field, Salted};

/// The projection fields by the token ids spell them in, `proj=<token>`
/// and `mapproj=<token>`.
pub const PROJECTIONS: [(&str, Field); 3] = [
    ("fp130", Field::Fp130),
    ("fp127", Field::Fp127),
    ("gf2_127", Field::Gf2_127),
];

type Single<G> = Box<dyn Fn(&G, &Salted, &[u8]) -> <G as Accumulate>::Addend>;
type Batch<G> = Box<dyn Fn(&G, &Salted, &[&[u8]]) -> Vec<<G as Accumulate>::Addend>>;

/// One way from an item to a group's addend.
pub struct Construction<G: Accumulate> {
    pub token: &'static str,
    pub single: Single<G>,
    pub batch: Batch<G>,
}

impl<G: Accumulate> Construction<G> {
    fn new(token: &'static str, single: Single<G>, batch: Batch<G>) -> Self {
        Self {
            token,
            single,
            batch,
        }
    }

    /// A map with no batched form: the batch maps each item.
    fn each(
        token: &'static str,
        f: impl Fn(&G, &Salted, &[u8]) -> G::Addend + Clone + 'static,
    ) -> Self {
        let g = f.clone();
        Self::new(
            token,
            Box::new(f),
            Box::new(move |c, h, ms| ms.iter().map(|m| g(c, h, m)).collect()),
        )
    }
}

/// The group's own hash, then `prepare`.
fn ti<G: Accumulate + HashToCurve>() -> Construction<G> {
    Construction::new(
        "ti",
        Box::new(|g, h, m| g.prepare(&g.hash(h, m))),
        Box::new(|g, h, ms| g.prepare_batch(&g.hash_batch(h, ms))),
    )
}

/// The digest half every map takes, as `map1` does.
fn half(h: &Salted, m: &[u8]) -> u128 {
    h.half(m, 0, 0)
}

fn halves_of(h: &Salted, ms: &[&[u8]]) -> Vec<u128> {
    ms.iter().map(|m| half(h, m)).collect()
}

pub trait Constructions: Accumulate + HashToCurve {
    fn constructions(&self) -> Vec<Construction<Self>>;
}

impl<M: Pornin> Constructions for binary::Curve<M> {
    fn constructions(&self) -> Vec<Construction<Self>> {
        vec![
            ti(),
            Construction::new(
                "pornin",
                Box::new(map1),
                Box::new(|g, h, ms| g.hash_to_curve_map1_batch(h, ms)),
            ),
        ]
    }
}

impl<M: Pornin> Constructions for binary::lambda::Curve<M> {
    fn constructions(&self) -> Vec<Construction<Self>> {
        vec![
            ti(),
            Construction::new(
                "ti-addend",
                Box::new(|g, h, m| g.0.hash_to_lambda(h, m)),
                Box::new(|g, h, ms| g.0.hash_to_lambda_batch(h, ms)),
            ),
            Construction::new(
                "pornin-addend",
                Box::new(|g, h, m| g.0.map_to_lambda(half(h, m))),
                Box::new(|g, h, ms| g.0.map_to_lambda_batch(&halves_of(h, ms))),
            ),
        ]
    }
}

impl<M: Pornin> Constructions for binary::wcodec::Curve<M> {
    fn constructions(&self) -> Vec<Construction<Self>> {
        vec![
            ti(),
            Construction::new(
                "pornin-addend",
                Box::new(|g, h, m| g.c.map_to_lambda(half(h, m))),
                Box::new(|g, h, ms| g.c.map_to_lambda_batch(&halves_of(h, ms))),
            ),
        ]
    }
}

impl<M: Pornin> Constructions for binary::unscaled::Curve<M> {
    fn constructions(&self) -> Vec<Construction<Self>> {
        vec![
            ti(),
            Construction::new(
                "pornin-addend",
                Box::new(|g, h, m| g.map_to_addend(half(h, m))),
                Box::new(|g, h, ms| g.map_to_addend_batch(&halves_of(h, ms))),
            ),
        ]
    }
}

impl<F: OddField + Signed + 'static> Constructions for edwards::Curve<F> {
    fn constructions(&self) -> Vec<Construction<Self>> {
        let ell = Elligator2::new(*self).expect("fixture has a2 != 0");
        vec![
            ti(),
            Construction::each("ti-addend", |g: &Self, h, m| g.hash_to_edwards(h, m)),
            Construction::each("elligator2", move |g: &Self, h, m| {
                g.prepare(&map1(&ell, h, m))
            }),
        ]
    }
}

impl<F: Packed + 'static> Constructions for twisted::Curve<F> {
    fn constructions(&self) -> Vec<Construction<Self>> {
        let ell = Elligator2::new(MontgomeryModel::new(*self)).expect("fixture has a2 != 0");
        vec![
            ti(),
            Construction::each("elligator2", move |g: &Self, h, m| {
                g.prepare(&map1(&ell, h, m))
            }),
        ]
    }
}

/// SSWU's Z search per field, as the comparison suites set it up.
pub trait SswuSetup: SswuField + Packed {
    fn sswu(curve: weier::Curve<Self>) -> Sswu<Self>;
}

impl SswuSetup for fp107::Fp {
    fn sswu(curve: weier::Curve<Self>) -> Sswu<Self> {
        Sswu::search(curve, 256).expect("fixture admits SSWU setup")
    }
}

impl SswuSetup for fp127::Fp {
    fn sswu(curve: weier::Curve<Self>) -> Sswu<Self> {
        Sswu::search(curve, 256).expect("fixture admits SSWU setup")
    }
}

impl SswuSetup for fp61x2::Fq {
    fn sswu(curve: weier::Curve<Self>) -> Sswu<Self> {
        let i = fp61x2::Fq::new(fp61x2::Fp::ZERO, fp61x2::Fp::ONE);
        Sswu::search_from(curve, i, 64).expect("fixture admits SSWU setup")
    }
}

impl<F: OddField + Signed + SswuSetup + 'static> Constructions for weier::OddCurve<F> {
    fn constructions(&self) -> Vec<Construction<Self>> {
        let sswu = F::sswu(self.curve());
        vec![
            ti(),
            Construction::each("sswu", move |_, h, m| sswu.map(half(h, m))),
        ]
    }
}

impl<F: OddField + Signed + SswuSetup + 'static> Constructions for weier::jacobian::Curve<F> {
    fn constructions(&self) -> Vec<Construction<Self>> {
        let sswu = F::sswu(self.0.curve());
        vec![
            ti(),
            Construction::each("sswu", move |_, h, m| sswu.map(half(h, m))),
        ]
    }
}

/// How much of each swept dimension a family measures, widest first.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Scope {
    Full,
    Buffer,
    Spot,
}

impl Scope {
    fn parse(s: &str) -> Option<Self> {
        match s {
            "full" => Some(Self::Full),
            "buffer" => Some(Self::Buffer),
            "spot" => Some(Self::Spot),
            _ => None,
        }
    }
}

/// The projection fields of the checksum and of the mapping.
pub type Projections = (Field, Field);

/// A family's line of `RIBLT_PLAN`: (token, scope, projections).
pub type Line = (String, Scope, Option<Projections>);

/// `RIBLT_PLAN`'s lines by family. None when unset.
pub fn plan() -> Option<HashMap<String, Vec<Line>>> {
    let path = std::env::var_os("RIBLT_PLAN")?;
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("RIBLT_PLAN={}: {e}", path.to_string_lossy()));
    let mut out: HashMap<String, Vec<Line>> = HashMap::new();
    for line in text
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
    {
        let fields: Vec<&str> = line.split('\t').collect();
        let [family, token, scope, proj, map] = fields[..] else {
            panic!("RIBLT_PLAN: expected family, token, scope and projections in {line:?}");
        };
        let scope = Scope::parse(scope).unwrap_or_else(|| panic!("RIBLT_PLAN: scope {scope:?}"));
        let field = |s: &str| {
            (s != "-").then(|| {
                PROJECTIONS
                    .iter()
                    .find(|(t, _)| *t == s)
                    .unwrap_or_else(|| panic!("RIBLT_PLAN: projection {s:?}"))
                    .1
            })
        };
        let proj = match (field(proj), field(map)) {
            (Some(p), Some(m)) => Some((p, m)),
            (None, None) => None,
            _ => panic!("RIBLT_PLAN: projections {proj:?} and {map:?}, not both or neither"),
        };
        out.entry(family.to_string())
            .or_default()
            .push((token.to_string(), scope, proj));
    }
    Some(out)
}

/// The constructions, scopes and projection fields `family` runs under
/// `plan`. A planned token the group lacks is an error, not a silent
/// omission.
pub fn planned<G: Constructions>(
    group: &G,
    family: &str,
    plan: Option<&HashMap<String, Vec<Line>>>,
) -> Vec<(Construction<G>, Scope, Vec<Projections>)> {
    let all = group.constructions();
    let Some(plan) = plan else {
        let fields = PROJECTIONS.map(|(_, f)| (f, f)).to_vec();
        return all
            .into_iter()
            .map(|c| (c, Scope::Full, fields.clone()))
            .collect();
    };
    let Some(lines) = plan.get(family) else {
        return all.into_iter().map(|c| (c, Scope::Spot, vec![])).collect();
    };
    for (t, _, _) in lines {
        assert!(
            all.iter().any(|c| c.token == t),
            "RIBLT_PLAN: {family} has no construction {t}"
        );
    }
    all.into_iter()
        .flat_map(|c| {
            let scopes: Vec<(Scope, Option<Projections>)> = lines
                .iter()
                .filter(|(t, _, _)| *t == c.token)
                .map(|(_, s, f)| (*s, *f))
                .collect();
            // the construction is moved into its first planned line
            scopes
                .into_iter()
                .scan(Some(c), |c, (s, f)| {
                    c.take().map(|c| (c, s, f.into_iter().collect()))
                })
                .collect::<Vec<_>>()
        })
        .collect()
}

/// `G` hashing by one construction, straight to its addend: the adapter's
/// affine type is `G`'s addend, so `prepare` is a copy and the
/// construction's output is what the cells add.
#[derive(Clone, Copy)]
pub struct Hashed<'a, G: Accumulate> {
    pub group: G,
    pub h2c: &'a Construction<G>,
}

impl<G: Accumulate> Group for Hashed<'_, G> {
    type Affine = G::Addend;
    type Point = G::Point;
    fn identity(&self) -> G::Point {
        self.group.identity()
    }
    fn is_identity(&self, p: &G::Point) -> bool {
        self.group.is_identity(p)
    }
    fn to_affine(&self, p: &G::Point) -> G::Addend {
        self.group.prepare(&self.group.to_affine(p))
    }
    fn to_affine_batch(&self, ps: &[G::Point]) -> Vec<G::Addend> {
        self.group.prepare_batch(&self.group.to_affine_batch(ps))
    }
}

impl<G: Accumulate> Accumulate for Hashed<'_, G> {
    type Addend = G::Addend;
    fn prepare(&self, a: &G::Addend) -> G::Addend {
        *a
    }
    fn prepare_batch(&self, a: &[G::Addend]) -> Vec<G::Addend> {
        a.to_vec()
    }
    fn add(&self, p: &G::Point, a: &G::Addend) -> G::Point {
        self.group.add(p, a)
    }
    fn add_affine(&self, p: &G::Point, a: &G::Addend) -> G::Point {
        self.group.add(p, a)
    }
}

impl<G: Negate> Negate for Hashed<'_, G> {
    fn neg(&self, a: &G::Addend) -> G::Addend {
        self.group.neg_addend(a)
    }
    fn neg_addend(&self, a: &G::Addend) -> G::Addend {
        self.group.neg_addend(a)
    }
    fn equals_addend(&self, p: &G::Point, a: &G::Addend) -> bool {
        self.group.equals_addend(p, a)
    }
}

impl<G: Accumulate> HashToCurve for Hashed<'_, G> {
    fn hash(&self, h: &Salted, msg: &[u8]) -> G::Addend {
        (self.h2c.single)(&self.group, h, msg)
    }
    fn hash_batch(&self, h: &Salted, msgs: &[&[u8]]) -> Vec<G::Addend> {
        (self.h2c.batch)(&self.group, h, msgs)
    }
}
