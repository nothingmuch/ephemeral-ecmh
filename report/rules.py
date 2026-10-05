"""Classification of benchmark identifiers into report axes.

Rules are evaluated in order, with the first match selecting the category.
An identifier without a layer, a family or an operation is unclassified,
and the report refuses to render until a rule names it.
"""

import re
from typing import NamedTuple

OTHER = "other"

# (regex on the group id, layer); one figure per layer, in this order.
LAYERS = [
    (r"^field$", "field"),
    (r"^(negate|add)$", "group ops"),
    (r"^(hash_to_curve|h2c|h2c_parts|on_curve)$", "hash to curve"),
    (r"^digest$", "digest"),
]

# (Regex on function/parameter, the operation the map is timed as.) The
# comparison maps are timed beside the groups' hash, try-and-increment, and
# compatible measured maps also supply insertion recipes. Primary rows include
# try-and-increment, Ristretto's hash_from_bytes and ElligatorSwift decoding
# on uniform inputs.
# A fixed-sign lift H = s H', with s in {+1, -1}, maps each signed relation
# among H outputs to a relation among H' outputs by multiplying coefficients
# by the corresponding signs. Fixed-sign x lifts are retained for the signed
# multiset-hash comparison under this relation.

LAYER_ORDER = [layer for _, layer in LAYERS] + [OTHER]
# what each layer is, when its benches don't say
LAYER_NOTES = {}

# Every family by name: (base family, field width in bits, or None). A
# family other than its base is a variant: a smaller field (gf2_109,
# edwards107), gf2_127-lambda, the same curves under λ-projective
# accumulators, gf2_127-w, those with the w codec (its own hash and
# encoding), gf2_127-u, unscaled accumulators on the w codec's wire,
# gf2_122-gls, binary122's curves with their constant in GF(2^61) (and
# gf2_122-gls-lambda, both), or weier127-jacobian, the Weierstrass curves
# under Jacobian accumulators. xor is the SHA-256-only baseline.
BINARY = [(127, "")]
ACCUMULATORS = ["", "-lambda", "-w", "-u"]
# Families of benchmarks that involve no curve or field: the RIBLT
# mapping's index generators (benches/riblt.rs, riblt.mapping).
UNGROUPED = {}
FAMILIES = {
    "xor": ("xor", None),
    "gf2_127": ("gf2_127", 127),
    "ristretto255": ("ristretto255", None),
    "secp256k1": ("secp256k1", None),
}
# How ids spell a family: the first word of the function id's first
# component, e.g. "binary-u.122-gls" in group.add/binary-u.122-gls/...,
# "gf2_127" in h2c_parts/gf2_127 t&i/... Each family's name spells it.
# benches/group.rs and benches/riblt.rs spell the binary families
# binary[-lambda|-w|-u].<bits>[-gls] and the rest <model>.<field>. Earlier
# runs, which the 2026-09-30 fixture keeps, spelt the 127-bit families
# without the width (gf2, fp, edwards, weier) and the XOR baseline sha256 or
# xor-sha256.
SPELLINGS = {
    **{name: name for name in FAMILIES},
    "xor-sha256": "xor",
    "sha256": "xor",
    "gf2": "gf2_127",
}
# (regex on the group, spelling, family): spellings some groups use for
# another family: Ristretto's input digest is a step of its hash.
CONTEXT_SPELLINGS = [
    (r"^h2c_parts$", "digest", "ristretto255"),
]
# groups whose functions name no family: agm and zq only count binary
# curves; riblt.mapping's functions name index generators
GROUP_FAMILIES = {}


def curve(family: str) -> str:
    """The curve under a family's accumulators and codec: gf2_127's four
    representations share one. This projection only joins shared selection
    data and presentation groups; rows keep their original representation IDs.
    GLS remains distinct because its curve selection is different. Compatibility
    spellings are handled by SPELLINGS rather than erased by this projection.
    """
    for s in ("-lambda", "-w", "-u", "-jacobian"):
        family = family.removesuffix(s)
    return family


# Tableau 20 (the 2016 palette): a hue per base family, so that a family and
# its variants read as one group in every figure, in two shades: dark for an
# item hashed or prepared one at a time, light for a batch. Blue, orange,
# green, yellow, purple and brown stay distinct under the common colour
# vision deficiencies; the baseline is grey.
COLORS = {
    "gf2_127": "#4e79a7",
    "xor": "#79706e",
    "ristretto255": "#b07aa1",
    "secp256k1": "#9d7660",
    OTHER: "#79706e",
}
LIGHT = {
    "gf2_127": "#a0cbe8",
    "xor": "#bab0ac",
    "ristretto255": "#d4a6c8",
    "secp256k1": "#d7b5a6",
    OTHER: "#bab0ac",
}
# legend and bar order: the baseline first
FAMILY_ORDER = [
    "xor",
    "gf2_127",
    "ristretto255",
    "secp256k1",
    OTHER,
]

# (regex on the group, regex on "function/parameter", operation). The
# operation may use the function regex's named groups, and {group}. Facets
# follow this order, then cost.
OPERATIONS = [
    (r"^field", r"batch invert", "batch invert"),
    # GF(2^122) uses qsolve for the quadratic equation; half-trace requires
    # an odd extension degree.
    (r"^field", r"\bqsolve\b", "halftrace"),
    (
        r"^field",
        (
            r"\b(?P<op>mul(_base|_u2)?|square|invert|sqrt(_ratio)?|halftrace|trace"
            r"|pow_p34|add|sub|neg|normalize|pack|unpack)\b"
        ),
        "{op}",
    ),
    # benches/group.rs: <layer>.<op>/<family>.<bits>/<parameters>
    (r"^h2c$", r"", "hash to curve"),
    (r"^negate", r"", "negate"),
    # combine_keys sums the whole slice in Jacobian coordinates and
    # normalizes once: a reduction, like the tree sums, not an add
    # extended accumulators to affine, sharing one inversion: not a sum
    (r"^add$", r"\bnormalize\b", "normalize"),
    (r"^add$", r"tree sum|batch|combine_keys", "batch sum"),
    (r"^add$", r"-=", "subtract"),
    (r"^add$", r"", "add"),
    # the λ families' hash straight to (x, λ) (benches/compare122.rs): its
    # output is an addend, so it stands for a hash and a prepare together
    (r"^hash_to_curve", r"", "hash to curve"),
    (r"^h2c_parts", r"^(?P<algo>[^/]+)/", "steps: {algo}"),
    (r"^on_curve", r"", "x on curve"),
    (r"^digest", r"", "digest"),
    (r"^agm", r"order", "point count (Rust AGM)"),
    (r"^agm", r"", "AGM steps"),
    (r"^zq", r"", "Z_q ring op"),
]

# (regex on "function/parameter", mode); the fallback is per-element when
# criterion has an element count, else total (one iteration).
MODES = [
    # group.equals: sums equal to the addend (pure cells) or not
    (r"\bmode=match\b", "match"),
    (r"\bmode=mismatch\b", "mismatch"),
    (r"latency|dependent|\b1 accumulator\b", "latency"),
    (r"throughput|\b\d+ chains\b|\b([2-9]|\d{2,}) accumulators\b", "throughput"),
    # riblt.peel spells batch=true|false; only true is batched
    (r"batch(?!=false)|tree sum|product tree|combine_keys", "batch"),
    (r"streaming", "streaming"),
]


# The report's first section: each elementary operation's cost on each curve.
# (operation, variant, layer, bench operation, modes or None for any); a cell
# is the cheapest per-element bench of that layer, operation and mode on that
# curve, so it picks the best representation each curve has (extended +=
# extended on gf2_127, += cached on Edwards). Field rows are per field: Edwards
# and Weierstrass share one. An operation's variants share a facet: the first
# is how it runs alone, the second how it runs in the pipeline (a dependent
# chain, or per item of a batch).
ELEMENTARY = [
    ("field mul", "throughput", "field", "mul", ("throughput",)),
    ("field mul", "latency", "field", "mul", ("latency",)),
    ("field square", "throughput", "field", "square", ("throughput",)),
    ("field square", "latency", "field", "square", ("latency",)),
    ("field add", None, "field", "add", None),
    ("field sqrt", None, "field", "sqrt", None),
    ("field sqrt ratio", None, "field", "sqrt_ratio", None),
    ("field halftrace", None, "field", "halftrace", None),
    ("field invert", "one at a time", "field", "invert", None),
    ("field invert", "batched", "field", "batch invert", None),
    ("field pack", None, "field", "pack", None),
    ("field unpack", None, "field", "unpack", None),
    ("point add", "throughput", "group ops", "add", ("throughput", "per-element")),
    ("point add", "latency", "group ops", "add", ("latency",)),
    ("point subtract", None, "group ops", "subtract", ("throughput", "per-element")),
    ("point negate", None, "group ops", "negate", None),
    # The unsplit per-element measurement contains mostly mismatches.
    (
        "hash to curve",
        "one at a time",
        "hash to curve",
        "hash to curve",
        ("per-element",),
    ),
    ("hash to curve", "batched", "hash to curve", "hash to curve", ("batch",)),
]
# A hash's cell is the fastest construction that yields its output, so the
# comparison maps' rows compete with try-and-increment's.
ALSO = {
    ("hash to curve", "hash to curve"): ("comparison maps", "one map"),
    ("hash to curve", "hash to addend"): ("comparison maps", "one map to addend"),
}
# the section's columns: fields for the field rows, each coloured as its base
# field; curves for the rest, each -lambda beside the family it shares its
# curves with
FIELDS = {
    "gf2_127": "gf2_127",
}
CURVES = [
    "xor",
    "ristretto255",
    "secp256k1",
]

# The index generator of every RIBLT workload row (src/riblt.rs's default),
# whose walk the insertion estimate adds to the map digest.
MAPPING = "xoshiro256pp"

# The references whose hash output is the operand of their addition, so an
# insertion is a hash and k additions with nothing to prepare. Neither hash
# has a batched form.
HASH_IS_ADDEND = ("xor", "ristretto255")

# Comparison rows grouped by field: (heading, [(family, note)]).
# Annotations distinguish models, parameters, coordinates and codecs.
CURVE_GROUPS = [
    (
        "baseline",
        [
            ("xor", "SHA-256 digests XORed: no group, 32 bytes"),
        ],
    ),
    (
        "references, 32 and 33 bytes",
        [
            ("ristretto255", "curve25519-dalek"),
            ("secp256k1", "libsecp256k1"),
        ],
    ),
]
FIELD_GROUPS = [("binary", [("gf2_127", "F_2[z]/(z^127 + z^63 + 1)")])]
# each curve's field, whose inversions its batches share
FIELD_OF = {} | {c: "gf2_122" for c in CURVES if c.startswith("gf2_122")}

# Marks on rows a table compares with the rest though they don't do the
# same work, {mark: why}. λ-projective (and so the w codec's) and Jacobian
# accumulators trade completeness for speed (src/curve/binary/lambda.rs,
# src/curve/weier/jacobian.rs).
INCOMPLETE = {
    "incomplete": "its additions are not complete: the identity and doubling "
    "take branches"
}
INCOMPLETE_FAMILIES = {
    f for f in FAMILIES if f.endswith(("-lambda", "-w", "-jacobian"))
}
# the λ families' hash straight to (x, λ) is timed by the compare suites,
# not benches/group.rs, which times every other term of an insertion


def marks(family: str) -> dict[str, str]:
    return dict(INCOMPLETE) if family in INCOMPLETE_FAMILIES else {}


class Facets(NamedTuple):
    layer: str
    family: str
    base: str
    bits: int | None
    operation: str
    op_rank: int
    mode: str


def _first(table, text):
    for pattern, value in table:
        m = re.search(pattern, text, re.IGNORECASE)
        if m:
            return m, value
    return None, None


def family(text: str, group: str = "") -> tuple[str, str, int | None]:
    """(family, base family, field bits or None) of a function id in a
    group, or OTHER's when SPELLINGS has no such family."""
    spelled = text.split("/", 1)[0].split(" ", 1)[0]
    context = (
        f for g, s, f in CONTEXT_SPELLINGS if s == spelled and re.search(g, group)
    )
    name = next(context, SPELLINGS.get(spelled) or GROUP_FAMILIES.get(group))
    if name is None:
        return OTHER, OTHER, None
    return name, *FAMILIES[name]


def classify(
    group: str, function: str | None, parameter: str | None, elements
) -> Facets:
    rest = "/".join(p for p in (function, parameter) if p)
    _, layer = _first(LAYERS, group)
    fam, base, bits = family(rest, group)
    op, rank = group, len(OPERATIONS)  # unclaimed: one facet per group
    for i, (gpat, fpat, template) in enumerate(OPERATIONS):
        m = re.search(fpat, rest, re.IGNORECASE)
        if m and re.search(gpat, group, re.IGNORECASE):
            op, rank = template.format(group=group, **m.groupdict()), i
            break
    _, mode = _first(MODES, rest)
    if mode is None:
        mode = "per-element" if elements else "total"
    return Facets(layer or OTHER, fam, base, bits, op, rank, mode)
