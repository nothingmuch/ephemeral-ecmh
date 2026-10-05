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
    (r"^(negate|add|group\.\w+)$", "group ops"),
    (r"^(hash_to_curve|h2c|h2c\.id|h2c_parts|on_curve)$", "hash to curve"),
    (r"^riblt\.\w+$", "RIBLT workload"),
    (r"^digest$", "digest"),
    (r"^(curvegen/\w+|agm|zq|sieve/\w+)$", "curve generation"),
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
EXPERIMENTAL = "comparison maps"
MAPS = [
    # Pornin's map taken straight to a model's addend (benches/compare.rs,
    # compare109.rs): one map that also stands for the prepare
    (r"map x1 to \((x, λ|u, v)\)", "one map to addend"),
    (r"map x1|elligator2 x1|sswu x1", "one map"),
    (r"map x2", "two maps summed"),
]

LAYER_ORDER = [layer for _, layer in LAYERS] + [OTHER]
LAYER_ORDER.insert(LAYER_ORDER.index("hash to curve") + 1, EXPERIMENTAL)
# what each layer is, when its benches don't say
LAYER_NOTES = {
    "RIBLT workload": "Ratios against Go riblt-ecmh also include the mapping "
    "generator and the port to Rust, which are not attributed to this project "
    "(Workload, Comparison with the reference implementations).",
    EXPERIMENTAL: "Maps from digests to points other than try-and-increment, "
    "timed beside it: Pornin's map on binary curves, Elligator 2 on Edwards "
    "curves and simplified SWU on Weierstrass curves, once, and Pornin's map "
    "twice with the points summed; and Pornin's map taken straight to the "
    "unscaled or $\\lambda$-affine addend, where it stands for a hash and a prepare "
    "together, as try-and-increment to $(x, \\lambda)$ does. One map reaches at most "
    "$2^{m-1}$ points of $E[r]$, not uniformly. The checksum needs relations to be hard to "
    "find among hash outputs within the subset of the group that the hash "
    "reaches, not outputs uniform on the group (docs/problem.md, Adversary). "
    "Insertion estimates include each measured map with compatible preparation and "
    "addition; Pornin already returns an extended binary addend.",
}

# Every family by name: (base family, field width in bits, or None). A
# family other than its base is a variant: a smaller field (gf2_109,
# edwards107), gf2_127-lambda, the same curves under λ-projective
# accumulators, gf2_127-w, those with the w codec (its own hash and
# encoding), gf2_127-u, unscaled accumulators on the w codec's wire,
# gf2_122-gls, binary122's curves with their constant in GF(2^61) (and
# gf2_122-gls-lambda, both), or weier127-jacobian, the Weierstrass curves
# under Jacobian accumulators. xor is the SHA-256-only baseline.
BINARY = [(127, ""), (109, ""), (122, ""), (122, "-gls")]
ACCUMULATORS = ["", "-lambda", "-w", "-u"]
# Families of benchmarks that involve no curve or field: the RIBLT
# mapping's index generators (benches/riblt.rs, riblt.mapping).
UNGROUPED = {"mapping": ("mapping", None)}
FAMILIES = {
    **UNGROUPED,
    "xor": ("xor", None),
    "xor-siphash": ("xor", None),
    **{
        f"gf2_{bits}{gls}{acc}": ("gf2_127", bits)
        for bits, gls in BINARY
        for acc in ACCUMULATORS
    },
    "fp127": ("fp127", 127),
    "fp107": ("fp127", 107),
    "fp128": ("fp127", 128),
    "fp61x2": ("fp127", 122),
    "fp64x2": ("fp127", 128),
    "goldilocks2": ("fp127", 128),
    "edwards127": ("edwards127", 127),
    "edwards107": ("edwards127", 107),
    "twisted128": ("edwards127", 128),
    "edwards61x2": ("edwards127", 122),
    "twisted61x2": ("edwards127", 122),
    "twisted64x2": ("edwards127", 128),
    "twisted-goldilocks2": ("edwards127", 128),
    "weier127": ("weier127", 127),
    "weier127-jacobian": ("weier127", 127),
    "weier107": ("weier127", 107),
    "weier107-jacobian": ("weier127", 107),
    "weier61x2": ("weier127", 122),
    "weier61x2-jacobian": ("weier127", 122),
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
    **{
        f"binary{acc}.{bits}{gls}": f"gf2_{bits}{gls}{acc}"
        for bits, gls in BINARY
        for acc in ACCUMULATORS
    },
    **{f"edwards.{f}": f"edwards{f}" for f in ("107", "127", "61x2")},
    **{f"weier.{f}": f"weier{f}" for f in ("107", "127", "61x2")},
    **{f"weier-jacobian.{f}": f"weier{f}-jacobian" for f in ("107", "127", "61x2")},
    "twisted.128": "twisted128",
    # the a = -1 quotient before it was twisted128; the 2026-10-03 fixture
    # keeps edwards.128
    "edwards.128": "twisted128",
    "edwards128": "twisted128",
    "twisted.61x2": "twisted61x2",
    "twisted.64x2": "twisted64x2",
    "twisted.goldilocks2": "twisted-goldilocks2",
    "xor-sha256": "xor",
    "xor-sha256.64": "xor",
    "xor-siphash.64": "xor-siphash",
    "sha256": "xor",
    "gf2": "gf2_127",
    "fp": "fp127",
    "edwards": "edwards127",
    "weier": "weier127",
}
# (regex on the group, spelling, family): spellings some groups use for
# another family: Ristretto's input digest is a step of its hash.
CONTEXT_SPELLINGS = [
    (r"^h2c_parts$", "digest", "ristretto255"),
]
# groups whose functions name no family: agm and zq only count binary
# curves; riblt.mapping's functions name index generators
GROUP_FAMILIES = {"agm": "gf2_127", "zq": "gf2_127", "riblt.mapping": "mapping"}

# Available point-counting tools, not the provenance of a timing. The measured
# selection method comes from curvegen.csv; Rust also verifies every family's
# certificate. Sage scripts cover the families without a curvegen_times runner.
COUNTING_TOOLS = {
    "gf2_127": "Rust (AGM), PARI",
    "gf2_109": "Rust (AGM), PARI",
    "gf2_122": "Rust (AGM), PARI",
    "gf2_122-gls": "Rust (AGM), PARI",
    "edwards127": "PARI",
    "weier127": "PARI",
    "edwards107": "Sage",
    "weier107": "Sage",
    "twisted128": "Sage",
    "edwards61x2": "PARI",
    "weier61x2": "PARI",
    "twisted61x2": "PARI",
    "twisted64x2": "PARI",
    "twisted-goldilocks2": "PARI",
}


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
    "weier127": "#f28e2b",
    "edwards127": "#59a14f",
    "fp127": "#b6992d",
    "xor": "#79706e",
    "ristretto255": "#b07aa1",
    "secp256k1": "#9d7660",
    OTHER: "#79706e",
}
LIGHT = {
    "gf2_127": "#a0cbe8",
    "weier127": "#ffbe7d",
    "edwards127": "#8cd17d",
    "fp127": "#f1ce63",
    "xor": "#bab0ac",
    "ristretto255": "#d4a6c8",
    "secp256k1": "#d7b5a6",
    OTHER: "#bab0ac",
}
# legend and bar order: the baseline first
FAMILY_ORDER = [
    "xor",
    "gf2_127",
    "fp127",
    "edwards127",
    "weier127",
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
    (r"^group\.prepare$", r"", "prepare"),
    (r"^group\.add$", r"", "add"),
    (r"^group\.sub$", r"", "subtract"),
    (r"^group\.neg$", r"", "negate"),
    (r"^group\.is_identity$", r"", "is identity"),
    (r"^group\.equals$", r"", "equals addend"),
    (r"^group\.encode$", r"", "encode"),
    (r"^group\.decode$", r"", "decode"),
    (r"^h2c$", r"", "hash to curve"),
    # benches/group.rs's 32-byte IDs to the addend, under a construction
    # (h2c=) and a per-salt hash of the ID (proj=)
    (r"^h2c\.id$", r"", "ID to addend"),
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
    (r"^hash_to_curve", r"to \((x, λ|u, v)\)", "hash to addend"),
    (r"^hash_to_curve", r"", "hash to curve"),
    (r"^h2c_parts", r"^(?P<algo>[^/]+)/", "steps: {algo}"),
    (r"^on_curve", r"", "x on curve"),
    (r"^digest", r"", "digest"),
    (r"^riblt\.encode", r"", "encode (hash + cells)"),
    (r"^riblt\.cells", r"", "cell updates"),
    (r"^riblt\.peel", r"", "peel, per difference"),
    (r"^riblt\.stream", r"", "rateless encode and decode, per difference"),
    (r"^riblt\.mapping", r"/next\b", "mapping, per index"),
    (r"^riblt\.mapping", r"/item\b", "mapping, per item"),
    # the per-salt map digest of an ID (salted SHA-256 or a projection),
    # what each hash derives once per salt, and the ID they digest, which
    # outlives salts
    (r"^riblt\.mapping", r"/digest of id\b", "map digest of an ID"),
    (r"^riblt\.mapping", r"/keys\b", "per-salt keys"),
    (r"^riblt\.mapping", r"/id\b", "item ID"),
    (r"^riblt\.mapping", r"/digest\b", "map digest"),
    (r"^curvegen/verify_accept", r"", "accept certificate"),
    (r"^curvegen/verify_full", r"", "verify certificate"),
    (r"^curvegen/embedding", r"", "embedding-degree bound"),
    (r"^curvegen/find", r"agm", "find (Rust: AGM + sieve)"),
    (r"^curvegen/certify", r"agm", "certify (Rust: AGM + sieve + order witnesses)"),
    (r"^curvegen/count", r"", "point count (PARI)"),
    (r"^agm", r"order", "point count (Rust AGM)"),
    (r"^agm", r"", "AGM steps"),
    (r"^zq", r"", "Z_q ring op"),
    (r"^sieve/l", r"", "sieve, one l"),
    (r"^sieve/candidate", r"", "sieve, per candidate"),
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
    ("point prepare", "one at a time", "group ops", "prepare", ("per-element",)),
    ("point prepare", "batched", "group ops", "prepare", ("batch",)),
    ("is identity", None, "group ops", "is identity", None),
    ("equals addend", "match", "group ops", "equals addend", ("match",)),
    ("equals addend", "mismatch", "group ops", "equals addend", ("mismatch",)),
    # The unsplit per-element measurement contains mostly mismatches.
    ("equals addend", None, "group ops", "equals addend", ("per-element",)),
    ("encode", "one at a time", "group ops", "encode", ("per-element",)),
    ("encode", "batched", "group ops", "encode", ("batch",)),
    ("decode", "one at a time", "group ops", "decode", ("per-element",)),
    ("decode", "batched", "group ops", "decode", ("batch",)),
    (
        "hash to curve",
        "one at a time",
        "hash to curve",
        "hash to curve",
        ("per-element",),
    ),
    ("hash to curve", "batched", "hash to curve", "hash to curve", ("batch",)),
    (
        "hash to addend",
        "one at a time",
        "hash to curve",
        "hash to addend",
        ("per-element",),
    ),
    ("hash to addend", "batched", "hash to curve", "hash to addend", ("batch",)),
]
# A hash's cell is the fastest construction that yields its output, so the
# comparison maps' rows compete with try-and-increment's.
ALSO = {
    ("hash to curve", "hash to curve"): ("comparison maps", "one map"),
    ("hash to curve", "hash to addend"): ("comparison maps", "one map to addend"),
}
# The constructions benches/riblt.rs hashes by, as its h2c= parameter names
# them, and the rows that measure each: (layer, operation, id pattern).
H2C = {
    "ti": ("hash to curve", "hash to curve", r"^h2c/|try-and-increment(?! to)"),
    "ti-addend": ("hash to curve", "hash to addend", r"try-and-increment to"),
    "pornin": (EXPERIMENTAL, "one map", r"pornin"),
    "pornin-addend": (EXPERIMENTAL, "one map to addend", r"pornin"),
    "elligator2": (EXPERIMENTAL, "one map", r"elligator2"),
    "sswu": (EXPERIMENTAL, "one map", r"sswu"),
}
# How much of benches/riblt.rs a family gets, by its modeled insertion's ratio
# to the cheapest family's: every sweep for the leaders, a reduced sweep for
# the contenders behind them, and one point for the rest ("spot")
RIBLT_SCOPES = (("full", 1.15), ("buffer", 1.5))
# the section's columns: fields for the field rows, each coloured as its base
# field; curves for the rest, each -lambda beside the family it shares its
# curves with
FIELDS = {
    "gf2_127": "gf2_127",
    "gf2_109": "gf2_127",
    "gf2_122": "gf2_127",
    "fp127": "fp127",
    "fp107": "fp127",
    "fp128": "fp127",
    "fp61x2": "fp127",
    "fp64x2": "fp127",
    "goldilocks2": "fp127",
}
CURVES = [
    "xor",
    "xor-siphash",
    "gf2_127",
    "gf2_127-lambda",
    "gf2_127-w",
    "gf2_127-u",
    "gf2_109",
    "gf2_109-lambda",
    "gf2_109-w",
    "gf2_109-u",
    "gf2_122",
    "gf2_122-lambda",
    "gf2_122-w",
    "gf2_122-u",
    "gf2_122-gls",
    "gf2_122-gls-lambda",
    "gf2_122-gls-w",
    "gf2_122-gls-u",
    "edwards127",
    "edwards107",
    "twisted128",
    "weier127",
    "weier127-jacobian",
    "weier107",
    "weier107-jacobian",
    "edwards61x2",
    "weier61x2",
    "weier61x2-jacobian",
    "twisted61x2",
    "twisted64x2",
    "twisted-goldilocks2",
    "ristretto255",
    "secp256k1",
]

# The index generator of every RIBLT workload row (src/riblt.rs's default),
# whose walk the insertion estimate adds to the map digest.
MAPPING = "xoshiro256pp"

# The projection fields, as benches spell them in proj= and mapproj=
PROJECTIONS = {
    "fp130": "$\\mathbb{F}_{2^{130} - 5}$",
    "fp127": "$\\mathbb{F}_{2^{127} - 1}$",
    "gf2_127": "$\\mathrm{GF}(2^{127})$",
}

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
            ("xor-siphash", "SipHash-2-4 hashes XORed: no group, 8 bytes"),
        ],
    ),
    (
        "GF(2^127), 16 bytes",
        [
            ("gf2_127", "(X:S:Z:T) extended accumulators, Pornin's formulas"),
            ("gf2_127-lambda", "λ-projective (X:L:Z) accumulators"),
            ("gf2_127-w", "λ-projective, w codec: hashes and decodes to λ-affine"),
            ("gf2_127-u", "unscaled (X:S:Z), w codec, direct addend hashes"),
        ],
    ),
    (
        "GF(2^122) = GF(2^61)[u], dense constant, 16 bytes",
        [
            ("gf2_122", "(X:S:Z:T) extended"),
            ("gf2_122-lambda", "λ-projective"),
            ("gf2_122-w", "λ-projective, w codec"),
            ("gf2_122-u", "unscaled (X:S:Z), w codec"),
        ],
    ),
    (
        "GF(2^122), constant in GF(2^61) (GLS-shaped), 16 bytes",
        [
            ("gf2_122-gls", "(X:S:Z:T) extended"),
            ("gf2_122-gls-lambda", "λ-projective"),
            ("gf2_122-gls-w", "λ-projective, w codec"),
            ("gf2_122-gls-u", "unscaled (X:S:Z), w codec"),
        ],
    ),
    (
        "GF(2^109), 14 bytes",
        [
            ("gf2_109", "(X:S:Z:T) extended"),
            ("gf2_109-lambda", "λ-projective"),
            ("gf2_109-w", "λ-projective, w codec"),
            ("gf2_109-u", "unscaled (X:S:Z), w codec"),
        ],
    ),
    (
        "F_p, p = 2^127 − 1, 16 bytes",
        [
            ("edwards127", "a = 1 Edwards, extended += cached"),
            ("weier127", "short Weierstrass, projective += affine"),
            ("weier127-jacobian", "short Weierstrass, Jacobian += affine"),
        ],
    ),
    (
        "F_p, p = 2^107 − 1, 14 bytes",
        [
            ("edwards107", "a = 1 Edwards, extended += cached"),
            ("weier107", "short Weierstrass, projective += affine"),
            ("weier107-jacobian", "short Weierstrass, Jacobian += affine"),
        ],
    ),
    (
        "F_p, p = 2^128 − 275, 16 bytes",
        [
            (
                "twisted128",
                "a = −1 twisted Edwards modulo 2-torsion, extended += cached",
            ),
        ],
    ),
    (
        "F_{p²}, p = 2^61 − 1, 16 bytes",
        [
            ("edwards61x2", "a = 1 Edwards, signed-x codec, extended += cached"),
            ("weier61x2", "short Weierstrass, signed-x codec, projective += affine"),
            (
                "weier61x2-jacobian",
                "short Weierstrass, signed-x codec, Jacobian += affine",
            ),
            ("twisted61x2", "a = −1 Edwards modulo 2-torsion, extended += cached"),
        ],
    ),
    (
        "F_{p²}, p = 2^64 − 59, 16 bytes",
        [("twisted64x2", "a = −1 Edwards modulo 2-torsion, extended += cached")],
    ),
    (
        "F_{p²}, p = 2^64 − 2^32 + 1, 16 bytes",
        [
            (
                "twisted-goldilocks2",
                "a = −1 Edwards modulo 2-torsion, extended += cached",
            )
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
FIELD_GROUPS = [
    (
        "binary",
        [
            ("gf2_127", "F_2[z]/(z^127 + z^63 + 1)"),
            ("gf2_122", "GF(2^61)[u]/(u^2 + u + 1)"),
            ("gf2_109", "F_2[z]/(z^109 + z^5 + z^4 + z^2 + 1)"),
        ],
    ),
    (
        "prime",
        [
            ("fp127", "p = 2^127 − 1"),
            ("fp107", "p = 2^107 − 1"),
            ("fp128", "p = 2^128 − 275, for twisted128"),
        ],
    ),
    (
        "quadratic extensions of prime fields",
        [
            ("fp61x2", "GF(p^2) = F_p[i]/(i^2 + 1), p = 2^61 − 1"),
            ("fp64x2", "GF(p^2) = F_p[i]/(i^2 − 2), p = 2^64 − 59"),
            ("goldilocks2", "GF(p^2) = F_p[i]/(i^2 − 7), p = 2^64 − 2^32 + 1"),
        ],
    ),
]
# each curve's field, whose inversions its batches share
FIELD_OF = {
    "gf2_127": "gf2_127",
    "gf2_127-lambda": "gf2_127",
    "gf2_127-w": "gf2_127",
    "gf2_127-u": "gf2_127",
    "gf2_109": "gf2_109",
    "gf2_109-lambda": "gf2_109",
    "gf2_109-w": "gf2_109",
    "gf2_109-u": "gf2_109",
    "edwards127": "fp127",
    "weier127": "fp127",
    "edwards107": "fp107",
    "twisted128": "fp128",
    "weier107": "fp107",
    "weier127-jacobian": "fp127",
    "weier107-jacobian": "fp107",
    "edwards61x2": "fp61x2",
    "weier61x2": "fp61x2",
    "weier61x2-jacobian": "fp61x2",
    "twisted61x2": "fp61x2",
    "twisted64x2": "fp64x2",
    "twisted-goldilocks2": "goldilocks2",
} | {c: "gf2_122" for c in CURVES if c.startswith("gf2_122")}

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
OTHER_HARNESS = {
    "compare suite": "timed by benches/compare*.rs's hash_to_curve group, "
    "not benches/group.rs"
}


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
    if layer == "hash to curve":
        m, model = _first(MAPS, rest)
        if m:
            layer, op = EXPERIMENTAL, model or op
    _, mode = _first(MODES, rest)
    if mode is None:
        mode = "per-element" if elements else "total"
    return Facets(layer or OTHER, fam, base, bits, op, rank, mode)
