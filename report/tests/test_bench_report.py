import io
import json
import math
import re
import shutil
from pathlib import Path

import bench_report as br
import figures
import pandas as pd
import pytest
import rules
from criterion_fixture import read_rows, write_bench, write_tree

FIXTURE = Path(__file__).parent / "fixtures" / "2026-09-30.tsv"
ROWS = read_rows(FIXTURE)
LAYERS = [
    "field",
    "group ops",
    "hash to curve",
    "comparison maps",
    "digest",
    "curve generation",
]

# the workload, and the smaller-field variants the rules must name
EXTRA = [
    ("riblt.encode", "binary.109", "m=150,n=3500", 3500),
    ("add", "gf2_109/extended += affine, 8 accumulators", None, 1024),
    ("add", "edwards107/+= cached, 8 accumulators", None, 1024),
    ("hash_to_curve", "weier107/try-and-increment", None, 1024),
    ("hash_to_curve", "edwards107/montgomery try-and-increment", None, 1024),
    ("group.add", "binary-lambda.127/mode=throughput", None, 1024),
    ("group.add", "binary-w.127/mode=throughput", None, 1024),
    ("group.add", "binary.122/mode=throughput", None, 1024),
    ("group.add", "binary-lambda.122/mode=throughput", None, 1024),
    ("group.add", "binary.122-gls/mode=throughput", None, 1024),
    ("group.add", "binary-lambda.122-gls/mode=throughput", None, 1024),
    ("group.add", "binary-w.122/mode=throughput", None, 1024),
    ("group.add", "binary-w.122-gls/mode=throughput", None, 1024),
    # one family's pipeline, through the group traits
    ("h2c", "binary.127/mode=indep", None, 1024),
    ("h2c", "binary.127/mode=batch,n=1024", None, 1024),
    ("group.prepare", "binary.127/mode=indep", None, 1024),
    ("group.prepare", "binary.127/mode=batch,n=1024", None, 1024),
    ("group.add", "binary.127/mode=throughput", None, 1024),
    ("group.add", "binary.127/mode=latency", None, 1024),
    # a λ family's hash straight to its addend
    ("hash_to_curve", "gf2_122-lambda/try-and-increment to (x, λ)", None, 1024),
    (
        "hash_to_curve",
        "gf2_122-lambda/try-and-increment to (x, λ), batched",
        None,
        1024,
    ),
]


@pytest.fixture(scope="module")
def tree(tmp_path_factory):
    return write_tree(ROWS, tmp_path_factory.mktemp("criterion"))


@pytest.fixture(scope="module")
def table(tree):
    df, skipped = br.load(tree)
    assert skipped == []
    return br.tidy(df)


@pytest.fixture(scope="module")
def rendered(tmp_path_factory):
    root = write_tree(ROWS, tmp_path_factory.mktemp("criterion"))
    for g, f, p, n in EXTRA:
        write_bench(root, g, f, p, 90.0 * n, 100.0 * n, 110.0 * n, {"Elements": n})
    out = tmp_path_factory.mktemp("out")
    return out, br.report(root, out, formats=("svg",))


def test_per_element_values_match_the_summary(table):
    assert len(table) == len(ROWS)
    got = dict(zip(table.full_id, table.value_ns))
    for r in ROWS:
        full = "/".join(p for p in (r["group"], r["function"], r["parameter"]) if p)
        want = float(r["summary_ns"])
        # the summary is rounded to 2 decimals (ns) or 4 digits (ms)
        assert abs(got[full] - want) <= max(0.006, 5e-4 * want), full


def test_ids_survive_filename_mangling(table):
    # directories are criterion's filename-safe names; ids come from the json
    assert "on_curve/gf2/x: Tr(b/x) = 0 (1 I)" in set(table.full_id)
    row = table[table.full_id == "curvegen/verify_full/weier/256 rejections/3"].iloc[0]
    assert (row.function, row.parameter) == ("weier/256 rejections", "3")


def test_slope_else_mean_and_never_base(table):
    x = table.set_index("full_id")
    add = x.loc["add/gf2/batch affine (tree sum)"]
    assert add.stat == "slope"
    r = next(r for r in ROWS if r["function"] == "gf2/batch affine (tree sum)")
    assert add.estimate_ns == pytest.approx(float(r["point_ns"]))
    assert add.ci_lo_ns < add.estimate_ns < add.ci_hi_ns
    cert = x.loc["curvegen/verify_full/gf2/30 rejections/0"]
    assert (cert.stat, cert.unit, cert["mode"]) == ("mean", "ns/iter", "total")
    assert cert.value_ns == cert.estimate_ns


def test_stat_choice(tree):
    df, _ = br.load(tree, "median")
    assert set(df.stat) == {"median"}
    _, skipped = br.load(tree, "slope")  # flat-sampled benches have none
    assert len(skipped) == sum(not r["elements"] for r in ROWS)


@pytest.mark.parametrize(
    "t, want",
    [
        (None, (None, None)),
        ({"Elements": 64}, (64, None)),
        ({"ElementsAndBytes": {"elements": 4, "bytes": 144}}, (4, 144)),
        ({"Bytes": 36}, (None, 36)),
        ({"BytesDecimal": 36}, (None, 36)),
        ({"Bits": 256}, (None, 32)),
    ],
)
def test_throughput(t, want):
    assert br.throughput(t) == want


def test_current_suite_is_fully_classified(table):
    assert br.unclassified(table).empty
    assert set(table.layer) == set(LAYERS)
    assert set(table.family) == set(rules.FAMILY_ORDER) - {rules.OTHER}


@pytest.mark.parametrize(
    "full_id, want",
    [
        (
            "field/gf2/mul latency (dependent chain)",
            ("field", "gf2_127", "mul", "latency"),
        ),
        ("field/fp/mul throughput (8 chains)", ("field", "fp127", "mul", "throughput")),
        (
            "field/gf2/batch invert (product tree)",
            ("field", "gf2_127", "batch invert", "batch"),
        ),
        ("field/gf2/halftrace", ("field", "gf2_127", "halftrace", "per-element")),
        ("field/fp/add", ("field", "fp127", "add", "per-element")),
        (
            "add/gf2/extended += affine, 1 accumulator",
            ("group ops", "gf2_127", "add", "latency"),
        ),
        (
            "add/edwards/-= cached, 8 accumulators",
            ("group ops", "edwards127", "subtract", "throughput"),
        ),
        (
            "add/weier/batch affine (tree sum)",
            ("group ops", "weier127", "batch sum", "batch"),
        ),
        ("add/xor-sha256/xor 32B", ("group ops", "xor", "add", "per-element")),
        (
            "add/secp256k1/combine_keys (jacobian += affine)",
            ("group ops", "secp256k1", "batch sum", "batch"),
        ),
        (
            "negate/secp256k1/PublicKey::negate",
            ("group ops", "secp256k1", "negate", "per-element"),
        ),
        (
            "hash_to_curve/ristretto255/hash_from_bytes<Sha512>",
            ("hash to curve", "ristretto255", "hash to curve", "per-element"),
        ),
        # maps that are not uniform, or not proven so, stand apart
        (
            "hash_to_curve/gf2_127/pornin map x1, batched",
            ("comparison maps", "gf2_127", "one map", "batch"),
        ),
        (
            "hash_to_curve/gf2_127/pornin map x2",
            (
                "comparison maps",
                "gf2_127",
                "two maps summed",
                "per-element",
            ),
        ),
        (
            "hash_to_curve/gf2_109/pornin map x1",
            ("comparison maps", "gf2_109", "one map", "per-element"),
        ),
        # the map straight to an addend: a comparison map that is also a prepare
        (
            "hash_to_curve/gf2_127-u/pornin map x1 to (u, v), batched",
            ("comparison maps", "gf2_127-u", "one map to addend", "batch"),
        ),
        (
            "hash_to_curve/gf2_109-lambda/pornin map x1 to (x, λ)",
            ("comparison maps", "gf2_109-lambda", "one map to addend", "per-element"),
        ),
        (
            "hash_to_curve/gf2_122-gls/pornin map x1, batched",
            ("comparison maps", "gf2_122-gls", "one map", "batch"),
        ),
        (
            "hash_to_curve/edwards127/elligator2 x1",
            ("comparison maps", "edwards127", "one map", "per-element"),
        ),
        (
            "hash_to_curve/weier107/sswu x1",
            ("comparison maps", "weier107", "one map", "per-element"),
        ),
        (
            "hash_to_curve/weier127/sswu x1",
            ("comparison maps", "weier127", "one map", "per-element"),
        ),
        (
            "hash_to_curve/secp256k1/ellswift decode",
            ("hash to curve", "secp256k1", "hash to curve", "per-element"),
        ),
        (
            "hash_to_curve/sha256 (XOR baseline)",
            ("hash to curve", "xor", "hash to curve", "per-element"),
        ),
        (
            "h2c_parts/gf2 t&i/1. invert x, batched",
            ("hash to curve", "gf2_127", "steps: gf2 t&i", "batch"),
        ),
        (
            "on_curve/edwards/edwards u: Jacobi",
            ("hash to curve", "edwards127", "x on curve", "per-element"),
        ),
        (
            "digest/edwards/edwards-native streaming",
            ("digest", "edwards127", "digest", "streaming"),
        ),
        ("digest/gf2/batch", ("digest", "gf2_127", "digest", "batch")),
        # smaller fields: gf2_109, edwards107, weier107; 127 is the base
        (
            "add/gf2_109/extended += affine, 8 accumulators",
            ("group ops", "gf2_109", "add", "throughput"),
        ),
        ("digest/edwards107/batch", ("digest", "edwards107", "digest", "batch")),
        # GF(2^122): qsolve sits with the halftraces
        (
            "field/gf2_122/qsolve (z^2 + z = c, 2 base halftraces)",
            ("field", "gf2_122", "halftrace", "per-element"),
        ),
        (
            "field/gf2_122/mul_base (by a GF(2^61) constant) throughput (8 chains)",
            ("field", "gf2_122", "mul_base", "throughput"),
        ),
        (
            "hash_to_curve/weier107/try-and-increment",
            ("hash to curve", "weier107", "hash to curve", "per-element"),
        ),
        ("digest/weier127/batch", ("digest", "weier127", "digest", "batch")),
        (
            "field/gf2_122/mul_u2 (by a^2 = 1 + a) latency (dependent chain)",
            ("field", "gf2_122", "mul_u2", "latency"),
        ),
        (
            "field/gf2_109/normalize (to_u128)",
            ("field", "gf2_109", "normalize", "per-element"),
        ),
        # a batch to affine, not a sum
        (
            "add/gf2_122-gls/normalize, batched",
            ("group ops", "gf2_122-gls", "normalize", "batch"),
        ),
        (
            "agm/order",
            ("curve generation", "gf2_127", "point count (Rust AGM)", "per-element"),
        ),
        ("zq/sigma", ("curve generation", "gf2_127", "Z_q ring op", "per-element")),
        # benches/group.rs: layered ids, the family spelled <name>.<bits>
        (
            "group.add/binary.127/mode=throughput",
            ("group ops", "gf2_127", "add", "throughput"),
        ),
        (
            "group.sub/edwards.107/mode=latency",
            ("group ops", "edwards107", "subtract", "latency"),
        ),
        (
            "group.add/edwards.128/mode=throughput",
            ("group ops", "twisted128", "add", "throughput"),
        ),
        (
            "group.prepare/binary.109/mode=batch,n=1024",
            ("group ops", "gf2_109", "prepare", "batch"),
        ),
        ("group.neg/weier.107", ("group ops", "weier107", "negate", "per-element")),
        (
            "group.is_identity/weier.127",
            ("group ops", "weier127", "is identity", "per-element"),
        ),
        (
            "group.encode/edwards.127/mode=indep",
            ("group ops", "edwards127", "encode", "per-element"),
        ),
        (
            "group.encode/binary-w.122/mode=batch,n=1024",
            ("group ops", "gf2_122-w", "encode", "batch"),
        ),
        (
            "group.decode/binary.109/mode=indep",
            ("group ops", "gf2_109", "decode", "per-element"),
        ),
        (
            "h2c/binary.127/mode=batch,n=1024",
            ("hash to curve", "gf2_127", "hash to curve", "batch"),
        ),
        (
            "h2c/weier.107/mode=indep",
            ("hash to curve", "weier107", "hash to curve", "per-element"),
        ),
        # Identifiers with an explicit 127-bit field width.
        (
            "field/gf2_127/mul latency (dependent chain)",
            ("field", "gf2_127", "mul", "latency"),
        ),
        ("field/fp127/sqrt", ("field", "fp127", "sqrt", "per-element")),
        (
            "add/edwards127/+= cached, 8 accumulators",
            ("group ops", "edwards127", "add", "throughput"),
        ),
        ("digest/weier127/streaming", ("digest", "weier127", "digest", "streaming")),
        (
            "h2c_parts/gf2_127 t&i/1. invert x, batched",
            ("hash to curve", "gf2_127", "steps: gf2_127 t&i", "batch"),
        ),
        # binary127::lambda: gf2_127's curves, λ-projective accumulators
        (
            "group.add/binary-lambda.127/mode=throughput",
            ("group ops", "gf2_127-lambda", "add", "throughput"),
        ),
        (
            "h2c/binary-lambda.127/mode=batch,n=1024",
            ("hash to curve", "gf2_127-lambda", "hash to curve", "batch"),
        ),
        # riblt.peel's parameters say whether it peels in batches
        (
            "riblt.peel/binary.127/d=100,m=200,prefilter=true,batch=false",
            ("RIBLT workload", "gf2_127", "peel, per difference", "per-element"),
        ),
        (
            "riblt.peel/binary.127/d=100,m=200,prefilter=false,batch=true",
            ("RIBLT workload", "gf2_127", "peel, per difference", "batch"),
        ),
        # the two XOR baselines are told apart by their hash
        (
            "riblt.encode/xor-siphash.64/m=150,n=3500",
            ("RIBLT workload", "xor-siphash", "encode (hash + cells)", "per-element"),
        ),
        (
            "riblt.encode/edwards.128/m=150,n=3500",
            ("RIBLT workload", "twisted128", "encode (hash + cells)", "per-element"),
        ),
        (
            "riblt.encode/xor-sha256.64/m=150,n=3500",
            ("RIBLT workload", "xor", "encode (hash + cells)", "per-element"),
        ),
        # binary122: GLS constants, λ accumulators, and both
        (
            "group.add/binary.122/mode=latency",
            ("group ops", "gf2_122", "add", "latency"),
        ),
        (
            "group.add/binary.122-gls/mode=throughput",
            ("group ops", "gf2_122-gls", "add", "throughput"),
        ),
        (
            "group.add/binary-lambda.122/mode=throughput",
            ("group ops", "gf2_122-lambda", "add", "throughput"),
        ),
        (
            "h2c/binary-lambda.122-gls/mode=batch,n=1024",
            ("hash to curve", "gf2_122-gls-lambda", "hash to curve", "batch"),
        ),
        (
            "add/gf2_122-gls/extended += affine, 8 accumulators",
            ("group ops", "gf2_122-gls", "add", "throughput"),
        ),
        (
            "hash_to_curve/gf2_122-lambda/try-and-increment to (x, λ), batched",
            ("hash to curve", "gf2_122-lambda", "hash to addend", "batch"),
        ),
        # the F_{p^2} prototypes, and the codecs of the odd fields
        (
            "field/fp61x2/mul throughput (8 chains)",
            ("field", "fp61x2", "mul", "throughput"),
        ),
        ("field/fp64x2/sqrt", ("field", "fp64x2", "sqrt", "per-element")),
        (
            "field/fp64x2/batch invert (product tree)",
            ("field", "fp64x2", "batch invert", "batch"),
        ),
        ("field/fp61x2/pack", ("field", "fp61x2", "pack", "per-element")),
        (
            "field/fp64x2/square latency (dependent chain)",
            ("field", "fp64x2", "square", "latency"),
        ),
        # the Weierstrass curves' Jacobian families
        (
            "group.add/weier-jacobian.127/mode=throughput",
            ("group ops", "weier127-jacobian", "add", "throughput"),
        ),
        (
            "group.encode/weier-jacobian.107/mode=batch,n=64",
            ("group ops", "weier107-jacobian", "encode", "batch"),
        ),
        (
            "field/fp128/sqrt_ratio",
            ("field", "fp128", "sqrt_ratio", "per-element"),
        ),
        # Plonky3's fields
        (
            "field/goldilocks2/mul throughput (8 chains)",
            ("field", "goldilocks2", "mul", "throughput"),
        ),
        ("field/goldilocks2/unpack", ("field", "goldilocks2", "unpack", "per-element")),
        (
            "field/gf2_127/normalize (to_u128)",
            ("field", "gf2_127", "normalize", "per-element"),
        ),
        # binary109's λ and w families, and the λ one's hash to its addend
        (
            "group.add/binary-lambda.109/mode=throughput",
            ("group ops", "gf2_109-lambda", "add", "throughput"),
        ),
        (
            "group.decode/binary-w.109/mode=batch,n=64",
            ("group ops", "gf2_109-w", "decode", "batch"),
        ),
        (
            "hash_to_curve/gf2_109-lambda/try-and-increment to (x, λ)",
            ("hash to curve", "gf2_109-lambda", "hash to addend", "per-element"),
        ),
        (
            "h2c_parts/gf2_122-gls t&i/2. test Tr(b/x) (rejected)",
            ("hash to curve", "gf2_122-gls", "steps: gf2_122-gls t&i", "per-element"),
        ),
        # unclaimed: kept, under "other"
        ("mystery/thing", ("other", "other", "mystery", "per-element")),
    ],
)
def test_classify(full_id, want):
    group, function = full_id.split("/", 1)
    f = rules.classify(group, function, None, True)
    assert (f.layer, f.family, f.operation, f.mode) == want


@pytest.mark.parametrize(
    "old, new",
    [
        ("gf2/mul latency (dependent chain)", "gf2_127/mul latency (dependent chain)"),
        ("fp/sqrt (x^(2^125))", "fp127/sqrt (x^(2^125))"),
        ("edwards/+= cached, 1 accumulator", "edwards127/+= cached, 1 accumulator"),
        ("weier/x: Jacobi", "weier127/x: Jacobi"),
        ("gf2 pornin/2. invert m1 m2 m3", "gf2_127 pornin/2. invert m1 m2 m3"),
    ],
)
def test_both_spellings_name_one_family(old, new):
    assert rules.family(old) == rules.family(new)
    assert rules.family(new)[0] == rules.family(new)[1]  # the base, not a variant


@pytest.mark.parametrize(
    "group, function, parameter, operation",
    [
        ("curvegen/count", "edwards61x2", None, "point count (PARI)"),
        ("curvegen/verify_accept", "gf2_122-gls", "3", "accept certificate"),
        (
            "curvegen/verify_full",
            "twisted-goldilocks2/12 rejections",
            "0",
            "verify certificate",
        ),
        (
            "curvegen/find",
            "gf2_127 agm+sieve/index 30",
            "0",
            "find (Rust: AGM + sieve)",
        ),
        ("curvegen/embedding", "gf2_127", "0", "embedding-degree bound"),
    ],
)
def test_curvegen_names_each_family(group, function, parameter, operation):
    f = rules.classify(group, function, parameter, False)
    fam = function.split("/")[0].split(" ")[0]
    assert (f.layer, f.family, f.operation) == ("curve generation", fam, operation)


QUADRATIC_CURVES = [
    ("edwards.61x2", "edwards61x2", "edwards127", "fp61x2", 122),
    ("weier.61x2", "weier61x2", "weier127", "fp61x2", 122),
    ("weier-jacobian.61x2", "weier61x2-jacobian", "weier127", "fp61x2", 122),
    ("twisted.61x2", "twisted61x2", "edwards127", "fp61x2", 122),
    ("twisted.64x2", "twisted64x2", "edwards127", "fp64x2", 128),
    ("twisted.goldilocks2", "twisted-goldilocks2", "edwards127", "goldilocks2", 128),
]


@pytest.mark.parametrize("spelling, curve, base, field, bits", QUADRATIC_CURVES)
@pytest.mark.parametrize("group", ["group.add", "group.decode", "h2c"])
def test_quadratic_curve_identifiers(spelling, curve, base, field, bits, group):
    f = rules.classify(group, spelling + "/mode=indep", None, True)
    assert (f.family, f.base, f.bits) == (curve, base, bits)
    assert rules.family(curve) == (curve, base, bits)
    assert rules.FIELD_OF[curve] == field


@pytest.mark.parametrize(
    "spelling, curve, field, bits",
    [
        ("binary-u.127", "gf2_127-u", "gf2_127", 127),
        ("binary-u.109", "gf2_109-u", "gf2_109", 109),
        ("binary-u.122", "gf2_122-u", "gf2_122", 122),
        ("binary-u.122-gls", "gf2_122-gls-u", "gf2_122", 122),
    ],
)
@pytest.mark.parametrize("group", ["group.add", "group.decode", "h2c", "riblt.cells"])
def test_unscaled_binary_families_are_distinct(spelling, curve, field, bits, group):
    f = rules.classify(group, spelling + "/mode=throughput", None, True)
    assert (f.family, f.base, f.bits) == (curve, "gf2_127", bits)
    assert rules.family(curve) == (curve, "gf2_127", bits)
    assert rules.FIELD_OF[curve] == field
    assert curve in rules.CURVES
    assert any(name == curve for _, rows in rules.CURVE_GROUPS for name, _ in rows)


@pytest.mark.parametrize("acc", rules.ACCUMULATORS)
@pytest.mark.parametrize("bits, gls", rules.BINARY)
@pytest.mark.parametrize("group", ["group.add", "h2c", "riblt.peel"])
def test_binary_spellings_keep_their_width_and_representation(acc, bits, gls, group):
    f = rules.classify(group, f"binary{acc}.{bits}{gls}/mode=throughput", None, True)
    assert (f.family, f.base, f.bits) == (f"gf2_{bits}{gls}{acc}", "gf2_127", bits)


@pytest.mark.parametrize(
    "spelling", ["binary.131", "binary-v.127", "binary", "gf2_128", "edwards.61x3"]
)
def test_unlisted_spellings_name_no_family(spelling):
    assert rules.family(spelling + "/mode=indep") == (rules.OTHER, rules.OTHER, None)


def test_every_spelling_names_a_listed_family():
    for spelled, name in rules.SPELLINGS.items():
        assert rules.family(spelled + "/x") == (name, *rules.FAMILIES[name])
    assert set(rules.CURVES) <= set(rules.FAMILIES)
    assert set(rules.FAMILIES) <= set(rules.CURVES) | set(rules.FIELDS) | set(
        rules.UNGROUPED
    )


@pytest.mark.parametrize(
    "group, function, family, operation",
    [
        (
            "riblt.encode",
            "ristretto255/m=150,n=3500",
            "ristretto255",
            "encode (hash + cells)",
        ),
        (
            "riblt.peel",
            "ristretto255/d=4,m=10,prefilter=false,batch=false",
            "ristretto255",
            "peel, per difference",
        ),
        (
            "riblt.stream",
            "ristretto255/d=1000",
            "ristretto255",
            "rateless encode and decode, per difference",
        ),
        (
            "riblt.stream",
            "binary.127/d=10,map=mcg64",
            "gf2_127",
            "rateless encode and decode, per difference",
        ),
        (
            "riblt.encode",
            "xor-siphash.64/m=150,n=3500,map=mcg64",
            "xor-siphash",
            "encode (hash + cells)",
        ),
        ("riblt.mapping", "chacha8/next", "mapping", "mapping, per index"),
        ("riblt.mapping", "sha256-ctr/item,m=150", "mapping", "mapping, per item"),
    ],
)
def test_rateless_and_mapping_ids_are_classified(group, function, family, operation):
    f = rules.classify(group, function, None, True)
    assert (f.layer, f.family, f.operation) == ("RIBLT workload", family, operation)
    assert f.mode == "per-element"
    assert br.check_ids([f"{group}/{function}"]) == []


def test_explicit_square_throughput_is_present_in_elementary_table(tmp_path):
    root = tmp_path / "criterion"
    write_bench(
        root,
        "field",
        "gf2_127/square throughput (8 chains)",
        None,
        64.0,
        128.0,
        192.0,
        {"Elements": 64},
    )
    df, skipped = br.load(root)
    assert not skipped
    table = br.tidy(df)
    e = br.elementary(table)
    row = e[(e.facet == "field square") & (e.variant == "throughput")]
    assert len(row) == 1
    assert row.iloc[0].value_ns == pytest.approx(2.0)


def test_addend_equality_is_present_in_elementary_table(tmp_path):
    root = tmp_path / "criterion"
    write_bench(
        root,
        "group.equals",
        "binary.127",
        None,
        64.0,
        128.0,
        192.0,
        {"Elements": 64},
    )
    df, skipped = br.load(root)
    assert not skipped
    table = br.tidy(df)
    assert br.unclassified(table).empty
    e = br.elementary(table)
    row = e[e.facet == "equals addend"]
    assert len(row) == 1
    assert row.iloc[0].value_ns == pytest.approx(2.0)
    grid = br.elementary_grid(e, fields=False)
    key = ("equals addend", "")
    assert [column.key for column in grid.cols] == [key]
    assert grid.rows()[0].values[key] == pytest.approx(2.0)
    rendered = br.to_html([("grid", grid)])
    assert "equals addend" in rendered
    assert "group.equals/binary.127" in rendered


def test_quadratic_curve_headings_describe_structure_without_certification():
    for _, curve, _, _, _ in QUADRATIC_CURVES:
        assert curve in rules.CURVES
        rows = [
            (heading, note)
            for heading, rows in rules.CURVE_GROUPS
            for name, note in rows
            if name == curve
        ]
        assert len(rows) == 1
        assert "certified" not in " ".join(rows[0]).lower()
    headings = [
        heading
        for heading, rows in rules.FIELD_GROUPS
        if any(name == "fp61x2" for name, _ in rows)
    ]
    assert all("no curve yet" not in heading for heading in headings)


@pytest.mark.parametrize("name, bits", [("fp61x2", 122), ("fp64x2", 128)])
def test_quadratic_field_width_is_the_extension_width(name, bits):
    assert rules.family(name) == (name, "fp127", bits)


def test_variants_keep_their_family_and_base():
    lam = rules.classify("group.add", "binary-lambda.127/mode=latency", None, True)
    assert (lam.family, lam.base, lam.bits) == ("gf2_127-lambda", "gf2_127", 127)
    w = rules.classify("group.add", "binary-w.127/mode=latency", None, True)
    assert (w.family, w.base, w.bits) == ("gf2_127-w", "gf2_127", 127)
    gw = rules.classify("group.add", "binary-w.122-gls/mode=latency", None, True)
    assert (gw.family, gw.bits) == ("gf2_122-gls-w", 122)
    assert rules.family("gf2_122-gls-w/add")[0] == "gf2_122-gls-w"


def test_every_base_has_a_dark_and_a_light_shade():
    assert rules.COLORS.keys() == rules.LIGHT.keys()
    assert (
        len(set(rules.COLORS.values()) - {rules.COLORS[rules.OTHER]})
        == len(rules.COLORS) - 2
    )  # the baseline and OTHER share grey


def test_curve_columns_read_back_as_their_family():
    # the elementary section's column names read back as their family
    for curve in rules.CURVES:
        assert rules.family(curve)[0] == curve, curve


def test_classify_total_without_throughput():
    assert (
        rules.classify("curvegen/verify_full", "gf2/30 rejections", "0", False).mode
        == "total"
    )


def test_unclassified_ids_are_an_error(tmp_path):
    root = write_tree(ROWS, tmp_path / "criterion")
    odd = [
        ("mystery", "gf2/thing"),  # no layer
        ("field", "fq/mul"),  # no family
        ("riblt.encode", "binary.131/m=150,n=3500"),  # no such width
        ("field", "gf2_127/frobnicate"),  # no operation
    ]
    for g, f in odd:
        write_bench(root, g, f, None, 90.0, 100.0, 110.0, {"Elements": 1})
    with pytest.raises(SystemExit, match="4 benchmark ids match no rule") as err:
        br.report(root, tmp_path / "out", formats=("svg",))
    for g, f in odd:
        assert f"{g}/{f}" in str(err.value)
    assert not (tmp_path / "out").exists()


IDS = Path(__file__).parent / "fixtures" / "ids-2026-10-03.tsv"


def ids() -> list[list[str]]:
    lines = IDS.read_text().splitlines()
    return [line.split("\t") for line in lines if not line.startswith("#")]


def test_every_id_of_a_full_run_is_classified():
    rows = ids()
    assert len(rows) > 700
    for group, function in rows:
        f = rules.classify(group, function, None, True)
        assert rules.OTHER not in (f.layer, f.family), (group, function)
        assert f.op_rank < len(rules.OPERATIONS), (group, function)


def test_listed_ids_are_checked_before_a_run(monkeypatch):
    # criterion --list joins group and function with '/', as the groups do
    listed = ["/".join(row) for row in ids()]
    assert br.check_ids(listed) == []
    odd = ["mystery/gf2/thing", "field/fq/mul", "field/gf2_127/frobnicate"]
    assert br.check_ids(["", *listed[:3], *odd]) == odd
    monkeypatch.setattr("sys.stdin", io.StringIO("\n".join(listed[:3] + odd[:1])))
    with pytest.raises(SystemExit, match="1 benchmark ids match no rule") as err:
        br.main(["--check-ids"])
    assert "mystery/gf2/thing" in str(err.value)
    monkeypatch.setattr("sys.stdin", io.StringIO("\n".join(listed)))
    assert br.main(["--check-ids"]) is None
    monkeypatch.setattr("sys.stdin", io.StringIO("\n"))
    with pytest.raises(SystemExit, match="no benchmark ids"):
        br.main(["--check-ids"])


def test_coverage_shows_what_was_not_measured(rendered):
    out, t = rendered
    c = br.coverage(t)
    assert c.loc["gf2_127", "prepare"] == 2  # alone and batched
    assert c.loc["gf2_127", "equals addend"] == 0
    assert c.loc["gf2_109", "RIBLT, encode"] == 1
    assert (c.loc["twisted64x2"] == 0).all()
    rows = {r.label: r for r in br.coverage_grid(c).rows()}
    assert rows["gf2_127"].values["prepare"] == 2
    assert "equals addend" not in rows["gf2_127"].values
    assert rows["twisted64x2"].values == {}
    md = (out / "report.md").read_text()
    assert "## Coverage" in md and "| twisted64x2 (" in md
    assert (out / "coverage.csv").read_text().startswith("family,hash to curve,")


def test_report_files(rendered):
    out, t = rendered
    md, page = (out / "report.md").read_text(), (out / "report.html").read_text()
    # field and group operations are drawn once, in the elementary figures
    drawn = [layer for layer in LAYERS if layer not in ("field", "group ops")]
    for name in [figures.slug(layer) for layer in drawn + ["RIBLT workload"]] + [
        "elementary-field",
        "elementary-group",
    ]:
        assert (out / f"{name}.svg").exists(), name
        assert f"]({name}.svg)" in md and f'src="{name}.svg"' in page
    assert not (out / "field.svg").exists() and not (out / "group-ops.svg").exists()
    # every benchmark appears in the tables
    for label in t.label:
        assert br._md_cell(label) in md
    # a raw <Sha512> would render as an unknown HTML tag, i.e. nothing
    assert "hash_from_bytes&lt;Sha512&gt;" in md and "<Sha512>" not in md


def test_svg_is_deterministic(rendered, tmp_path):
    out, t = rendered
    br.plot_layer(t, "digest", tmp_path, ("svg",))
    assert (tmp_path / "digest.svg").read_bytes() == (out / "digest.svg").read_bytes()


def test_riblt_figure_with_a_family_missing_from_a_panel(rendered, tmp_path):
    # a spot-checked family lacks the sweeps the encode leaders have
    _, t = rendered
    w = t.layer == "RIBLT workload"
    enc = t[w & (t.operation == "encode (hash + cells)")]
    other = t[~t.family.isin(enc.family) & t.family.str.startswith("gf2_")].iloc[0]
    cells = enc.assign(operation="cell updates", family=other.family, base=other.base)
    br.plot_layer(pd.concat([t, cells]), "RIBLT workload", tmp_path, ("svg",))
    assert (tmp_path / f"{figures.slug('RIBLT workload')}.svg").exists()


def test_riblt_cost(rendered):
    out, t = rendered
    c = br.riblt_cost(t).set_index(["family", "hashing"])
    gf2 = c.loc[("gf2_127", "batched")]
    assert gf2["hash"] == "gf2/try-and-increment, batched"
    assert gf2["add"] == "gf2/extended += extended, 8 accumulators"
    k20 = br.mapping_degree(20)
    assert gf2["m=20"] == pytest.approx(gf2.hash_ns + k20 * gf2.add_ns)
    assert gf2["m=5"] == pytest.approx(96.94 + 2.9 * 14.36, abs=0.1)
    # The minimum is chosen across TAI and comparison maps.
    candidates = t[
        (t.family == "gf2_127")
        & t.layer.isin(["hash to curve", rules.EXPERIMENTAL])
        & (t["mode"] == "per-element")
    ]
    assert c.loc[("gf2_127", "one at a time"), "hash_ns"] == candidates.value_ns.min()
    assert (
        c.loc[("edwards127", "one at a time"), "add"]
        == "edwards/+= cached, 8 accumulators"
    )
    # combine_keys sums the whole batch: not an add, and secp256k1 has no other
    assert "secp256k1" not in {f for f, _ in c.index}
    assert ("edwards107", "one at a time") in c.index
    assert ("weier107", "one at a time") not in c.index  # no add bench
    assert {f for f, _ in c.index} >= {
        "xor",
        "gf2_127",
        "edwards127",
        "weier127",
        "ristretto255",
    }
    csv = (out / "riblt_cost_lower_bound.csv").read_text()
    assert csv.splitlines()[0].endswith("m=5,m=20,m=150,m=1350,m=12150")


def test_mapping_degree_is_twice_a_harmonic_tail():
    ks = [br.mapping_degree(m) for m in br.DEFAULT_MS]
    assert ks == pytest.approx([2.90, 5.29, 9.20, 13.57, 17.96], abs=0.005)
    assert [ld.key for ld in br.DEFAULT_LOADS] == [f"m={m}" for m in br.DEFAULT_MS]


def test_pipeline_cost(rendered):
    out, t = rendered
    c = br.pipeline_cost(t).set_index(["family", "hashing", "recipe"])
    # only families with group benches: the fixture's EXTRA binary.127, and
    # gf2_122-lambda, whose only recipe here is its hash to (x, λ); and
    # ristretto255, a reference. The fixture's XOR add has one accumulator,
    # no throughput, so XOR has no recipe.
    assert list(c.index) == [
        ("gf2_122-lambda", "one at a time", "hash to addend"),
        ("gf2_122-lambda", "batched", "hash to addend"),
        ("gf2_127", "one at a time", "hash + prepare"),
        ("gf2_127", "one at a time", "Pornin x1: hash to addend"),
        ("gf2_127", "one at a time", "Pornin x2: hash to addend"),
        ("gf2_127", "batched", "hash + prepare"),
        ("gf2_127", "batched", "Pornin x1: hash to addend"),
        ("ristretto255", "one at a time", "hash is the addend"),
    ]
    row = c.loc[("gf2_127", "batched", "hash + prepare")]
    # EXTRA benches take 100 ns per element
    assert (row.hash_ns, row.prepare_ns, row.add_ns) == (100, 100, 100)
    assert row["m=20"] == pytest.approx(100 + 100 + br.mapping_degree(20) * 100)
    fused = c.loc[("gf2_122-lambda", "batched", "hash to addend")]
    assert pd.isna(fused.prepare_ns)
    assert fused["m=20"] == pytest.approx(100 + br.mapping_degree(20) * 100)
    csv = (out / "riblt_cost.csv").read_text().splitlines()[0]
    assert csv.startswith("family,hashing,recipe,hash_ns,prepare_ns,add_ns,")
    assert "in the hash" in br.cost_table(c.reset_index()).prepare.tolist()
    md = (out / "report.md").read_text()
    assert "<details><summary>Lower bound, every family</summary>" in md


def test_pipeline_cost_reads_one_benchmark_per_term(tmp_path):
    root = tmp_path / "criterion"
    for g, f in [
        ("h2c", "binary.127/mode=indep"),
        ("group.prepare", "binary.127/mode=indep"),
        ("group.add", "binary.127/mode=throughput"),
        ("group.add", "binary.127/mode=throughput, again"),
    ]:
        write_bench(root, g, f, None, 90.0, 100.0, 110.0, {"Elements": 1})
    df, _ = br.load(root)
    with pytest.raises(ValueError, match="gf2_127 add: 2 benchmarks") as err:
        br.pipeline_cost(br.tidy(df))
    assert "group.add/binary.127/mode=throughput, again" in str(err.value)


def test_elementary(rendered):
    out, t = rendered
    e = br.elementary(t).set_index(["operation", "curve"])
    # the cheapest representation, in the named mode
    assert (
        e.loc[("point add, throughput", "gf2_127"), "bench"]
        == "add/gf2/extended += extended, 8 accumulators"
    )
    assert (
        e.loc[("point add, latency", "edwards127"), "bench"]
        == "add/edwards/+= cached, 1 accumulator"
    )
    assert (
        e.loc[("field mul, throughput", "fp127"), "bench"]
        == "field/fp/mul throughput (8 chains)"
    )
    assert (
        e.loc[("hash to curve, batched", "gf2_127"), "bench"]
        == "hash_to_curve/gf2/try-and-increment, batched"
    )
    assert e.loc[("hash to curve, one at a time", "xor"), "bench"].startswith(
        "hash_to_curve/sha256"
    )
    # a variant is its own row, and shares its operation's facet
    assert e.loc[("hash to curve, batched", "gf2_127"), "facet"] == "hash to curve"
    # a hash to (x, λ) is a hash and a prepare, not a hash to compare with
    # the others
    assert ("hash to curve, batched", "gf2_122-lambda") not in e.index
    assert (
        e.loc[("hash to addend, batched", "gf2_122-lambda"), "bench"]
        == "hash_to_curve/gf2_122-lambda/try-and-increment to (x, λ), batched"
    )
    # field rows are per field, not repeated per curve
    assert ("field mul, throughput", "edwards127") not in e.index
    assert ("point add, throughput", "fp127") not in e.index
    # λ: its own column, beside gf2_127 and in its colour
    lam = e.loc[("point add, throughput", "gf2_127-lambda")]
    assert lam.bench == "group.add/binary-lambda.127/mode=throughput"
    assert lam.base == "gf2_127"
    # rows grouped by field, each λ and w variant under its plain family
    g = br.elementary_grid(br.elementary(t), fields=False)
    curves = [r.label for r in g.rows()]
    assert curves.index("gf2_127-lambda") == curves.index("gf2_127") + 1
    assert curves.index("gf2_127-w") == curves.index("gf2_127") + 2
    assert curves[0] == "xor" and not g.rows()[0].compare
    fields = br.elementary_grid(br.elementary(t), fields=True)
    # the fixture has no gf2_109 or fp107 field benches
    assert [r.label for r in fields.rows()] == ["gf2_127", "fp127"]
    assert [h for h, _ in fields.groups] == ["binary", "prime"]
    # a binary field's add is an XOR, not timed; a prime field's is
    gf2, fp = fields.rows()
    assert gf2.values[br.FIELD_ADD] == "xor" and fp.values[br.FIELD_ADD] > 0
    assert not e.bench.str.contains("pornin|elligator2").any()
    md = (out / "report.md").read_text()
    order = ["## Insertion", "## Batching"]
    order += ["## Group operations", "## Field operations"]
    order += ["## All benchmarks", "### field", "### hash to curve"]
    order += ["### comparison maps", "### RIBLT workload", "### curve generation"]
    order += ["<summary>Lower bound, every family</summary>"]
    assert [md.index(h) for h in order] == sorted(md.index(h) for h in order)
    # the layers' figures are reference, folded away like their tables
    fig = re.search(
        r"<details><summary>Figure</summary>\n\n!\[([^]]+)\]\(digest.svg\)", md
    )
    assert fig and fig.group(1).startswith("digest: a bar per benchmark")
    for kind in ("field", "group"):
        assert (out / f"elementary-{kind}.svg").exists()
        assert f"](elementary-{kind}.svg)" in md
    # the long tables fold away
    assert "<details><summary>Which bench each cell is</summary>" in md
    assert (
        (out / "elementary.csv")
        .read_text()
        .startswith("operation,facet,variant,layer,")
    )


def test_unreadable_results_are_listed(tmp_path):
    root = tmp_path / "criterion"
    write_bench(root, "digest", "gf2/batch", None, 1.0, 2.0, 3.0, {"Elements": 2})
    d = write_bench(
        root, "digest", "gf2/streaming", None, 1.0, 2.0, 3.0, {"Elements": 2}
    )
    (d / "new" / "estimates.json").unlink()
    (root / "report").mkdir()  # criterion's own html report: no benchmark.json
    t = br.report(root, tmp_path / "out", loads=br.fixed_loads([2]), formats=("svg",))
    assert list(t.full_id) == ["digest/gf2/batch"]
    md = (tmp_path / "out" / "report.md").read_text()
    assert "## Skipped" in md and "gf2_streaming" in md
    assert str(tmp_path) not in md and f"benchmarks from run {tmp_path.name}." in md
    assert br.riblt_cost(t, (2,)).empty  # no hash-to-curve bench


@pytest.mark.parametrize(
    "damage",
    [
        lambda b: {k: v for k, v in b.items() if k != "group_id"},
        lambda b: [],
        lambda b: {**b, "throughput": {"Elements": "ten"}},
    ],
    ids=["missing field", "not an object", "bad throughput"],
)
def test_malformed_results_are_listed(tmp_path, damage):
    root = tmp_path / "criterion"
    write_bench(root, "digest", "gf2/batch", None, 1.0, 2.0, 3.0, {"Elements": 2})
    d = write_bench(root, "digest", "gf2/streaming", None, 1.0, 2.0, 3.0)
    bj = d / "new" / "benchmark.json"
    bj.write_text(json.dumps(damage(json.loads(bj.read_text()))))
    t = br.report(root, tmp_path / "out", loads=br.fixed_loads([2]), formats=("svg",))
    assert list(t.full_id) == ["digest/gf2/batch"]
    md = (tmp_path / "out" / "report.md").read_text()
    assert "## Skipped" in md and "gf2_streaming" in md
    assert str(tmp_path) not in md and f"benchmarks from run {tmp_path.name}." in md


def test_point_outside_its_ci(tmp_path):
    # criterion's slope is a least-squares fit, its CI bootstrap percentiles
    root = tmp_path / "criterion"
    write_bench(root, "digest", "gf2/x2", None, 1000.0, 950.0, 1100.0, {"Elements": 10})
    write_bench(
        root, "digest", "gf2/x3", None, 1000.0, 1200.0, 1100.0, {"Elements": 10}
    )
    t = br.report(root, tmp_path / "out", formats=("svg",))
    assert len(t) == 2 and (tmp_path / "out" / "digest.svg").exists()


def test_empty_tree_is_an_error(tmp_path):
    with pytest.raises(SystemExit):
        br.report(tmp_path, tmp_path / "out")


def test_cli(tmp_path):
    root = tmp_path / "criterion"
    write_bench(
        root,
        "add",
        "gf2/extended += affine, 8 accumulators",
        None,
        900.0,
        1000.0,
        1100.0,
        {"Elements": 10},
    )
    write_bench(
        root,
        "hash_to_curve",
        "gf2/try-and-increment",
        None,
        9e3,
        1e4,
        1.1e4,
        {"Elements": 10},
    )
    br.main([str(root), str(tmp_path / "out"), "--k", "3", "--formats", "png"])
    assert (tmp_path / "out" / "hash-to-curve.png").exists()
    assert not (tmp_path / "out" / "hash-to-curve.svg").exists()
    c = (tmp_path / "out" / "riblt_cost_lower_bound.csv").read_text().splitlines()
    assert c[0].endswith(",k=3") and c[1].endswith(",1300.0")


def test_no_machine_is_flagged(rendered):
    md = (rendered[0] / "report.md").read_text()
    assert "## Machine\n\nUnknown: no meta.json" in md


GLS122 = {"r": str(2**121 + 5), "cofactor": 2, "automorphisms": 4}


def test_runs_before_r_was_recorded_still_render(tmp_path):
    root = tmp_path / "criterion"
    write_bench(
        root, "group.add", "binary.127/mode=latency", None, 1, 2, 3, {"Elements": 1}
    )
    table = br.tidy(br.load(root)[0])
    meta = {"group_fixtures": {"binary.127": "certified"}}
    md = br.to_markdown(br.blocks(table, pd.DataFrame(), [], {}, [], "old", meta=meta))
    assert "## Machine" in md
    assert br.fixtures(table, meta) == {"binary.127": {"status": "certified"}}


@pytest.mark.parametrize("record", [GLS122, None])
def test_curve_parameters_are_not_benchmark_results(tmp_path, record):
    root = tmp_path / "criterion"
    write_bench(
        root, "group.add", "binary.122-gls/mode=latency", None, 1, 2, 3, {"Elements": 1}
    )
    df, skipped = br.load(root)
    assert not skipped
    table = br.tidy(df)
    meta = {"group_fixtures": {"binary.122-gls": record}} if record else None
    blocks = br.blocks(table, pd.DataFrame(), [], {}, [], "fixture", meta=meta)
    md, html = br.to_markdown(blocks), br.to_html(blocks)
    assert "Curve parameters" not in md and "Curve parameters" not in html
    assert "Certificate verified" not in md


def test_rho_is_sqrt_pi_r_over_2a():
    r = 2**126
    assert br.rho_bits(r, 2) == pytest.approx(math.log2(math.sqrt(math.pi * r / 4)))
    # a further sqrt 2 for an automorphism group of order 4
    assert br.rho_bits(r, 2) - br.rho_bits(r, 4) == pytest.approx(0.5)
    # docs/ecc_security.md: 59.8 for the GLS family
    assert f"{br.rho_bits(int(GLS122['r']), GLS122['automorphisms']):.1f}" == "59.8"


def test_fixtures_do_not_cover_unlisted_families(tmp_path):
    root = tmp_path / "criterion"
    for family in ["binary.127", "twisted.61x2"]:
        write_bench(root, "group.add", family, None, 1, 2, 3, {"Elements": 1})
    df, skipped = br.load(root)
    assert not skipped
    table = br.tidy(df)
    rec = {"r": str(2**126 + 1), "cofactor": 2, "automorphisms": 2}
    fixed = br.fixtures(table, {"group_fixtures": {"binary.127": rec}})
    assert fixed == {"binary.127": rec, "twisted.61x2": None}


def test_machine(tmp_path):
    run = tmp_path / "run"
    write_bench(run / "criterion", "digest", "gf2/batch", None, 1.0, 2.0, 3.0, None)
    meta = {
        "name": "box",
        "host": "box.local",
        "cpu": "Some CPU",
        "cores": 8,
        "memory": 16 * 2**30,
        "rustflags": "-C target-cpu=native",
        "target_features": ["avx", "pclmulqdq", "sse4.1"],
        "commit": "abc123",
        "rustc": "rustc 1.98.1 (48a229cea 2026-09-01)",
        "started": "2026-10-03T00:22:23Z",
        "profile": "quick",
    }
    (run / "meta.json").write_text(json.dumps(meta))
    br.report(run, tmp_path / "out", formats=("svg",))
    md = (tmp_path / "out" / "report.md").read_text()
    assert "## Machine" in md
    assert "| CPU | Some CPU, 8 cores, 16 GiB |" in md
    assert "| field backend features | pclmulqdq, sse4.1 |" in md
    assert "| commit | abc123 |" in md
    # published without the hostname, and named by the run
    published = json.loads((tmp_path / "out" / "meta.json").read_text())
    assert published == {k: v for k, v in meta.items() if k != "host"}
    assert "box.local" not in md and "benchmarks from run box." in md
    # every page and figure says which run, machine, build and sources
    line = "box · Some CPU · rustc 1.98.1 · RUSTFLAGS=-C target-cpu=native"
    line += " · abc123 · quick profile · 2026-10-03"
    assert md.endswith(f"---\n\n{line}\n") and "| profile | quick |" in md
    assert f"<footer>{line}</footer>" in (tmp_path / "out" / "report.html").read_text()
    assert line in (tmp_path / "out" / "digest.svg").read_text()
    dirty = run.parent / "dirty"
    write_bench(dirty / "criterion", "digest", "gf2/batch", None, 1.0, 2.0, 3.0, None)
    (dirty / "meta.json").write_text(json.dumps(meta | {"dirty": True}))
    br.report(dirty, tmp_path / "out3", formats=("svg",))
    md = (tmp_path / "out3" / "report.md").read_text()
    assert "| commit | abc123 (uncommitted edits) |" in md
    # so does the criterion directory, beside it; a bare one has none
    br.report(run / "criterion", tmp_path / "out1", formats=("svg",))
    assert "| CPU | Some CPU" in (tmp_path / "out1" / "report.md").read_text()
    bare = tmp_path / "bare"
    write_bench(bare, "digest", "gf2/batch", None, 1.0, 2.0, 3.0, None)
    br.report(bare, tmp_path / "out2", formats=("svg",))
    assert "Unknown: no meta.json" in (tmp_path / "out2" / "report.md").read_text()


def test_index_lists_the_runs_newest_first(tmp_path):
    runs = tmp_path / "runs"
    for name, started, cpu in [
        ("mac", "2026-10-02T09:00:00Z", "Apple M4"),
        ("x86", "2026-10-03T09:00:00Z", "AMD EPYC"),
    ]:
        d = runs / name / "report"
        d.mkdir(parents=True)
        (d / "report.html").write_text("<html></html>")
        meta = {
            "name": name,
            "started": started,
            "cpu": cpu,
            "commit": "abcdef0123456789",
        }
        (d / "meta.json").write_text(json.dumps(meta | {"dirty": name == "mac"}))
    br.main(["--index", str(runs)])
    page = (runs / "index.html").read_text()
    assert page.index('<a href="x86/report/report.html">x86</a>') < page.index(">mac<")
    assert "<td>abcdef012345 + edits</td>" in page and "<td>AMD EPYC</td>" in page
    md = (runs / "index.md").read_text()
    assert "| [x86](x86/report/report.md) | 2026-10-03T09:00:00Z | AMD EPYC |" in md
    with pytest.raises(SystemExit, match="no report.html"):
        br.index(tmp_path / "empty")


def test_labels_wrap_between_parameters_not_inside_them():
    s = "binary-lambda.122-gls/d=100,m=200,prefilter=true,batch=false"
    lines = figures.wrap_label(s, 30).split("\n")
    assert lines == [
        "binary-lambda.122-gls/d=100,",
        "m=200,prefilter=true,",
        "batch=false",
    ]
    assert (
        figures.wrap_label("gf2/batch affine (tree sum)")
        == "gf2/batch affine (tree sum)"
    )


def test_raw_tables_tell_operations_apart(tmp_path):
    # group.rs's ids are <op>/<family>/<parameters>: without the operation
    # a family's add and subtract rows read alike
    for group in ("group.add", "group.sub"):
        f = "binary.127/mode=throughput"
        write_bench(tmp_path, group, f, None, 1.0, 2.0, 3.0, {"Elements": 1})
    d = br.display_table(br.tidy(br.load(tmp_path)[0]))
    assert len(set(d.benchmark)) == 1
    assert list(d.operation) == ["add", "subtract"]


def test_fmt_time():
    assert [br.fmt_time(v) for v in (0.378, 23.11, 538.8, 4797.0, 1.2723e6)] == [
        "0.38 ns",
        "23.1 ns",
        "539 ns",
        "4.80 µs",
        "1.27 ms",
    ]
    assert [figures.fmt_tick(v) for v in (1, 100, 1e3, 1e7)] == [
        "1 ns",
        "100 ns",
        "1 µs",
        "10 ms",
    ]


def test_batch_change_says_which_way():
    assert figures.batch_change(300.0, 100.0) == "66.7% saved"
    assert figures.batch_change(100.0, 91.7) == "8.3% saved"
    assert figures.batch_change(100.0, 105.0) == "5.0% more"
    assert figures.batch_change(100.0, 100.0) == "0.0% saved"


def test_grid_marks_the_best_and_spares_the_baseline():
    from tables import FACTOR, Col, Grid, Row

    cols = [Col("t", "time"), Col("f", "speedup", kind=FACTOR)]
    rows = [
        Row("xor", values={"t": 1.0, "f": 99.0}, compare=False),
        Row("a", values={"t": 10.0, "f": 2.0}),
        Row("b", values={"t": 40.0, "f": 3.0}),
    ]
    g = Grid(cols, [("baseline", rows[:1]), ("curves", rows[1:])])
    assert g.best(cols[0]) == 10.0 and g.best(cols[1]) == 3.0
    md = tables_md(g)
    assert "| a | **10.0 ns** | 2.0× |" in md
    assert "| b | 40.0 ns | **3.0×** |" in md
    assert "| xor | 1.00 ns | 99.0× |" in md
    assert "| **curves** | | |" in md
    page = __import__("tables").to_html(g, br.fmt_time)
    assert (
        '<td class="time best" style="background:rgb(252,255,164);color:#000"' in page
    )
    assert "4.00× the best" in page
    assert '<td class="time baseline"' in page
    # what the marks mean, as text: a key under the grid, and in each cell
    assert page.endswith(f'<p class="key">{__import__("tables").KEY}</p>')
    assert '>10.0 ns<span class="sr"> (best)</span></td>' in page
    assert '>40.0 ns<span class="sr"> (4.00× the best)</span></td>' in page
    assert ">1.00 ns</td>" in page  # the baseline says nothing


def test_shading_keeps_the_text_legible_on_every_cell():
    tables = __import__("tables")

    def contrast(a, b):
        hi, lo = max(a, b), min(a, b)
        return (hi + 0.05) / (lo + 0.05)

    styles = [tables._shade(1.0 + i / 10) for i in range(200)]
    for style in styles:
        bg, text = re.fullmatch(
            r"background:rgb\((.*)\);color:(#000|#fff)", style
        ).groups()
        lum = tables._luminance([int(v) for v in bg.split(",")])
        assert contrast(lum, 0.0 if text == "#000" else 1.0) >= 4.5
    # saturated at SHADE_RATIO, and dark enough there to need white text
    assert styles[-1] == tables._shade(1e6) and styles[-1].endswith("#fff")


def test_grid_calls_overlapping_intervals_a_tie():
    from tables import Col, Grid, Row

    col = Col("t", "time")
    a = Row("a", values={"t": 10.0}, cis={"t": (9.0, 11.0)})
    b = Row("b", values={"t": 10.5}, cis={"t": (10.8, 12.0)})
    c = Row("c", values={"t": 12.0}, cis={"t": (11.5, 12.5)})
    g = Grid([col], [("", [a, b, c])])
    assert g.leaders(col) == {id(a), id(b)}
    md = tables_md(g)
    assert "| a | *10.0 ns* |" in md and "| b | *10.5 ns* |" in md
    assert "| c | 12.0 ns |" in md and "**" not in md
    page = __import__("tables").to_html(g, br.fmt_time)
    assert page.count('class="time tie"') == 2 and '"time best"' not in page
    assert '<span class="sr"> (tied for the best, within noise, 1.05×)</span>' in page
    # the interval behind the verdict, in the tooltip
    assert 'title="tied for the best, within noise; 9.00 ns to 11.0 ns"' in page
    assert 'title="1.20× the best; 11.5 ns to 12.5 ns"' in page
    # apart, the best wins alone
    b.cis["t"] = (11.2, 12.0)
    assert g.leaders(col) == {id(a)}
    assert "| a | **10.0 ns** |" in tables_md(g)


def test_sum_ci_adds_independent_half_widths_in_quadrature():
    def bench(v, lo, hi):
        return pd.Series({"value_ns": v, "value_lo_ns": lo, "value_hi_ns": hi})

    hsh, add = bench(100, 97, 104), bench(10, 9, 11)
    v, lo, hi = br.sum_ci([(hsh, 1), (None, 1), (add, 4)])
    assert v == 140
    assert (v - lo, hi - v) == pytest.approx((5.0, 5.657), abs=1e-3)
    assert all(map(math.isnan, br.sum_ci([(None, 1)])))


def tables_md(g):
    return __import__("tables").to_markdown(g, br.fmt_time, br._md_cell)


def test_batching_counts_the_inversions_shared():
    p = pd.DataFrame(
        [
            ("gf2_127", "one at a time", "hash + prepare", 900.0, 400.0, 14.0),
            ("gf2_127", "batched", "hash + prepare", 100.0, 4.0, 14.0),
        ],
        columns=["family", "hashing", "recipe", "hash_ns", "prepare_ns", "add_ns"],
    )
    e = pd.DataFrame(
        {
            "layer": ["field", "field"],
            "facet": ["field invert", "field invert"],
            "variant": ["one at a time", "batched"],
            "curve": ["gf2_127", "gf2_127"],
            "value_ns": [404.0, 4.0],
        }
    )
    (row,) = br.batching_grid(p, e).rows()
    assert row.values[("hash_ns", "inversions")] == 2.0
    assert row.values[("prepare_ns", "inversions")] == 1.0  # 396 / 400, to a tenth
    assert row.values[("hash_ns", "factor")] == 9.0
    # no codec benches, no codec columns
    assert not any(c.key[0] == "decode" for c in br.batching_grid(p, e).cols)
    codec = pd.DataFrame(
        {
            "layer": ["group ops"] * 2,
            "facet": ["decode"] * 2,
            "variant": ["one at a time", "batched"],
            "curve": ["gf2_127"] * 2,
            "value_ns": [450.0, 50.0],
        }
    )
    (row,) = br.batching_grid(p, pd.concat([e, codec])).rows()
    assert row.values[("decode", "inversions")] == 1.0
    assert row.values[("decode", "factor")] == 9.0


def test_insert_grid_sums_the_parts(rendered):
    out, t = rendered
    loads = br.fixed_loads([1, 8])
    g = br.insert_grid(br.pipeline_cost(t, loads), loads, "batched")
    rows = {(r.label, r.note): r for r in g.rows()}
    plain = next(r for (c, _), r in rows.items() if c == "gf2_127")
    assert plain.values["k=8"] == pytest.approx(100 + 100 + 8 * 100)
    # a λ family's hash to (x, λ) has no prepare
    fused = next(r for (c, _), r in rows.items() if c == "gf2_122-lambda")
    assert fused.values["prepare_ns"] == "in the hash"
    # unlike rows say how: λ's formulas, and the suite its fused hash is from
    assert set(fused.marks) == {"incomplete", "compare suite"}
    assert plain.marks == {}
    # each part names its benchmark and batch, the compare suites' too
    batch = "h2c/binary.127/mode=batch,n=1024, in batches of 1024"
    assert plain.tips["hash_ns"] == batch
    assert plain.tips["add_ns"] == "group.add/binary.127/mode=throughput"
    assert fused.tips["hash_ns"].endswith("to (x, λ), batched, in batches of 1024")
    assert "prepare_ns" not in fused.tips
    assert (out / "insert.svg").exists()
    md = (out / "report.md").read_text()
    assert "](insert.svg)" in md
    assert "| gf2_122-lambda [incomplete] [compare suite] (hashed straight" in md
    page = (out / "report.html").read_text()
    assert '<span class="mark" title="its additions are not complete' in page


def test_batch_sizes_sweep_but_the_largest_stands_for_batched(tmp_path):
    root = write_tree(ROWS, tmp_path / "criterion")
    sweep = [
        ("h2c", "binary.127/mode=indep", 1024, 800.0),
        ("h2c", "binary.127/mode=batch,n=8", 1024, 50.0),
        ("h2c", "binary.127/mode=batch,n=1024", 1024, 100.0),
        ("group.prepare", "binary.127/mode=batch,n=1024", 1024, 4.0),
        ("group.add", "binary.127/mode=throughput", 1024, 14.0),
    ]
    for g, f, n, per in sweep:
        write_bench(
            root, g, f, None, 0.9 * per * n, per * n, 1.1 * per * n, {"Elements": n}
        )
    df, _ = br.load(root)
    t = br.tidy(df)
    # the smaller batch is cheaper here, but "batched" means the largest
    (row,) = br.pipeline_cost(t)[
        lambda p: (
            (p.family == "gf2_127")
            & (p.hashing == "batched")
            & (p.recipe == "hash + prepare")
        )
    ].itertuples()
    assert row.hash_ns == pytest.approx(100.0)
    g = br.batch_size_grid(t)
    (r,) = [r for r in g.rows() if r.label == "gf2_127"]
    assert r.values[("h2c", 8)] == pytest.approx(50.0)
    assert r.values[("h2c", 1024)] == pytest.approx(100.0)
    assert r.values[("h2c", "alone")] == pytest.approx(800.0)


def test_references_insert_their_hash_output(tmp_path):
    root = tmp_path / "criterion"
    rows = [
        ("hash_to_curve", "ristretto255/hash_from_bytes<Sha512>", 5000.0),
        ("add", "ristretto255/+=, 1 accumulator", 70.0),
        ("add", "ristretto255/+=, 8 accumulators", 60.0),
        ("hash_to_curve", "sha256 (XOR baseline)", 20.0),
        # a single accumulator times no throughput: XOR stays out
        ("add", "xor-sha256/xor 32B", 0.5),
        ("h2c", "binary.127/mode=batch,n=1024", 100.0),
        ("group.prepare", "binary.127/mode=batch,n=1024", 4.0),
        ("group.add", "binary.127/mode=throughput", 14.0),
    ]
    for g, f, per in rows:
        n = 1024
        write_bench(
            root, g, f, None, 0.9 * per * n, per * n, 1.1 * per * n, {"Elements": n}
        )
    loads = [br.Load("k=10", "k = 10", 10)]
    p = br.pipeline_cost(br.tidy(br.load(root)[0]), loads)
    ref = p[p.recipe == "hash is the addend"]
    assert list(ref.family.astype(str)) == ["ristretto255"]
    r = ref.iloc[0]
    assert r.hashing == "one at a time" and pd.isna(r.prepare_ns)
    assert r["k=10"] == pytest.approx(5000.0 + 10 * 60.0)
    write_bench(
        root,
        "add",
        "xor-sha256/xor 32B, 8 accumulators",
        None,
        0.9 * 512,
        512,
        1.1 * 512,
        {"Elements": 1024},
    )
    p = br.pipeline_cost(br.tidy(br.load(root)[0]), loads)
    xor = p[p.family.astype(str) == "xor"].set_index("recipe")
    assert xor.loc["hash is the addend", "k=10"] == pytest.approx(20.0 + 10 * 0.5)
    fig = figures.plot_insert(p, loads, tmp_path, ["png"])
    assert "one item at a time" in fig.alt and "logarithmic" in fig.alt


def test_insertion_counts_the_mapping_only_where_every_load_has_it(tmp_path):
    root = tmp_path / "criterion"
    rows = [
        ("h2c", "binary.127/mode=batch,n=1024", 100.0),
        ("group.prepare", "binary.127/mode=batch,n=1024", 4.0),
        ("group.add", "binary.127/mode=throughput", 14.0),
        ("riblt.mapping", "salted-sha256/digest", 50.0),
        ("riblt.mapping", "xoshiro256pp/item,m=5", 7.0),
        # another generator's walk is not the workload's
        ("riblt.mapping", "chacha8/item,m=20", 1.0),
    ]
    for g, f, per in rows:
        n = 1024
        write_bench(
            root, g, f, None, 0.9 * per * n, per * n, 1.1 * per * n, {"Elements": n}
        )
    loads = br.riblt_loads([5, 20])
    p = br.pipeline_cost(br.tidy(br.load(root)[0]), loads)
    r = p.iloc[0]
    assert pd.isna(r["m=5 map_ns"]) and pd.isna(r["m=20 map_ns"])
    assert r["m=5"] == pytest.approx(104.0 + loads[0].k * 14.0)
    assert "were not timed" in br.insert_blocks(p, loads)[1][1]
    write_bench(
        root,
        "riblt.mapping",
        "xoshiro256pp/item,m=20",
        None,
        0.9 * 9 * 1024,
        9 * 1024,
        1.1 * 9 * 1024,
        {"Elements": 1024},
    )
    p = br.pipeline_cost(br.tidy(br.load(root)[0]), loads)
    r = p.iloc[0]
    assert r["m=5 map_ns"] == pytest.approx(57.0)
    assert r["m=20"] == pytest.approx(104.0 + loads[1].k * 14.0 + 59.0)
    assert "xoshiro256pp/item,m=20" in r["m=20 map bench"]
    assert "xoshiro256pp walk" in br.insert_blocks(p, loads)[1][1]
    fig = figures.plot_insert(p, loads, tmp_path, ["png"])
    assert "map digest" in fig.alt


def test_addend_equality_splits_matches_from_mismatches(tmp_path):
    root = tmp_path / "criterion"
    for mode, ns in (("match", 3.0), ("mismatch", 1.0)):
        write_bench(
            root,
            "group.equals",
            f"edwards.127/mode={mode}",
            None,
            64 * ns,
            128 * ns,
            192 * ns,
            {"Elements": 64},
        )
    df, skipped = br.load(root)
    assert not skipped
    table = br.tidy(df)
    assert br.unclassified(table).empty
    e = br.elementary(table)
    rows = e[e.facet == "equals addend"].set_index("variant").value_ns
    assert rows.to_dict() == pytest.approx({"match": 6.0, "mismatch": 2.0})
    grid = br.elementary_grid(e, fields=False)
    assert [column.key for column in grid.cols] == [
        ("equals addend", "match"),
        ("equals addend", "mismatch"),
    ]


# per-element ns: add, hash, prepare, encode, decode (None: not measured)
DECISION_COSTS = {
    "binary.127": (10, 100, 20, 30, 40),
    # slower adds, cheaper hashing: neither dominates the other
    "binary-w.127": (20, 90, 5, 30, 40),
    # slower everywhere, and lower rho
    "edwards.127": (30, 300, 30, 50, 60),
    # no encoding measured
    "weier.127": (10, 100, 20, None, 40),
}
DECISION_FIXTURES = {
    "binary.127": {"r": str(2**126 + 1), "cofactor": 2, "automorphisms": 2},
    "binary-w.127": {"r": str(2**126 + 1), "cofactor": 2, "automorphisms": 2},
    "edwards.127": {"r": str(2**125 + 1), "cofactor": 4, "automorphisms": 2},
    "weier.127": {"r": str(2**127 - 1), "cofactor": 1, "automorphisms": 2},
}


@pytest.fixture
def decided(tmp_path):
    return decide(tmp_path, DECISION_FIXTURES)


def test_a_run_without_r_omits_the_security_column(tmp_path):
    # runs before 2026-10-03 recorded only each fixture's status
    d = decide(tmp_path, {f: "certified" for f in DECISION_FIXTURES})
    assert d["rho"].isna().all() and d["dominated by"].isna().all()
    blocks = br.decision_blocks(d)
    md = br.to_markdown(blocks)
    assert "No fixture's $r$ was recorded, so no family is compared on security" in md
    assert "$\\log_2$ rho" not in md.split("\n| curve")[1]
    page = br._html_body(blocks)
    assert "security" not in page.split("<table")[1]
    # every row is "not compared", which stays, as it says something
    assert page.count(">not compared</td>") == len(d)


def decide(tmp_path, fixtures):
    run = tmp_path / "run"
    n = 1024
    batch = f"mode=batch,n={n}"
    for fam, costs in DECISION_COSTS.items():
        groups = ["group.add", "h2c", "group.prepare", "group.encode", "group.decode"]
        for g, v in zip(groups, costs):
            if v is not None:
                mode = "mode=throughput" if g == "group.add" else batch
                write_bench(
                    run / "criterion",
                    g,
                    f"{fam}/{mode}",
                    None,
                    0.9 * v * n,
                    v * n,
                    1.1 * v * n,
                    {"Elements": n},
                )
    (run / "curvegen.csv").write_text(
        "family,method,seed,find_s,verify_s\n"
        "gf2_127,agm+sieve,0,0.01,0.001\n"
        "gf2_127,agm+sieve,1,0.03,0.003\n"
        "gf2_127,pari,0,0.1,0.001\n"
        "gf2_127,pari,1,0.3,0.003\n"
    )
    t = br.tidy(br.load(run / "criterion")[0])
    meta = {"group_fixtures": fixtures}
    return br.decision(
        br.pipeline_cost(t),
        br.elementary(t),
        br.fixtures(t, meta),
        br.selection(run / "criterion"),
    )


def test_decision_dominance(decided):
    dom = decided["dominated by"]
    assert dom["gf2_127"] == [] and dom["gf2_127-w"] == []
    assert sorted(dom["edwards127"]) == ["gf2_127", "gf2_127-w"]
    # a missing measurement takes no part, and leaves its sums unknown
    assert dom["weier127"] is None
    assert math.isnan(decided.loc["weier127", "round"])


def test_decision_regimes(decided):
    r = decided.loc["gf2_127"]
    k = br.mapping_degree(br.ROUND_M)
    assert br.ROUND_M == 1350
    assert r["hash + prepare"] == pytest.approx(120)
    assert r.retained == pytest.approx(120 + k * 10)
    assert r["round"] == pytest.approx(
        br.ROUND_N * k * 10 + br.ROUND_D * 120 + br.ROUND_M * (30 + 40)
    )
    assert r["round lo"] < r["round"] < r["round hi"]
    assert r.rho == pytest.approx(62.8, abs=0.05)
    assert decided.loc["edwards127", "rho"] == pytest.approx(62.3, abs=0.05)


def test_decision_selection_is_the_faster_method(decided):
    r = decided.loc["gf2_127"]
    assert r.find == pytest.approx(0.02e9) and r.verify == pytest.approx(0.002e9)
    assert r.seeds == 2 and r.method == "agm+sieve"
    assert r["counting tools"] == "Rust (AGM), PARI"
    # shared by the curve's representations; none timed for edwards127
    assert decided.loc["gf2_127-w", "find"] == pytest.approx(0.02e9)
    e = decided.loc["edwards127"]
    assert math.isnan(e.find) and e["counting tools"] == "PARI"
    assert pd.isna(e.method)


def test_decision_leads_the_report(decided):
    blocks = br.decision_blocks(decided)
    md = br.to_markdown(blocks)
    assert md.startswith("## Families compared")
    assert "measured method" in md and "available tools" in md
    assert "agm+sieve" in md and "Rust (AGM), PARI" in md
    assert "gf2_127, gf2_127-w" in md and "not compared" in md and "none" in md
    # selection means are shown, never the best: they carry no interval
    assert "20.0 ms" in md and "**20.0 ms**" not in md
    page = br._html_body(blocks)
    assert (
        '<td class="plain-time" title="agm+sieve, mean of 2 seeds">20.0 ms</td>' in page
    )
    # the summed estimates name the benchmarks of their terms
    assert 'title="best; ' in page
    assert "k × (group.add/binary.127/mode=throughput)" in page
    assert "m × (group.encode/binary.127/mode=batch,n=1024; group.decode/" in page
    assert "d × (not measured)" not in page


def test_decision_states_the_share_of_addition(decided):
    md = br.to_markdown(br.decision_blocks(decided))
    share = (
        br.ROUND_N * br.mapping_degree(br.ROUND_M) * decided["add"] / decided["round"]
    )
    assert f"{100 * share.min():.0f} to {100 * share.max():.0f}%" in md


def test_decision_without_batched_rows_is_empty():
    p = pd.DataFrame({"hashing": ["one at a time"], "family": ["gf2_127"]})
    assert br.decision(p, pd.DataFrame(), {}, br.selection(Path("/nonexistent"))).empty


def test_dominance_needs_evidence_on_every_axis(decided):
    d = decided.copy()
    # Overlapping decode intervals prevent a dominance claim.
    d.loc["edwards127", "decode lo"] = 1.0
    assert br.dominators(d)["edwards127"] == []
    # unless both read one measurement
    d.loc["edwards127", "decode bench"] = d.loc["gf2_127", "decode bench"]
    assert br.dominators(d)["edwards127"] == ["gf2_127"]


@pytest.mark.parametrize(
    "spelling", ["edwards128", "edwards.128", "twisted128", "twisted.128"]
)
def test_twisted128_historical_names(spelling):
    f = rules.classify("group.add", f"{spelling}/mode=throughput", None, True)
    assert f.family == "twisted128"
    assert f.base == "edwards127"
    assert f.bits == 128


def test_historical_selection_name_is_canonical(tmp_path):
    (tmp_path / "curvegen.csv").write_text(
        "family,method,seed,find_s,verify_s\nedwards128,pari,0,0.1,0.001\n"
    )
    row = br.selection(tmp_path).iloc[0]
    assert row.family == "twisted128"
    assert row.find_s == 0.1


def test_curve_selection_shows_every_method_and_the_fastest_drawn(tmp_path):
    phases = "candidate,quick,sieve,count,factor,witness,accept"
    cols = ",".join(f"prove_{p}_s" for p in phases.split(","))
    (tmp_path / "curvegen.csv").write_text(
        "family,method,seed,candidates,prove_counts,find_s,prove_s,verify_s,"
        + cols
        + "\n"
        "gf2_127,agm+sieve,0,30,8,0.01,0.02,0.0001,0,0,0.005,0.012,0,0,0.001\n"
        "gf2_127,agm+sieve,1,50,10,0.03,0.04,0.0003,0,0,0.01,0.025,0,0,0.001\n"
        "gf2_127,pari,0,30,30,0.2,0.3,0.0001,0,0,0,0.29,0,0,0.001\n"
        "weier127,pari+sieve,0,200,20,1.0,1.3,0.0002,0.01,0,0.1,1.1,0.05,0,0.01\n"
    )
    runs = br.selection_runs(tmp_path)
    assert set(br.fastest_method(runs).method) == {"agm+sieve", "pari+sieve"}
    _, _, table = br.selection_blocks(runs)
    t = table[1].set_index(["family", "method"])
    assert t.loc[("gf2_127", "agm+sieve"), "shown"] == "yes"
    assert t.loc[("gf2_127", "pari"), "shown"] == ""
    assert t.loc[("gf2_127", "agm+sieve"), "prove"] == "30.0 ms"
    # 18.5 ms of 30 ms counting
    assert t.loc[("gf2_127", "agm+sieve"), "counting, of prove"] == "62%"
    fig = figures.plot_selection(br.fastest_method(runs), tmp_path, ["png"])
    assert "certify" in fig.alt and (tmp_path / "selection.png").is_file()


def test_embedding_share_relates_the_bound_to_acceptance(tmp_path):
    root = tmp_path / "criterion"
    for fam, embedding, accept in [("gf2_127", 30e3, 60e3), ("edwards127", 20e3, 80e3)]:
        write_bench(
            root, "curvegen/embedding", fam, "0", embedding, embedding, embedding
        )
        write_bench(root, "curvegen/verify_accept", fam, "0", accept, accept, accept)
    # timed only for acceptance: not part of the comparison
    write_bench(root, "curvegen/verify_accept", "weier127", "0", 1e6, 1e6, 1e6)
    t = br.tidy(br.load(root)[0])
    text = br.embedding_share(t)
    assert "25% to 50%" in text
    assert "1 ms" not in text and "1.00 ms" not in text
    assert br.embedding_share(t[t.group != "curvegen/embedding"]) == ""


RUN_META = {
    "name": "fixture-run",
    "host": "box.local",
    "started": "2026-10-03T11:55:20Z",
    "cpu": "Apple M4",
    "os": "macOS 26.6.2",
    "arch": "arm64",
    "rustc": "rustc 1.98.1 (48a229cea 2026-09-01)",
    "rustflags": "-C target-cpu=apple-m4",
    "commit": "ebe05c50203bacba65853a298692e4e5a1c4537a",
    "dirty": False,
    "profile": "full",
    "group_fixtures": DECISION_FIXTURES,
}
CURVEGEN = (
    "family,method,seed,find_s,verify_s\n"
    "gf2_127,agm+sieve,0,0.01,0.001\n"
    "gf2_127,pari,0,0.1,0.001\n"
)


def write_run(run: Path) -> Path:
    """A bench-run directory: the fixture rows, EXTRA, the decision
    table's families, meta.json and curvegen.csv."""
    root = write_tree(ROWS, run / "criterion")
    for g, f, p, n in EXTRA:
        write_bench(root, g, f, p, 90.0 * n, 100.0 * n, 110.0 * n, {"Elements": n})
    n = 1024
    groups = ["group.add", "h2c", "group.prepare", "group.encode", "group.decode"]
    for fam, costs in DECISION_COSTS.items():
        for g, v in zip(groups, costs):
            if v is not None:
                mode = "mode=throughput" if g == "group.add" else f"mode=batch,n={n}"
                b = (0.9 * v * n, v * n, 1.1 * v * n)
                write_bench(root, g, f"{fam}/{mode}", None, *b, {"Elements": n})
    (run / "meta.json").write_text(json.dumps(RUN_META))
    (run / "curvegen.csv").write_text(CURVEGEN)
    return run


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    return write_run(tmp_path_factory.mktemp("bench-runs") / "fixture-run")


def test_export_keeps_the_raw_estimates(run, tmp_path):
    dest = tmp_path / "results" / "fixture-run"
    br.main(["--export", str(run), str(dest)])
    assert sorted(p.name for p in dest.iterdir()) == [
        "benchmarks.csv",
        "curvegen.csv",
        "meta.json",
    ]
    df, _ = br.load(run / "criterion")
    raw = pd.read_csv(dest / "benchmarks.csv", keep_default_na=False)
    assert list(raw.columns) == br.RAW
    assert list(raw.group) == list(df.group)
    assert list(raw.estimate_ns) == list(df.estimate_ns)
    # element counts stay integers, and a benchmark without one has none
    lines = (dest / "benchmarks.csv").read_text().splitlines()
    assert ",1024," in lines[1] or ",1024," in lines[-1]
    assert ".0," not in "".join(ln.split(",")[-2] for ln in lines[1:])
    assert json.loads((dest / "meta.json").read_text()) == {
        k: v for k, v in RUN_META.items() if k != "host"
    }
    assert (dest / "curvegen.csv").read_text() == CURVEGEN
    with pytest.raises(SystemExit, match="exists"):
        br.export(run, dest)


LOAD_TSV = (
    "time\tphase\tload1\tload5\tload15\trunnable\n"
    "100\tidle before\t0.50\t0.5\t0.5\t\n"
    "105\tidle before\t0.50\t0.5\t0.5\t\n"
    "110\tgroup\t0.58\t0.5\t0.5\t\n"
    "115\tgroup\t0.65\t0.5\t0.5\t\n"
)


def test_machine_load_recovers_the_runnable_threads(tmp_path):
    (tmp_path / br.LOAD).write_text(LOAD_TSV)
    load = br.machine_load(tmp_path)
    assert list(load.index) == ["idle before", "group"]
    assert load.seconds.tolist() == [5, 10]
    d = math.exp(-5 / 60)
    assert load.runnable["idle before"] == pytest.approx(0.5)
    first, second = ((b - a * d) / (1 - d) for a, b in [(0.5, 0.58), (0.58, 0.65)])
    assert load.runnable["group"] == pytest.approx((first + second) / 2)
    assert load.runnable_max["group"] == pytest.approx(max(first, second))
    (tmp_path / br.LOAD).write_text(LOAD_TSV.replace("\t0.5\t\n", "\t0.5\t3\n"))
    assert br.machine_load(tmp_path).runnable.tolist() == [3, 3]


def test_export_keeps_the_load_and_the_report_shows_it(run, tmp_path):
    (run / br.LOAD).write_text(LOAD_TSV)
    dest = tmp_path / "results" / "fixture-run"
    br.export(run, dest)
    assert (dest / br.LOAD).read_text() == LOAD_TSV
    br.report(dest, tmp_path / "out", formats=("svg",))
    page = (tmp_path / "out" / "report.md").read_text()
    assert "runnable, mean" in page and "idle before" in page


def test_export_refuses_what_its_report_would_lack(tmp_path):
    run = tmp_path / "run"
    write_bench(run / "criterion", "digest", "gf2/batch", None, 1.0, 2.0, 3.0, None)
    with pytest.raises(SystemExit, match="no meta.json"):
        br.export(run, tmp_path / "a")
    (run / "meta.json").write_text(json.dumps(RUN_META))
    d = write_bench(run / "criterion", "digest", "gf2/streaming", None, 1, 2, 3)
    (d / "new" / "estimates.json").unlink()
    with pytest.raises(SystemExit, match="gf2_streaming"):
        br.export(run, tmp_path / "b")
    with pytest.raises(SystemExit, match="no criterion results"):
        br.export(tmp_path / "empty", tmp_path / "c")
    shutil.rmtree(d)
    (run / "meta.json").write_text(json.dumps({**RUN_META, "dirty": True}))
    with pytest.raises(SystemExit, match="committed source"):
        br.export(run, tmp_path / "d")
    assert not any((tmp_path / n).exists() for n in "abcd")


def test_compact_form_renders_the_same_report(run, tmp_path):
    dest = tmp_path / "results" / "fixture-run"
    br.export(run, dest)
    a, b = tmp_path / "from-criterion", tmp_path / "from-compact"
    ta = br.report(run, a, formats=("svg",))
    tb = br.report(dest, b, formats=("svg",))
    assert ta.elements.isna().any() and len(ta) == len(tb)
    names = sorted(p.name for p in a.iterdir())
    assert names == sorted(p.name for p in b.iterdir())
    assert {"decision.csv", "insert.svg", "meta.json"} <= set(names)
    for name in names:
        assert (a / name).read_bytes() == (b / name).read_bytes(), name
    md = (b / "report.md").read_text()
    assert "## Families compared" in md and "fixture-run · Apple M4" in md


def test_compact_form_holds_one_statistic(run, tmp_path):
    dest = tmp_path / "fixture-run"
    br.export(run, dest)
    with pytest.raises(SystemExit, match="typical estimate only"):
        br.report(dest, tmp_path / "out", stat="median")


@pytest.mark.parametrize(
    "family,map_name,group_name,recipe,needs_prepare",
    [
        ("gf2_127", "pornin map x1", "binary.127", "Pornin x1: hash to addend", False),
        ("gf2_127", "pornin map x2", "binary.127", "Pornin x2: hash to addend", False),
        (
            "edwards127",
            "elligator2 x1",
            "edwards.127",
            "Elligator 2 x1: hash + prepare",
            True,
        ),
        ("weier127", "sswu x1", "weier.127", "SSWU x1: hash + prepare", True),
    ],
)
def test_comparison_recipe_respects_output_representation(
    tmp_path, family, map_name, group_name, recipe, needs_prepare
):
    root = tmp_path / "criterion"
    for group, function, value in [
        ("hash_to_curve", f"{family}/{map_name}", 20),
        ("group.add", f"{group_name}/mode=throughput", 3),
    ]:
        write_bench(
            root,
            group,
            function,
            None,
            value * 0.9,
            value,
            value * 1.1,
            {"Elements": 1},
        )
    raw, _ = br.load(root)
    before = br.pipeline_cost(br.tidy(raw))
    if needs_prepare:
        assert before.empty
    else:
        assert before.iloc[0].recipe == recipe
    write_bench(
        root,
        "group.prepare",
        f"{group_name}/mode=indep",
        None,
        6,
        7,
        8,
        {"Elements": 1},
    )
    raw, _ = br.load(root)
    result = br.pipeline_cost(br.tidy(raw))
    assert len(result) == 1
    r = result.iloc[0]
    assert r.recipe == recipe
    assert r["m=20"] == pytest.approx(
        20 + (7 if needs_prepare else 0) + br.mapping_degree(20) * 3
    )
    assert (pd.notna(r.prepare_ns)) == needs_prepare
    assert map_name in r["hash bench"]


def test_comparison_maps_do_not_borrow_another_models_addition(tmp_path):
    root = tmp_path / "criterion"
    for group, function in [
        ("hash_to_curve", "gf2_127/pornin map x1"),
        ("group.add", "binary-u.127/mode=throughput"),
        ("group.prepare", "binary-u.127/mode=indep"),
    ]:
        write_bench(root, group, function, None, 9, 10, 11, {"Elements": 1})
    raw, _ = br.load(root)
    assert br.pipeline_cost(br.tidy(raw)).empty


@pytest.mark.parametrize(
    "family,suffix,group_name",
    [
        ("gf2_127-u", "(u, v)", "binary-u.127"),
        ("gf2_109-lambda", "(x, λ)", "binary-lambda.109"),
    ],
)
def test_pornin_direct_addend_contract_keeps_its_model(
    tmp_path, family, suffix, group_name
):
    root = tmp_path / "criterion"
    for mode, ending in [("per-element", ""), ("batch", ", batched")]:
        write_bench(
            root,
            "hash_to_curve",
            f"{family}/pornin map x1 to {suffix}{ending}",
            None,
            19,
            20,
            21,
            {"Elements": 1},
        )
    write_bench(
        root,
        "group.add",
        f"{group_name}/mode=throughput",
        None,
        2,
        3,
        4,
        {"Elements": 1},
    )
    raw, _ = br.load(root)
    cost = br.pipeline_cost(br.tidy(raw))
    assert len(cost) == 2
    assert set(cost.family) == {family}
    assert set(cost.recipe) == {"Pornin x1: hash to addend"}
    assert cost.prepare_ns.isna().all()
    assert all(suffix in bench for bench in cost["hash bench"])


def test_insert_grid_preserves_each_map_recipe(rendered):
    _, t = rendered
    p = br.pipeline_cost(t)
    grid = br.insert_grid(p, br.DEFAULT_LOADS, "one at a time")
    rows = [r for r in grid.rows() if r.label == "gf2_127"]
    assert len(rows) == 3
    assert len({r.tips["hash_ns"] for r in rows}) == 3
    mapped = [r for r in rows if "pornin map" in r.tips["hash_ns"]]
    assert len(mapped) == 2
    assert all(r.values["prepare_ns"] == "in the hash" for r in mapped)
    assert all("already an addend" in r.note for r in mapped)


def test_insert_plot_accepts_only_a_direct_addend_recipe(rendered, tmp_path):
    _, t = rendered
    p = br.pipeline_cost(t)
    p = p[(p.family == "gf2_127") & (p.recipe == "Pornin x1: hash to addend")]
    figure = br.plot_insert(p, br.DEFAULT_LOADS, tmp_path, ("svg",))
    assert (tmp_path / "insert.svg").is_file()
    assert "lowest compatible batched recipe" in figure.alt


def test_interval_boxes_preserve_bounds_and_an_estimate_outside_them():
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    ax.set_xscale("log")
    figures.interval_boxes(
        ax, [0, 1], [10, 20], [10.01, 30], [9, 25], [figures.SURFACE, "blue"]
    )
    assert len(ax.patches) == 2
    assert ax.patches[0].get_x() == 10
    assert ax.patches[0].get_width() == pytest.approx(0.01)
    assert ax.patches[1].get_x() == 20
    assert ax.patches[1].get_width() == 10
    assert len(ax.collections) == 2
    assert ax.collections[0].get_segments()[0][:, 0].tolist() == [9, 9]
    assert ax.collections[1].get_segments()[0][:, 0].tolist() == [25, 25]
    assert not ax.lines  # no dot obscures the bounds
    plt.close(fig)


def test_riblt_plan_follows_the_modeled_insertion(tmp_path):
    root = tmp_path / "criterion"
    per = {"Elements": 1024}

    def bench(group, function, ns):
        write_bench(root, group, function, None, ns * 1014, ns * 1024, ns * 1034, per)

    for fam, h, prep, add in [
        ("binary-u.127", 80, 10, 10),
        ("binary.109", 75, 10, 12),
        ("weier.127", 1000, 100, 30),
    ]:
        bench("h2c", f"{fam}/mode=batch,n=1024", h)
        bench("group.prepare", f"{fam}/mode=batch,n=1024", prep)
        bench("group.add", f"{fam}/mode=throughput", add)
    bench("hash_to_curve", "gf2_127-u/pornin map x1 to (u, v), batched", 56)
    listed = [
        f"riblt.encode/{fam}/m=1350,n=3500,h2c={h}"
        for fam, hs in [
            ("binary-u.127", ["ti", "pornin-addend"]),
            ("binary.109", ["ti", "pornin"]),
            ("weier.127", ["ti", "sswu"]),
            ("edwards.127", ["ti", "elligator2"]),
        ]
        for h in hs
    ] + ["riblt.encode/xor-siphash.64/m=1350,n=3500"]

    def plan():
        return br.riblt_plan(br.tidy(br.load(root)[0]), listed)

    # Pornin to (u, v) needs no prepare: 56 + 13.57 * 10 against 80 + 10 +
    # 13.57 * 10; binary.109 is within the leaders' margin, weier.127 far
    # behind. Unmeasured constructions do not compete, and edwards.127, with
    # no model, is left to its spot check.
    assert plan() == {
        "binary-u.127": ("pornin-addend", "full"),
        "binary.109": ("ti", "buffer"),
        "weier.127": ("ti", "spot"),
    }
    # a faster measurement changes the plan; no default does
    bench("hash_to_curve", "gf2_109/pornin map x1, batched", 20)
    assert plan()["binary.109"] == ("pornin", "full")


def test_group_plan_leads_and_contrasts_by_the_modeled_insertion(tmp_path):
    root = tmp_path / "criterion"
    per = {"Elements": 1024}

    def bench(group, function, ns):
        write_bench(root, group, function, None, ns * 1014, ns * 1024, ns * 1034, per)

    def cost(fam, h, add):
        bench("h2c", f"{fam}/mode=batch,n=1024", h)
        bench("group.prepare", f"{fam}/mode=batch,n=1024", 10)
        bench("group.add", f"{fam}/mode=throughput", add)

    cost("binary-u.127", 80, 10)
    cost("binary.109", 75, 12)
    cost("weier.127", 1000, 30)
    cost("edwards.127", 2000, 30)
    listed = [
        f"{g}/{fam}/mode=throughput"
        for g in ("group.add", "group.sub")
        for fam in ("binary-u.127", "binary.109", "weier.127", "edwards.127")
    ]

    def plan():
        return br.group_plan(br.tidy(br.load(root)[0]), listed)

    # the binary families lead; the cheaper odd-characteristic curve contrasts
    assert plan() == {
        "binary-u.127": "lead",
        "binary.109": "lead",
        "weier.127": "contrast",
    }
    # an odd-characteristic leader needs no contrast
    cost("edwards.127", 70, 10)
    assert plan() == {
        "binary-u.127": "lead",
        "binary.109": "lead",
        "edwards.127": "lead",
    }
    # families the suite does not offer are not planned
    assert "binary.109" not in br.group_plan(
        br.tidy(br.load(root)[0]), [i for i in listed if "binary.109" not in i]
    )
