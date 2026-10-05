from types import SimpleNamespace

import pandas as pd
import pytest
import summary_pairs


def run(name, generic=False, rows=None, **meta):
    target = "generic" if generic else "apple-m4"
    return SimpleNamespace(
        name=name,
        meta={
            "cpu": "Apple M4",
            "arch": "arm64",
            "commit": "abc",
            "dirty": False,
            "rustc": "rustc 1.98.1",
            "profile": "quick",
            "target": target,
            "rustflags": f"-C target-cpu={target}",
        }
        | meta,
        table=pd.DataFrame([row()] if rows is None else rows),
    )


def row(**changes):
    return {
        "group": "field/gf2_127",
        "function": "mul throughput (8 chains)",
        "parameter": None,
        "full_id": "field/gf2_127/mul throughput (8 chains)",
        "stat": "slope",
        "family": "gf2_127",
        "layer": "field",
        "operation": "mul",
        "mode": "throughput",
        "unit": "ns/elem",
        "elements": 1024,
        "batch_n": float("nan"),
        "value_ns": 10.0,
        "value_lo_ns": 8.0,
        "value_hi_ns": 12.0,
    } | changes


def test_ratio_direction_bounds_and_provenance():
    native = run("native")
    generic = run(
        "generic", True, [row(value_ns=30.0, value_lo_ns=24.0, value_hi_ns=36.0)]
    )
    result = summary_pairs.compare([generic, native])
    assert len(result) == 1
    pair = result[0]
    assert (pair["value"], pair["lo"], pair["hi"]) == (3.0, 2.0, 4.5)
    assert pair["run"] == "native / generic"
    assert pair["family"] == "gf2_127"
    assert pair["native_id"] == pair["generic_id"] == row()["full_id"]
    assert (pair["native_run"], pair["generic_run"]) == ("native", "generic")


@pytest.mark.parametrize(
    "change",
    [
        {"cpu": "Apple M3"},
        {"arch": "x86_64"},
        {"commit": "def"},
        {"rustc": "rustc 1.97"},
        {"profile": "full"},
        {"dirty": True},
        {"commit": ""},
        {"rustflags": "-C target-cpu=apple-m4"},
        {"rustflags": "-C target-cpu=generic -C opt-level=1"},
    ],
)
def test_rejects_incompatible_metadata(change):
    assert summary_pairs.compare([run("n"), run("g", True, **change)]) == []


def test_rejects_different_recorded_hosts():
    assert (
        summary_pairs.compare([run("n", host="one"), run("g", True, host="two")]) == []
    )


@pytest.mark.parametrize(
    "change",
    [
        {"function": "mul throughput (4 chains)"},
        {"parameter": "different"},
        {"group": "another suite"},
        {"stat": "mean"},
        {"elements": 512},
        {"batch_n": 512},
        {"family": "gf2_109"},
        {"unit": "ns/iter"},
    ],
)
def test_rejects_mismatched_benchmark_identity(change):
    assert summary_pairs.compare([run("n"), run("g", True, [row(**change)])]) == []


def test_missing_or_ambiguous_rows_do_not_pair():
    assert summary_pairs.compare([run("n")]) == []
    assert summary_pairs.compare([run("n"), run("g", True, [])]) == []
    assert summary_pairs.compare([run("n", rows=[row(), row()]), run("g", True)]) == []


def test_selected_operations_and_batch_size():
    records = [
        row(operation="invert", function="invert"),
        row(family="gf2_122", group="field/gf2_122"),
        row(family="edwards127", layer="group ops", group="group.add", operation="add"),
        row(
            family="gf2_122-gls",
            layer="group ops",
            group="group.prepare",
            function="binary.122-gls/mode=batch,n=1024",
            operation="prepare",
            mode="batch",
            batch_n=1024,
        ),
        row(operation="square"),
        row(family="fp107"),
        row(
            layer="group ops",
            group="group.prepare",
            operation="prepare",
            mode="batch",
            batch_n=512,
        ),
    ]
    result = summary_pairs.compare([run("n", rows=records), run("g", True, records)])
    assert len(result) == 4
    assert {r["operation"] for r in result} == {
        "field invert",
        "field mul throughput",
        "point add throughput",
        "point prepare batch 1024",
    }


def test_other_cpu_target_is_not_native():
    n = run("n", target="apple-m3", rustflags="-C target-cpu=apple-m3")
    assert summary_pairs.compare([n, run("g", True)]) == []


def test_explicit_native_target_and_equivalent_flag_spelling():
    n = run("n", target="native", rustflags="-Ctarget-cpu=native -C opt-level=3")
    g = run("g", True, rustflags="-C target-cpu=generic -Copt-level=3")
    assert len(summary_pairs.compare([n, g])) == 1


@pytest.mark.parametrize(
    "key", ["cpu", "arch", "commit", "rustc", "profile", "dirty", "rustflags"]
)
def test_missing_provenance_rejected(key):
    n, g = run("n"), run("g", True)
    del n.meta[key]
    del g.meta[key]
    assert summary_pairs.compare([n, g]) == []


@pytest.mark.parametrize(
    "key,value", [("value_ns", 0), ("value_lo_ns", -1), ("value_hi_ns", float("inf"))]
)
def test_invalid_estimates_rejected(key, value):
    assert (
        summary_pairs.compare([run("n"), run("g", True, [row(**{key: value})])]) == []
    )


def test_unknown_estimator_does_not_create_a_pair():
    n, g = run("n"), run("g", True)
    n.table = n.table.drop(columns="stat")
    g.table = g.table.drop(columns="stat")
    assert summary_pairs.compare([n, g]) == []


@pytest.mark.parametrize("base", ["gf2_109", "gf2_122-gls", "gf2_127"])
@pytest.mark.parametrize("variant", ["", "-lambda", "-u", "-w"])
@pytest.mark.parametrize("operation", ["add", "prepare"])
def test_binary_model_variants_pair_only_with_the_same_model(base, variant, operation):
    family = base + variant
    mode = "throughput" if operation == "add" else "batch"
    function = f"{family}/mode={mode}" + (",n=1024" if mode == "batch" else "")
    native = row(
        family=family,
        layer="group ops",
        group=f"group.{operation}",
        operation=operation,
        mode=mode,
        function=function,
        full_id=f"group.{operation}/{function}",
        batch_n=1024 if mode == "batch" else float("nan"),
    )
    generic = {**native, "value_ns": 30.0, "value_lo_ns": 24.0, "value_hi_ns": 36.0}
    result = summary_pairs.compare([run("n", rows=[native]), run("g", True, [generic])])
    assert len(result) == 1
    assert result[0]["family"] == family
    assert result[0]["native_id"] == result[0]["generic_id"] == native["full_id"]
    assert result[0]["value"] == 3.0
    other_model = {**generic, "family": base + ("-w" if variant != "-w" else "-u")}
    assert (
        summary_pairs.compare([run("n", rows=[native]), run("g", True, [other_model])])
        == []
    )


@pytest.mark.parametrize("family", ["gf2_122-lambda", "gf2_107-u", "gf2_127-unknown"])
def test_unselected_curve_models_remain_excluded(family):
    r = row(family=family, layer="group ops", group="group.add", operation="add")
    assert summary_pairs.compare([run("n", rows=[r]), run("g", True, [r])]) == []


def test_references_pair_their_additions_at_throughput_only():
    def add(family, function, mode):
        return row(
            family=family,
            layer="group ops",
            group="add",
            operation="add",
            mode=mode,
            function=function,
            full_id=f"add/{function}",
        )

    rows = [
        add("ristretto255", "ristretto255/+=, 8 accumulators", "throughput"),
        add("ristretto255", "ristretto255/+=, 1 accumulator", "latency"),
        add("xor", "xor-sha256/xor 32B", "per-element"),
    ]
    result = summary_pairs.compare([run("native", rows=rows), run("g", True, rows)])
    assert [(r["family"], r["operation"]) for r in result] == [
        ("ristretto255", "point add throughput")
    ]
    rows.append(add("xor", "xor-sha256/xor 32B, 8 accumulators", "throughput"))
    result = summary_pairs.compare([run("native", rows=rows), run("g", True, rows)])
    assert ("xor", "XOR add throughput") in {
        (r["family"], r["operation"]) for r in result
    }
