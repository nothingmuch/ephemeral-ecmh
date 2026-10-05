"""Selection stages are compared with measured insertion work, never added twice."""

from types import SimpleNamespace

import bench_report as br
import pandas as pd
import pytest
from summary_selection import references, selection


def bench(group, function, parameter, ns, elements=None):
    return {
        "group": group,
        "function": function,
        "parameter": parameter,
        "full_id": "/".join(
            str(v) for v in (group, function, parameter) if v is not None
        ),
        "estimate_ns": ns,
        "ci_lo_ns": ns / 2,
        "ci_hi_ns": ns * 2,
        "elements": elements,
    }


def run(*rows, name="run"):
    table = br.tidy(br.per_element(pd.DataFrame(rows))) if rows else pd.DataFrame()
    return SimpleNamespace(name=name, meta={}, table=table)


def insertion(family="binary.127", per_item=100, parameter="m=1350,n=3500"):
    return bench("riblt.encode", family, parameter, per_item * 3500, 3500)


def search(stage, seed, ns, index=30):
    return bench(f"curvegen/{stage}", f"gf2_127 agm+sieve/index {index}", str(seed), ns)


def test_certify_already_includes_search_and_uses_seed_points_not_intervals():
    r = run(
        insertion(),
        search("find", 0, 900),
        search("find", 1, 1900),
        search("certify", 0, 1000),
        search("certify", 1, 2000),
    )
    rows = selection([r])
    assert len(rows) == 1
    got = rows[0]
    assert got["run"] == "run" and got["family"] == "gf2_127"
    assert got["method"] == "agm+sieve"
    assert got["stage"] == "search + certificate"
    assert got["complete"] is True
    assert (got["lo_items"], got["hi_items"], got["seeds"]) == (10, 20, 2)
    assert len(got["raw_ids"]) == 2
    assert all("/certify/" in i for i in got["raw_ids"])
    assert got["insert_id"] == "riblt.encode/binary.127/m=1350,n=3500"
    assert got["insert_ns_per_item"] == 100


def test_partial_search_and_point_counts_are_separate_and_labeled():
    r = run(
        insertion(),
        insertion("edwards127", 200),
        insertion("weier127", 400),
        search("find", 0, 1000),
        bench("curvegen/count", "gf2_127", None, 500),
        bench("curvegen/count", "edwards127", None, 1000),
        bench("curvegen/count", "weier127", None, 4000),
    )
    rows = selection([r])
    assert len(rows) == 4
    assert all(not row["complete"] for row in rows)
    assert {row["stage"] for row in rows} == {
        "search only (partial)",
        "one candidate point count (partial)",
    }
    counts = {row["family"]: row for row in rows if row["method"] == "PARI"}
    assert {fam: row["lo_items"] for fam, row in counts.items()} == {
        "gf2_127": 5,
        "edwards127": 5,
        "weier127": 10,
    }
    assert all(row["lo_items"] == row["hi_items"] for row in counts.values())


def test_no_cross_run_or_base_family_substitution():
    r = run(
        insertion("binary.109"),
        insertion("binary-lambda.127"),
        insertion(parameter="m=150,n=3500"),
        search("certify", 0, 1000),
        name="selection-only",
    )
    other = run(insertion(), name="insertion-only")
    assert selection([r, other]) == []


def test_absent_or_unsupported_data_has_no_rows():
    assert selection([]) == []
    assert selection([run()]) == []
    assert selection([run(insertion())]) == []
    assert (
        selection(
            [run(insertion(), bench("curvegen/verify_full", "gf2_127", "0", 100))]
        )
        == []
    )


def test_missing_seed_is_not_backfilled_with_find_or_count():
    r = run(
        insertion(),
        search("find", 0, 900),
        search("find", 1, 1900),
        search("certify", 0, 1000),
    )
    rows = selection([r])
    assert len(rows) == 1
    assert rows[0]["seeds"] == 1
    assert rows[0]["lo_items"] == rows[0]["hi_items"] == 10


@pytest.mark.parametrize("bad", [0, -1, float("nan"), float("inf")])
def test_invalid_insertion_cost_is_not_a_denominator(bad):
    assert selection([run(insertion(per_item=bad), search("certify", 0, 1000))]) == []


def test_ambiguous_insertion_or_duplicate_seed_does_not_choose_a_minimum():
    assert (
        selection(
            [run(insertion(), insertion(per_item=90), search("certify", 0, 1000))]
        )
        == []
    )
    assert (
        selection(
            [
                run(
                    insertion(),
                    search("certify", 0, 1000),
                    search("certify", 0, 2000, index=54),
                )
            ]
        )
        == []
    )


def test_mapping_variant_is_not_a_second_insertion():
    mapped = insertion(per_item=90, parameter="m=1350,n=3500,map=mcg64")
    rows = selection([run(insertion(), mapped, search("certify", 0, 1000))])
    assert [r["insert_ns_per_item"] for r in rows] == [100]


def test_malformed_certify_group_is_not_replaced_with_search_only():
    r = run(insertion(), search("find", 0, 900), search("certify", 0, 1000))
    r.table.loc[r.table.group == "curvegen/certify", "function"] = float("nan")
    assert selection([r]) == []


def test_invalid_seed_measurement_does_not_silently_shrink_seed_range():
    r = run(insertion(), search("certify", 0, 1000), search("certify", 1, float("nan")))
    assert selection([r]) == []


def test_insertion_is_timed_under_the_cheapest_hash():
    ti = insertion(per_item=100, parameter="m=1350,n=3500,h2c=ti")
    pornin = insertion(per_item=90, parameter="m=1350,n=3500,h2c=pornin")
    rows = selection([run(ti, pornin, search("certify", 0, 1000))])
    assert [r["insert_ns_per_item"] for r in rows] == [90]


def test_references_are_the_baselines_insertion_alone():
    r = run(
        insertion(),
        insertion("xor-sha256.64", 50),
        insertion("ristretto255", 5000),
        search("certify", 0, 1000),
    )
    assert [(x["family"], x["insert_ns_per_item"]) for x in references([r])] == [
        ("xor", 50),
        ("ristretto255", 5000),
    ]
    assert references([run(insertion())]) == []


def selected(monkeypatch):
    import summary as sm

    r = run(
        insertion(),
        insertion("weier127", 400),
        insertion("xor-sha256.64", 50),
        insertion("ristretto255", 5000),
        search("certify", 0, 1000),
        search("certify", 1, 4000),
        bench("curvegen/count", "weier127", None, 8000),
    )
    figs = []
    monkeypatch.setattr(sm, "_save", lambda fig, out, name: figs.append(fig) or name)
    return sm, selection([r]), references([r]), figs


def test_cost_per_item_is_insertion_plus_amortised_selection(tmp_path, monkeypatch):
    sm, sel, refs, figs = selected(monkeypatch)
    assert sm.per_item(sel, refs, tmp_path) == "selection-per-item"
    ax = figs[0].axes[0]
    gf2 = next(r for r in sel if r["family"] == "gf2_127")
    fastest = ax.get_lines()[0]
    x, y = fastest.get_xdata()[0], fastest.get_ydata()[0]
    assert y == pytest.approx(100 * (1 + gf2["lo_items"] / x))
    flat = [line for line in ax.get_lines() if len(set(line.get_ydata())) == 1]
    assert sorted(line.get_ydata()[0] for line in flat) == [50, 5000]
    legend = [t.get_text().replace("\n", " ") for t in ax.get_legend().get_texts()]
    assert any(t.startswith("xor (baseline): insertion alone") for t in legend)


def test_selection_time_is_absolute_and_states_the_baselines(tmp_path, monkeypatch):
    sm, sel, refs, figs = selected(monkeypatch)
    assert sm.selection_times(sel, refs, tmp_path) == "selection-time"
    ax = figs[0].axes[0]
    ticks = [t.get_text() for t in ax.get_yticklabels()]
    assert ticks[:2] == ["xor (baseline)", "ristretto255"]
    assert {t.get_text() for t in ax.texts} == {
        "none: no curve to select",
        "none: a fixed group",
    }
    estimates = sorted(
        x for c in ax.collections for seg in c.get_segments() for x in {seg[0][0]}
    )
    assert estimates == [1000, 4000, 8000]


def test_normalised_selection_has_no_baseline_lines(tmp_path, monkeypatch):
    sm, sel, _, figs = selected(monkeypatch)
    sm.amortization(sel, tmp_path)
    legend = [t.get_text() for t in figs[0].axes[0].get_legend().get_texts()]
    assert not any("baseline" in t or "ristretto255" in t for t in legend)
