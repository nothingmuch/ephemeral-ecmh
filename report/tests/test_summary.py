import pandas as pd
import summary as sm


def row(family, value, lo=None, hi=None, **kw):
    return dict(
        family=family,
        value_ns=value,
        value_lo_ns=lo if lo is not None else value * 0.9,
        value_hi_ns=hi if hi is not None else value * 1.1,
        unit="ns/elem",
        group="riblt.encode",
        function=family,
        parameter="m=1350,n=3500",
        full_id=f"riblt.encode/{family}/m=1350,n=3500",
        mode="per-element",
        elements=3500,
        stat="slope",
        confidence=0.95,
        base="gf2_127",
        layer="RIBLT workload",
        operation="encode",
        batch_n=float("nan"),
        **kw,
    )


def run(name, rows):
    t = pd.DataFrame(rows)
    t["family"] = pd.Categorical(t.family)
    t["partial_batch"] = False
    return sm.Run(
        name, {"commit": "abc", "dirty": False, "cpu": "test", "arch": "arm64"}, t
    )


def test_subset_retains_boundary_overlap_and_fixed_representatives():
    d = pd.DataFrame(
        [
            row("a", 10, 9, 11),
            row("b", 20, 19, 21),
            row("c", 30, 29, 32),
            row("d", 33, 31, 34),
            row("e", 100),
            row("edwards127", 200),
            row("xor", 1),
        ]
    )
    assert sm.subset(d, ["edwards127"]) == ["a", "b", "c", "d", "edwards127"]


def test_insert_selection_requires_same_workload_and_keeps_union():
    a = run("a", [row("gf2_127", 10), row("edwards127", 20), row("xor", 1)])
    b = run("b", [row("gf2_109", 9), row("weier127", 21)])
    wrong = row("gf2_122", 0.1)
    wrong["parameter"] = "m=150,n=3500"
    wrong["full_id"] = "riblt.encode/x/m=150,n=3500"
    b.table = pd.concat([b.table, pd.DataFrame([wrong])], ignore_index=True)
    data, families = sm.insertion([a, b])
    assert set(families) == {"gf2_127", "gf2_109", "edwards127", "weier127"}
    assert not any(r["family"] == "gf2_122" for r in data)
    assert next(r for r in data if r["family"] == "xor")["value"] == 1


def test_batch_pairs_preserve_interval_endpoints_and_missing_modes():
    a = row("gf2_127", 100, 90, 110)
    a.update(group="h2c", mode="per-element", full_id="h2c/a/indep", elements=1024)
    b = {
        **a,
        "mode": "batch",
        "full_id": "h2c/a/n=1024",
        "batch_n": 1024,
        "value_ns": 10,
        "value_lo_ns": 8,
        "value_hi_ns": 12,
    }
    small = {**b, "batch_n": 8, "value_ns": 1}
    lone = {**a, "family": "edwards127", "full_id": "h2c/edwards127/indep"}
    data = sm.batching([run("r", [a, b, small, lone])])
    z = next(r for r in data if r["family"] == "gf2_127")
    assert (z["value"], z["lo"], z["hi"]) == (10, 90 / 12, 110 / 8)
    assert z["saved_percent"] == 90
    assert z["saved_lo_percent"] == 100 * (1 - 12 / 90)
    assert z["saved_hi_percent"] == 100 * (1 - 8 / 110)
    missing = next(r for r in data if r["family"] == "edwards127")
    assert missing["batch"] is None and missing["value"] is None


def test_no_runs_has_no_summary(tmp_path):
    assert sm.chapter([], tmp_path / "summary") is None
    assert not (tmp_path / "summary").exists()


def test_sparse_run_chapter_states_missing_views_and_has_light_dark(tmp_path):
    r = run("one", [row("gf2_127", 10), row("xor", 1)])
    page = sm.chapter([r], tmp_path / "summary")
    assert "Summary across runs" in page
    assert "not measured" in page.lower()
    assert "fig-light" in page and "fig-dark" in page
    assert (tmp_path / "summary" / "insertion.svg").is_file()
    assert (tmp_path / "summary" / "insertion-dark.svg").is_file()
    assert (tmp_path / "summary" / "evidence.csv").is_file()


def test_book_places_summary_before_runs(tmp_path, monkeypatch):
    import book
    from test_book import repository

    root = repository(tmp_path / "repo", ["a"])

    def chapter(runs, out):
        assert len(runs) == 1
        out.mkdir(parents=True)
        return "# Summary across runs\n"

    monkeypatch.setattr(sm, "chapter", chapter)
    book.assemble(root, tmp_path / "book", "https://example.org/r", "abc")
    contents = (tmp_path / "book" / "src" / "SUMMARY.md").read_text()
    assert contents.index("[Summary across runs]") < contents.index("[Benchmark run a]")


def test_components_never_treat_missing_prepare_as_zero():
    add = row("gf2_127", 10)
    add.update(group="group.add", mode="throughput", operation="add", elements=1024)
    h = {
        **add,
        "group": "h2c",
        "mode": "batch",
        "batch_n": 1024,
        "operation": "hash to curve",
        "value_ns": 50,
    }
    data, _ = sm.components([run("r", [add, h])])
    assert data == []


def test_three_run_panels_have_no_text_overlap(tmp_path):
    runs = [
        run(name, [row("gf2_127", 10), row("edwards127", 20), row("weier127", 25)])
        for name in ["m4-native", "m4-generic", "x86-native"]
    ]
    data, families = sm.insertion(runs)
    assert sm.dots(data, tmp_path, "insertion", families) == "insertion.svg"


def test_batching_keeps_legacy_references_as_individual_only():
    a = row("ristretto255", 5000)
    a.update(
        group="hash_to_curve",
        function="ristretto255/hash_from_bytes<Sha512>",
        full_id="hash_to_curve/ristretto255/hash_from_bytes<Sha512>",
    )
    b = row("secp256k1", 6000)
    b.update(
        group="hash_to_curve",
        function="secp256k1/ellswift decode",
        full_id="hash_to_curve/secp256k1/ellswift decode",
    )
    wrong = {**b, "function": "secp256k1/try-and-increment (parse)", "value_ns": 1}
    data = sm.batching([run("r", [a, b, wrong])])
    assert [r["individual"] for r in data] == [5000, 6000]
    assert all(r["batch"] is None for r in data)


def test_components_choose_1024_even_when_a_larger_batch_was_measured():
    add = row("gf2_127", 10)
    add.update(group="group.add", mode="throughput", operation="add", elements=1024)
    h = row("gf2_127", 50)
    h.update(
        group="h2c",
        mode="batch",
        batch_n=1024,
        operation="hash to curve",
        elements=1024,
    )
    prep = row("gf2_127", 5)
    prep.update(
        group="group.prepare",
        mode="batch",
        batch_n=1024,
        operation="prepare",
        elements=1024,
    )
    big = {**h, "batch_n": 4096, "value_ns": 1}
    r = run("r", [add, h, prep, big])
    r.table.loc[r.table.batch_n == 1024, "partial_batch"] = True
    data, _ = sm.components([r])
    assert len(data) == 1 and data[0]["hash_prepare"] == 55


def test_alt_describes_the_figure_without_repeating_caption():
    from html.parser import HTMLParser

    class Images(HTMLParser):
        def __init__(self):
            super().__init__()
            self.alt = []

        def handle_starttag(self, tag, attrs):
            if tag == "img":
                self.alt.append(dict(attrs)["alt"])

    caption = "The visible caption contains the uncertainty qualifications."
    page = sm._figure("insertion.svg", caption)
    parser = Images()
    parser.feed(page)
    assert len(parser.alt) == 2 and parser.alt[0] == parser.alt[1]
    assert "RIBLT" in parser.alt[0]
    assert caption not in parser.alt[0]
    assert page.count(caption) == 1


def test_mapping_variants_are_not_workload_measurements():
    mapped = row("gf2_127", 5)
    mapped.update(
        parameter="m=1350,n=3500,map=mcg64",
        full_id="riblt.encode/gf2_127/m=1350,n=3500,map=mcg64",
    )
    r = run("r", [row("gf2_127", 10), row("edwards127", 20), mapped])
    data, _ = sm.insertion([r])
    assert {(d["family"], d["value"]) for d in data} == {
        ("gf2_127", 10),
        ("edwards127", 20),
    }
    encode = sm.scaling([r], ["gf2_127"])
    assert [(d["family"], d["value"]) for d in encode] == [("gf2_127", 10)]


def test_scaling_lines_keep_each_size_and_run(tmp_path):
    def at(family, m, value):
        r = row(family, value)
        r.update(
            parameter=f"m={m},n=3500", full_id=f"riblt.encode/{family}/m={m},n=3500"
        )
        return r

    def hashed(family, n, value):
        r = row(family, value)
        r.update(
            group="h2c",
            mode="batch" if n > 1 else "per-element",
            batch_n=float(n) if n > 1 else float("nan"),
            parameter="",
            full_id=f"h2c/{family}/{n}",
        )
        return r

    rows = [
        at(f, m, v * m**0.1)
        for f, v in [("gf2_127", 10), ("gf2_109", 9)]
        for m in [5, 150, 1350]
    ]
    rows += [hashed("gf2_127", n, 100 / n**0.5) for n in [1, 8, 1024]]

    def mapped(function, batched, value):
        r = row("gf2_127", value)
        r.update(
            group="hash_to_curve",
            layer="comparison maps",
            operation="two maps summed" if "x2" in function else "one map",
            function=f"gf2_127/{function}" + (", batched" if batched else ""),
            mode="batch" if batched else "per-element",
            elements=1024,
            parameter="",
            full_id=f"hash_to_curve/gf2_127/{function}/{batched}",
        )
        return r

    rows += [
        mapped("pornin map x1", False, 80),
        mapped("pornin map x1", True, 12),
        mapped("pornin map x2", False, 160),
    ]
    runs = [run(name, rows) for name in ["m4-native", "x86-native"]]
    data = sm.scaling(runs, ["gf2_127", "gf2_109"])
    encode = [r for r in data if r["series"] == sm.SCALING[0][0]]
    assert len(encode) == 12
    assert sorted({r["x"] for r in encode}) == [5, 150, 1350]
    hashes = [r for r in data if r["series"] == sm.SCALING[2][0]]
    m4 = [r for r in hashes if r["run"] == "m4-native"]
    assert sorted(r["x"] for r in m4 if r["map"] == "try-and-increment") == [
        1,
        8,
        1024,
    ]
    assert sorted((r["x"], r["value"]) for r in m4 if r["map"] == "Pornin's map") == [
        (1, 80),
        (1024, 12),
    ]
    assert len(m4) == 5
    assert sm.curves(data, ["gf2_127", "gf2_109"], tmp_path) == "scaling.svg"
    assert (tmp_path / "scaling-dark.svg").is_file()


def test_each_line_has_its_own_colour_and_bases_keep_theirs():
    colors = sm._colors(["gf2_127-u", "gf2_109-u", "gf2_127", "edwards127", "xor"])
    assert colors["gf2_127-u"] == sm.figures.color("gf2_127")
    assert colors["edwards127"] == sm.figures.color("edwards127")
    assert len(set(colors.values())) == 5


def test_selection_beyond_the_spare_hues_has_distinct_colours():
    # seven of the twelve share a base with an earlier one, and six hues
    # are spare
    families = [
        "gf2_109-u",
        "gf2_127-u",
        "gf2_109-lambda",
        "gf2_109-w",
        "gf2_109",
        "gf2_127",
        "gf2_122-gls",
        "edwards127",
        "weier127",
        "ristretto255",
        "weier127-jacobian",
        "xor",
    ]
    colors = sm._colors(families)
    assert colors["gf2_109-u"] == sm.figures.color("gf2_127")
    assert colors["weier127"] == sm.figures.color("weier127")
    assert len(set(colors.values())) == len(families)
    assert len(sm._colors([f"f{i}" for i in range(40)])) == 40


def test_insertion_findings_identify_overlaps_and_separation():
    r = run(
        "r",
        [
            row("a", 10, 9, 11),
            row("b", 11, 10, 12),
            row("c", 20, 19, 21),
            row("d", 21, 20, 22),
        ],
    )
    data, _ = sm.insertion([r])
    finding = sm.insertion_findings(data)[0]
    assert "first two recorded intervals overlap" in finding
    assert "100.0% more time" in finding
    assert "Also overlapping the third interval: d" in finding


def test_layout_warning_does_not_discard_measurements(tmp_path, monkeypatch):
    import pytest

    monkeypatch.setattr(sm, "_overlaps", lambda fig: [("long label", "long title")])
    data, families = sm.insertion([run("r", [row("a", 10)])])
    with pytest.warns(UserWarning, match="overlapping figure text"):
        sm.dots(data, tmp_path, "insertion", families)
    assert (tmp_path / "insertion.svg").is_file()


def hashed(family, function, mode, value, group="hash_to_curve", **kw):
    r = row(family, value)
    r.update(
        group=group,
        layer="hash to curve",
        operation="hash to curve",
        function=function,
        parameter="",
        full_id=f"{group}/{function}/{mode}",
        mode=mode,
        elements=1024,
        batch_n=1024.0 if mode == "batch" else float("nan"),
        **kw,
    )
    return r


BASELINE_HASHES = [
    hashed("xor", "sha256 (XOR baseline)", "per-element", 20),
    hashed("ristretto255", "ristretto255/hash_from_bytes<Sha512>", "per-element", 5000),
]


def test_baseline_label_names_only_the_xor_baseline():
    assert sm._label("xor") == "xor (baseline)"
    assert sm._label("xor: SHA-256 digest, no curve").startswith("xor (baseline):")
    assert sm._label("xor-siphash") == "xor-siphash"
    assert sm._label("ristretto255") == "ristretto255"


def test_baseline_hashes_are_named_and_one_at_a_time_only():
    t_and_i = [
        hashed("gf2_127", "gf2_127/t&i", mode, v, group="h2c")
        for mode, v in [("per-element", 800), ("batch", 80)]
    ]
    r = run("r", [*t_and_i, *BASELINE_HASHES])
    bat = {d["family"]: d for d in sm.batching([r])}
    assert bat["xor"]["map"] == "SHA-256 digest, no curve"
    assert bat["ristretto255"]["map"] == "upstream map"
    assert all(bat[f]["batch"] is None for f in ("xor", "ristretto255"))
    labels = [d["label"] for d in sm.map_choices([r])]
    assert labels[0] == "xor: SHA-256 digest, no curve"
    assert "ristretto255: upstream map" in labels
    assert not any(label.startswith("ristretto255: try") for label in labels)


def test_hash_charts_join_a_constructions_two_modes(tmp_path, monkeypatch):
    guides = []
    monkeypatch.setattr(sm, "_guide", lambda ax, y, a, b: guides.append((y, a, b)))
    t_and_i = [
        hashed("gf2_127", "gf2_127/t&i", mode, v, group="h2c")
        for mode, v in [("per-element", 800), ("batch", 80)]
    ]
    r = run("r", [*t_and_i, *BASELINE_HASHES])
    assert sm.dots(sm.map_choices([r]), tmp_path, "map-choice") == "map-choice.svg"
    # the baselines are rows of their own, with no batch to join
    assert guides == [(1, 800, 80)]
    guides.clear()
    assert sm.paired(sm.batching([r]), tmp_path) == "batching.svg"
    assert guides == [(1, 800, 80)]


def test_insertion_draws_the_baseline_as_the_first_row(tmp_path, monkeypatch):
    figs = []
    monkeypatch.setattr(sm, "_save", lambda fig, out, name: figs.append(fig) or name)
    data, families = sm.insertion(
        [run("r", [row("gf2_127", 10), row("ristretto255", 50), row("xor", 1)])]
    )
    sm.dots(data, tmp_path, "insertion", families)
    ax = figs[0].axes[0]
    assert ax.get_yticklabels()[0].get_text() == "xor (baseline)"
    assert len(ax.patches) == 3


def test_scaling_hashes_baselines_at_batch_size_one():
    r = run("r", BASELINE_HASHES)
    hashes = sm.scaling([r], ["xor", "ristretto255"])
    assert {(d["family"], d["x"], d["map"]) for d in hashes} == {
        ("xor", 1, None),
        ("ristretto255", 1, None),
    }


def test_components_take_baselines_one_at_a_time_once_their_add_is_throughput():
    def added(family, mode, value):
        r = row(family, value)
        r.update(
            group="add",
            layer="group ops",
            operation="add",
            mode=mode,
            elements=1024,
            full_id=f"add/{family}/{mode}",
        )
        return r

    curve = row("gf2_127", 10)
    curve.update(group="group.add", mode="throughput", operation="add", elements=1024)
    h = hashed("gf2_127", "gf2_127/t&i", "batch", 50, group="h2c")
    prep = row("gf2_127", 5)
    prep.update(
        group="group.prepare",
        mode="batch",
        batch_n=1024,
        operation="prepare",
        elements=1024,
    )
    rows = [curve, h, prep, *BASELINE_HASHES, added("ristretto255", "throughput", 60)]
    one = added("xor", "per-element", 0.4)
    data, families = sm.components([run("r", [*rows, one])])
    assert families == ["gf2_127", "ristretto255"]
    ristretto = next(d for d in data if d["family"] == "ristretto255")
    assert ristretto["hashing"] == "one at a time" and ristretto["hash_prepare"] == 5000
    data, families = sm.components(
        [run("r", [*rows, one, added("xor", "throughput", 0.1)])]
    )
    assert families[0] == "xor"
    assert sm.stacked_log(data)
