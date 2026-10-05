import io
import json
import math
import re
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
    "digest",
]

# the workload, and the smaller-field variants the rules must name
EXTRA = [
    # one family's pipeline, through the group traits
    # a λ family's hash straight to its addend
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
        (
            "field/gf2/batch invert (product tree)",
            ("field", "gf2_127", "batch invert", "batch"),
        ),
        ("field/gf2/halftrace", ("field", "gf2_127", "halftrace", "per-element")),
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
        # the map straight to an addend: a comparison map that is also a prepare
        (
            "hash_to_curve/secp256k1/ellswift decode",
            ("hash to curve", "secp256k1", "hash to curve", "per-element"),
        ),
        (
            "hash_to_curve/sha256 (XOR baseline)",
            ("hash to curve", "xor", "hash to curve", "per-element"),
        ),
        # smaller fields: gf2_109, edwards107, weier107; 127 is the base
        # GF(2^122): qsolve sits with the halftraces
        # a batch to affine, not a sum
        # benches/group.rs: layered ids, the family spelled <name>.<bits>
        # Identifiers with an explicit 127-bit field width.
        (
            "field/gf2_127/mul latency (dependent chain)",
            ("field", "gf2_127", "mul", "latency"),
        ),
        # binary127::lambda: gf2_127's curves, λ-projective accumulators
        # riblt.peel's parameters say whether it peels in batches
        # the two XOR baselines are told apart by their hash
        # binary122: GLS constants, λ accumulators, and both
        # the F_{p^2} prototypes, and the codecs of the odd fields
        # the Weierstrass curves' Jacobian families
        # Plonky3's fields
        (
            "field/gf2_127/normalize (to_u128)",
            ("field", "gf2_127", "normalize", "per-element"),
        ),
        # binary109's λ and w families, and the λ one's hash to its addend
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
    [],
)
def test_both_spellings_name_one_family(old, new):
    assert rules.family(old) == rules.family(new)
    assert rules.family(new)[0] == rules.family(new)[1]  # the base, not a variant


@pytest.mark.parametrize(
    "group, function, parameter, operation",
    [],
)
def test_curvegen_names_each_family(group, function, parameter, operation):
    f = rules.classify(group, function, parameter, False)
    fam = function.split("/")[0].split(" ")[0]
    assert (f.layer, f.family, f.operation) == ("curve generation", fam, operation)


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


IDS = Path(__file__).parent / "fixtures" / "ids-2026-10-03.tsv"


def ids() -> list[list[str]]:
    lines = IDS.read_text().splitlines()
    return [line.split("\t") for line in lines if not line.startswith("#")]


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


def test_svg_is_deterministic(rendered, tmp_path):
    out, t = rendered
    br.plot_layer(t, "digest", tmp_path, ("svg",))
    assert (tmp_path / "digest.svg").read_bytes() == (out / "digest.svg").read_bytes()


def test_empty_tree_is_an_error(tmp_path):
    with pytest.raises(SystemExit):
        br.report(tmp_path, tmp_path / "out")


def test_no_machine_is_flagged(rendered):
    md = (rendered[0] / "report.md").read_text()
    assert "## Machine\n\nUnknown: no meta.json" in md


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


def tables_md(g):
    return __import__("tables").to_markdown(g, br.fmt_time, br._md_cell)


# per-element ns: add, hash, prepare, encode, decode (None: not measured)


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
}


def write_run(run: Path) -> Path:
    root = write_tree(ROWS, run / "criterion")
    for g, f, p, n in EXTRA:
        write_bench(root, g, f, p, 90.0 * n, 100.0 * n, 110.0 * n, {"Elements": n})
    (run / "meta.json").write_text(json.dumps(RUN_META))
    return run


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    return write_run(tmp_path_factory.mktemp("bench-runs") / "fixture-run")


def test_export_keeps_the_raw_estimates(run, tmp_path):
    dest = tmp_path / "results" / "fixture-run"
    br.main(["--export", str(run), str(dest)])
    assert sorted(p.name for p in dest.iterdir()) == [
        "benchmarks.csv",
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


def test_compact_form_holds_one_statistic(run, tmp_path):
    dest = tmp_path / "fixture-run"
    br.export(run, dest)
    with pytest.raises(SystemExit, match="typical estimate only"):
        br.report(dest, tmp_path / "out", stat="median")


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
