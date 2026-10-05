"""Criterion results -> tidy table, facet-grid plots and a markdown/HTML report.

  nix run .#bench-report -- target/criterion out/
  nix run .#bench-report -- --export bench-runs/ID results/ID
  nix run .#bench-report -- --index runs/
  BENCH --list | sed 's/: benchmark$//' | nix run .#bench-report -- --check-ids

Reads criterion's <group>/<function>[/<parameter>]/new/{benchmark,estimates}.json,
or in its place a run's compact form (--export), maps each id to facet axes
with the rules in rules.py and writes, into the output directory: tidy.csv, elementary.csv (the operations-by-curve
tables), riblt_cost.csv (each family's pipeline), riblt_cost_lower_bound.csv,
coverage.csv (which family has each hot-path operation measured), a PNG + SVG each for the per-item insert cost, the elementary field and group
operations and every suite layer, report.md and report.html. With --index,
it lists the reports under a directory of runs in its index.html and
index.md instead. With --check-ids, it reads benchmark ids and fails if any
matches no rule, which bench-run does before it times anything. With
--export, it writes a run's compact form (export) instead of a report.
"""

import argparse
import html
import json
import math
import re
import shutil
import sys
from pathlib import Path

import figures
import pandas as pd
import rules
import tables
from figures import (
    GRID,
    SURFACE,
    TEXT,
    TEXT_2,
    dark_twin,
    fmt_time,
    plot_elementary,
    plot_layer,
    plot_selection,
)
from tables import PLAIN, Col, Grid, Row

# the table sizes benches/riblt.rs encodes into (its MS)


def throughput(t) -> tuple[int | None, float | None]:
    """(elements, bytes) per iteration from criterion's serde Throughput."""
    if not isinstance(t, dict) or len(t) != 1:
        return None, None
    ((kind, v),) = t.items()
    if kind == "Elements":
        return v, None
    if kind == "ElementsAndBytes":
        return v["elements"], v["bytes"]
    if kind in ("Bytes", "BytesDecimal"):
        return None, v
    if kind == "Bits":
        return None, v / 8
    return None, None


def estimate(e: dict, stat: str) -> tuple[str, dict]:
    """Resolve "typical" to the slope estimate when available, otherwise the mean."""
    if stat == "typical":
        stat = "slope" if e.get("slope") else "mean"
    if not e.get(stat):
        raise KeyError(f"no {stat} estimate")
    return stat, e[stat]


# What a published run keeps of each benchmark (results/<run>/benchmarks.csv,
# written by export): criterion's ids and its typical estimate with its
# interval and throughput. per_element and the rules derive the rest, and
# criterion's full id is its group, function and parameter joined by "/".
RAW = [
    "group",
    "function",
    "parameter",
    "stat",
    "estimate_ns",
    "ci_lo_ns",
    "ci_hi_ns",
    "confidence",
    "elements",
    "bytes",
]
COMPACT = "benchmarks.csv"


def _row(bj: Path, stat: str) -> dict:
    b = json.loads(bj.read_text())
    which, est = estimate(json.loads((bj.parent / "estimates.json").read_text()), stat)
    elements, nbytes = throughput(b.get("throughput"))
    ci = est["confidence_interval"]
    # per_element divides the whole table, where a value that is not a
    # number would fail every row rather than skip this one
    numbers = [est["point_estimate"], ci["lower_bound"], ci["upper_bound"]]
    if not all(isinstance(v, (int, float)) for v in [*numbers, elements or 0]):
        raise TypeError(f"not a number among {numbers + [elements]}")
    return {
        "group": b["group_id"],
        "function": b.get("function_id"),
        "parameter": b.get("value_str"),
        "full_id": b["full_id"],
        "stat": which,
        "estimate_ns": est["point_estimate"],
        "ci_lo_ns": ci["lower_bound"],
        "ci_hi_ns": ci["upper_bound"],
        "confidence": ci["confidence_level"],
        "elements": elements,
        "bytes": nbytes,
    }


def per_element(df: pd.DataFrame) -> pd.DataFrame:
    """Each estimate and its interval per element where the benchmark
    declares an element count, else per iteration."""
    n = pd.to_numeric(df.elements)
    has = n.fillna(0) != 0
    per = n.where(has, 1)
    return df.assign(
        unit=has.map({True: "ns/elem", False: "ns/iter"}),
        value_ns=df.estimate_ns / per,
        value_lo_ns=df.ci_lo_ns / per,
        value_hi_ns=df.ci_hi_ns / per,
    )


def load(root: Path, stat: str = "typical") -> tuple[pd.DataFrame, list[str]]:
    """One row per benchmark, and the directories that could not be read,
    from criterion's tree or a run's compact form."""
    if (root / COMPACT).is_file():
        return load_compact(root / COMPACT, stat), []
    rows, skipped = [], []
    for bj in sorted(root.glob("**/new/benchmark.json")):
        try:
            rows.append(_row(bj, stat))
        except (OSError, ValueError, KeyError, TypeError, AttributeError) as err:
            # Record malformed input while retaining the remaining benchmarks,
            # by paths under root: the report names no place on this machine
            err = repr(err).replace(f"{root}/", "")
            skipped.append(f"{bj.parent.parent.relative_to(root)}: {err}")
    return (per_element(pd.DataFrame(rows)) if rows else pd.DataFrame()), skipped


def load_compact(path: Path, stat: str) -> pd.DataFrame:
    """The table load() reads from criterion's tree, from export's
    benchmarks.csv: equal to it, so that both render the same report."""
    if stat != "typical":
        raise SystemExit(f"{path} holds the typical estimate only, not --stat {stat}")
    ids = {c: str for c in ("group", "function", "parameter", "stat")}
    # an empty field is a missing id or throughput; no id is read as a number
    # or a NaN spelling, and every estimate as the float it was written from
    df = pd.read_csv(
        path,
        dtype=ids,
        keep_default_na=False,
        na_values=[""],
        float_precision="round_trip",
    )
    parts = df[["group", "function", "parameter"]].astype(object)
    df.insert(
        3,
        "full_id",
        ["/".join(p for p in r if isinstance(p, str)) for r in parts.itertuples(False)],
    )
    return per_element(df)


def tidy(df: pd.DataFrame) -> pd.DataFrame:
    """Add the facet axes (layer, family, operation, mode) and a bar label."""
    facets = [
        rules.classify(
            r.group, _str(r.function), _str(r.parameter), pd.notna(r.elements)
        )
        for r in df.itertuples()
    ]
    t = pd.concat([df, pd.DataFrame(facets, columns=rules.Facets._fields)], axis=1)
    t["label"] = [
        "/".join(p for p in (f, v) if isinstance(p, str)) or g
        for g, f, v in zip(t.group, t.function, t.parameter)
    ]
    t["layer"] = pd.Categorical(t.layer, rules.LAYER_ORDER, ordered=True)
    fams = sorted(
        t.family.unique(),
        key=lambda f: (rules.FAMILY_ORDER.index(base_of(t, f)), f),
    )
    t["family"] = pd.Categorical(t.family, fams, ordered=True)
    # facets in rule order, then cheapest first (a name may recur across layers)
    by = t.groupby(["layer", "operation"], observed=True)
    t["op_rank"] = by.op_rank.transform("min")
    # benches/group.rs times each batched step at several batch sizes; the
    # largest stands for "batched", the smaller ones only for the sweep
    t["batch_n"] = t.full_id.str.extract(r"\bn=(\d+)", expand=False).astype(float)
    largest = t.groupby(["group", "family"], observed=True).batch_n.transform("max")
    t["partial_batch"] = t.batch_n < largest
    t["op_cost"] = by.value_ns.transform("median")
    keys = ["layer", "op_rank", "op_cost", "operation", "family", "label"]
    return (
        t.sort_values(keys, kind="stable")
        .drop(columns="op_cost")
        .reset_index(drop=True)
    )


def _str(v) -> str | None:
    return v if isinstance(v, str) else None


def base_of(t: pd.DataFrame, family: str) -> str:
    base = t.loc[t.family == family, "base"].iloc[0]
    return base if base in rules.COLORS else rules.OTHER


def unclassified(t: pd.DataFrame) -> pd.DataFrame:
    """Rows without a layer, a family, or an operation rule."""
    unclaimed = t.op_rank == len(rules.OPERATIONS)
    return t[(t.layer == rules.OTHER) | (t.family == rules.OTHER) | unclaimed]


def classifiable(full_id: str) -> bool:
    """Whether a listed id has a layer, a family and an operation rule, as
    `unclassified` asks of a run's rows. criterion's --list joins group,
    function and parameter with '/', which a group may contain too, so any
    split of the group from the rest will do."""
    parts = full_id.split("/")
    for i in range(1, len(parts)):
        f = rules.classify("/".join(parts[:i]), "/".join(parts[i:]), None, True)
        if rules.OTHER not in (f.layer, f.family) and f.op_rank < len(rules.OPERATIONS):
            return True
    return False


def check_ids(lines) -> list[str]:
    """The ids, one per line, that `classifiable` refuses."""
    return [i for i in (ln.strip() for ln in lines) if i and not classifiable(i)]


# The point stated in docs/workload.md (Cost in repeated reconciliation): a
# namespace of n items, d of them new
# each round, reconciled through m = 1.35 d coded symbols
# the costs dominance compares; the regimes' estimates are sums of them


def selection_runs(root: Path) -> pd.DataFrame:
    """curvegen.csv's rows, one per family, method and seed, under the
    canonical family names; empty without curvegen.csv."""
    found = beside(root, "curvegen.csv")
    if found is None:
        return pd.DataFrame(columns=["family", "method", "seed", "find_s", "verify_s"])
    c = pd.read_csv(found)
    c["family"] = c.family.map(lambda name: rules.SPELLINGS.get(name, name))
    return c


def fastest_method(runs: pd.DataFrame) -> pd.DataFrame:
    """Of runs, those of each family's method that finds fastest on average."""
    if runs.empty:
        return runs
    mean = runs.groupby(["family", "method"]).find_s.mean()
    best = set(mean.groupby(level="family").idxmin())
    return runs[[k in best for k in zip(runs.family, runs.method)]]


def selection(root: Path) -> pd.DataFrame:
    """Mean wall-clock seconds to find a curve and to verify its
    certificate, per curve, over curvegen.csv's seeds, by the method that
    finds fastest; empty without curvegen.csv."""
    return (
        fastest_method(selection_runs(root))
        .groupby(["family", "method"])
        .agg(
            seeds=("seed", "count"),
            find_s=("find_s", "mean"),
            verify_s=("verify_s", "mean"),
        )
        .reset_index()
    )


def elementary(t: pd.DataFrame) -> pd.DataFrame:
    """One row per (elementary operation, field or curve): the cheapest
    bench that rules.ELEMENTARY selects."""
    per = t[(t.unit == "ns/elem") & ~t.partial_batch]
    rows = []
    for rank, (facet, variant, layer, op, modes) in enumerate(rules.ELEMENTARY):
        d = per[(per.layer == layer) & (per.operation == op)]
        if (layer, op) in rules.ALSO:
            also_layer, also_op = rules.ALSO[layer, op]
            also = per[(per.layer == also_layer) & (per.operation == also_op)]
            d = pd.concat([d, also])
        if modes is not None:
            d = d[d["mode"].isin(modes)]
        for curve in columns(layer):
            f = d[d.family == curve]
            if f.empty:
                continue
            best = f.loc[f.value_ns.idxmin()]
            rows.append(
                {
                    "operation": f"{facet}, {variant}" if variant else facet,
                    "facet": facet,
                    "variant": variant or "",
                    "layer": layer,
                    "op_rank": rank,
                    "curve": curve,
                    "base": base_family(layer, curve),
                    "bench": best.full_id,
                    "value_ns": best.value_ns,
                    "value_lo_ns": best.value_lo_ns,
                    "value_hi_ns": best.value_hi_ns,
                }
            )
    return pd.DataFrame(
        rows,
        columns=[
            "operation",
            "facet",
            "variant",
            "layer",
            "op_rank",
            "curve",
            "base",
            "bench",
            "value_ns",
            "value_lo_ns",
            "value_hi_ns",
        ],
    )


def columns(layer: str) -> list[str]:
    return list(rules.FIELDS) if layer == "field" else rules.CURVES


def base_family(layer: str, curve: str) -> str:
    return rules.FIELDS[curve] if layer == "field" else rules.family(curve)[1]


# The comparison tables: a row per curve (or field), grouped by field as in
# rules.CURVE_GROUPS, a column per operation and mode, in pipeline order.
GROUP_COLS = [
    Col(("hash to curve", "one at a time"), "hash to curve", "alone"),
    Col(("hash to curve", "batched"), "hash to curve", "batched"),
    Col(("hash to addend", "one at a time"), "hash to addend", "alone"),
    Col(("hash to addend", "batched"), "hash to addend", "batched"),
    Col(("point prepare", "one at a time"), "prepare", "alone"),
    Col(("point prepare", "batched"), "prepare", "batched"),
    Col(("point add", "throughput"), "add", "throughput"),
    Col(("point add", "latency"), "add", "latency"),
    Col(("point subtract", ""), "subtract"),
    Col(("point negate", ""), "negate"),
    Col(("is identity", ""), "is identity"),
    Col(("equals addend", "match"), "equals addend", "match"),
    Col(("equals addend", "mismatch"), "equals addend", "mismatch"),
    Col(("equals addend", ""), "equals addend"),
    Col(("encode", "one at a time"), "encode", "alone"),
    Col(("encode", "batched"), "encode", "batched"),
    Col(("decode", "one at a time"), "decode", "alone"),
    Col(("decode", "batched"), "decode", "batched"),
]
FIELD_ADD = ("field add", "")
FIELD_COLS = [
    Col(("field mul", "throughput"), "mul", "throughput"),
    Col(("field mul", "latency"), "mul", "latency"),
    Col(("field square", "throughput"), "square", "throughput"),
    Col(("field square", "latency"), "square", "latency"),
    Col(FIELD_ADD, "add"),
    Col(("field sqrt", ""), "sqrt"),
    Col(
        ("field sqrt ratio", ""),
        "sqrt ratio",
        note="√(u/v) in one exponentiation, as odd-field hashes take it",
    ),
    Col(("field halftrace", ""), "halftrace", note="qsolve over GF(2^122)"),
    Col(("field invert", "one at a time"), "invert", "alone"),
    Col(("field invert", "batched"), "invert", "batched"),
    Col(("field pack", ""), "pack", note="to the 128-bit wire form"),
    Col(("field unpack", ""), "unpack", note="from the 128-bit wire form"),
]


def grouped(groups, have, compare=lambda c: not c.startswith("xor")):
    """Retain rows with data and nonempty groups as [(heading, [Row])]."""
    out = []
    for heading, rows in groups:
        kept = [
            Row(c, note, compare=compare(c), marks=rules.marks(c))
            for c, note in rows
            if c in have
        ]
        if kept:
            out.append((heading, kept))
    return out


def elementary_grid(e: pd.DataFrame, fields: bool) -> Grid:
    """elementary()'s cells as a Grid: curves (or fields) down, operations
    across, each cell's bench and CI in its tooltip."""
    e = e[(e.layer == "field") == fields]
    cols = FIELD_COLS if fields else GROUP_COLS
    groups = grouped(rules.FIELD_GROUPS if fields else rules.CURVE_GROUPS, set(e.curve))
    cells = {(r.curve, (r.facet, r.variant)): r for r in e.itertuples()}
    for _, rows in groups:
        for row in rows:
            for c in cols:
                r = cells.get((row.label, c.key))
                if r is not None:
                    row.values[c.key] = r.value_ns
                    row.cis[c.key] = (r.value_lo_ns, r.value_hi_ns)
                    row.tips[c.key] = (
                        f"{r.bench}: {fmt_time(r.value_lo_ns)} .. {fmt_time(r.value_hi_ns)}"
                    )
    # a binary field adds by XOR, too fast to time on its own
    if fields:
        for _, rows in groups:
            for row in rows:
                if rules.FIELDS[row.label] == "gf2_127":
                    row.values.setdefault(FIELD_ADD, "xor")
                    row.tips.setdefault(FIELD_ADD, "one XOR per word, not timed")
    used = {k for _, rows in groups for r in rows for k in r.values}
    cols = [c for c in cols if c.key in used]
    return Grid(cols, groups, corner="field" if fields else "curve")


def _tex_int(n: int) -> str:
    # in math a bare comma is punctuation and spaces the digits after it
    return f"{n:,}".replace(",", "{,}")


# the steps benches/group.rs times at several batch sizes
SWEEP = [
    ("h2c", "hash to curve"),
    ("group.prepare", "prepare"),
    ("group.encode", "encode"),
    ("group.decode", "decode"),
]


def batch_size_grid(t: pd.DataFrame) -> Grid | None:
    """Each batched step's time per item at each batch size the benches
    sweep, beside its time alone; None without a sweep."""
    per = t[(t.unit == "ns/elem") & t.group.isin([g for g, _ in SWEEP])]
    sizes = sorted({int(n) for n in per.batch_n.dropna()})
    if len(sizes) < 2:
        return None
    cols, cells = [], {}
    for group, head in SWEEP:
        d = per[per.group == group]
        if d.empty:
            continue
        cols.append(Col((group, "alone"), head, "alone"))
        cols += [Col((group, n), head, f"$n = {n}$") for n in sizes]
        for r in d.itertuples():
            if r.mode == "per-element":
                cells[(str(r.family), (group, "alone"))] = r
            elif r.mode == "batch" and pd.notna(r.batch_n):
                cells[(str(r.family), (group, int(r.batch_n)))] = r
    groups = grouped(rules.CURVE_GROUPS, {f for f, _ in cells})
    for _, rows in groups:
        for row in rows:
            for c in cols:
                r = cells.get((row.label, c.key))
                if r is not None:
                    row.values[c.key] = r.value_ns
                    row.cis[c.key] = (r.value_lo_ns, r.value_hi_ns)
                    row.tips[c.key] = r.full_id
    return Grid(cols, groups, corner="curve")


# the dominant operations: (layer, operation, head, sub)
COVERAGE = [
    ("hash to curve", "hash to curve", "hash to curve", ""),
    ("group ops", "prepare", "prepare", ""),
    ("group ops", "add", "add", ""),
    ("group ops", "subtract", "subtract", ""),
    ("group ops", "negate", "negate", ""),
    ("group ops", "is identity", "is identity", ""),
    ("group ops", "equals addend", "equals addend", ""),
    ("group ops", "encode", "encode", ""),
    ("group ops", "decode", "decode", ""),
    ("RIBLT workload", "encode (hash + cells)", "RIBLT", "encode"),
    ("RIBLT workload", "cell updates", "RIBLT", "cells"),
    ("RIBLT workload", "peel, per difference", "RIBLT", "peel"),
]


def coverage(t: pd.DataFrame) -> pd.DataFrame:
    """The number of benchmarks of each hot-path operation on each family
    rules.CURVE_GROUPS lists, measured or not."""
    families = [c for _, rows in rules.CURVE_GROUPS for c, _ in rows]
    n = t.groupby([t.family.astype(str), "layer", "operation"], observed=True).size()
    return pd.DataFrame(
        [[n.get((f, layer, op), 0) for layer, op, _, _ in COVERAGE] for f in families],
        index=pd.Index(families, name="family"),
        columns=[f"{h}, {s}" if s else h for _, _, h, s in COVERAGE],
    )


def coverage_grid(c: pd.DataFrame) -> Grid:
    """coverage() as a Grid, a missing operation an empty cell."""
    cols = [Col(name, *name.split(", ", 1), kind=PLAIN) for name in c.columns]
    groups = grouped(rules.CURVE_GROUPS, set(c.index))
    for _, rows in groups:
        for row in rows:
            row.values = {k: int(v) for k, v in c.loc[row.label].items() if v}
    return Grid(cols, groups, corner="curve")


def display_table(d: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "operation": d.operation,
            "benchmark": d.label,
            "family": d.family.astype(str),
            "mode": d["mode"],
            "time": [fmt_time(v) for v in d.value_ns],
            "95% CI": [
                f"{fmt_time(a)} .. {fmt_time(b)}"
                for a, b in zip(d.value_lo_ns, d.value_hi_ns)
            ],
            "unit": d.unit.str.replace("ns/", "per "),
        }
    )


def beside(root: Path, name: str) -> Path | None:
    """A file bench-run writes, in the run directory or beside criterion's."""
    for d in (root, root.parent):
        if (d / name).is_file():
            return d / name
    return None


def machine(root: Path) -> dict | None:
    """bench-run's meta.json, named for a reader by label() when bench-run
    identified the run by its id alone."""
    found = beside(root, "meta.json")
    if not found:
        return None
    m = json.loads(found.read_text())
    if "name" not in m and "id" in m:
        m["name"] = label(m)
    return m


def label(m: dict) -> str:
    """A run's name in reports: machine, build target, profile and commit,
    since its UUIDv7 says nothing to a reader."""
    parts = [m.get("cpu"), m.get("target"), m.get("profile"), m.get("commit", "")[:8]]
    return " · ".join(p for p in parts if p)


# bench-run's samples of the load averages, taken every few seconds
LOAD = "load.tsv"


def machine_load(root: Path) -> pd.DataFrame | None:
    """For each phase of a run in bench-run's load samples, in order: its
    length, the mean of the 1-minute load average, and the threads runnable,
    on average and at most. Linux reports the runnable count at each sample;
    otherwise it is recovered from successive 1-minute averages, which the
    kernel decays by e^(-t/60) over t seconds: n = (L_k - L_(k-1) e^(-t/60))
    / (1 - e^(-t/60)), exact for samples on the kernel's 5-second updates,
    an estimate of the mean over the interval otherwise."""
    found = beside(root, LOAD)
    if not found:
        return None
    s = pd.read_csv(found, sep="\t")
    dt = s.time.diff()
    decay = (-dt / 60).map(math.exp)
    recovered = (s.load1 - s.load1.shift() * decay) / (1 - decay)
    s["runnable"] = s.runnable.where(s.runnable.notna(), recovered)
    s = s[dt > 0]
    if s.empty:
        return None
    g = s.groupby("phase", sort=False)
    return pd.DataFrame(
        {
            "seconds": g.size() * dt[dt > 0].median(),
            "load1": g.load1.mean(),
            "runnable": g.runnable.mean(),
            "runnable_max": g.runnable.max(),
        }
    )


def load_table(load: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "phase": load.index,
            "seconds": [f"{v:.0f}" for v in load.seconds],
            "1-minute average, mean": [f"{v:.2f}" for v in load.load1],
            "runnable, mean": [f"{v:.2f}" for v in load.runnable],
            "runnable, most": [f"{v:.1f}" for v in load.runnable_max],
        }
    )


# meta.json's keys that name this machine rather than describe it
UNPUBLISHED = {"host"}


def published(meta: dict) -> str:
    # machine() derives a name from an id; the export keeps what bench-run wrote
    derived = {"name"} if "id" in meta else set()
    return (
        json.dumps(
            {k: v for k, v in meta.items() if k not in UNPUBLISHED | derived},
            indent=2,
        )
        + "\n"
    )


def run_name(root: Path, meta: dict | None) -> str:
    """The run's name, not where it sits on this machine: bench-run's, else
    its directory's."""
    if meta and meta.get("name"):
        return meta["name"]
    return (root.parent if root.name == "criterion" else root).name


# the target features that pick a field backend, of all rustc lists
BACKEND_FEATURES = {"pclmulqdq", "sse4.1", "avx2", "aes", "neon"}


def machine_table(m: dict) -> pd.DataFrame:
    gib = m.get("memory", 0) / 2**30
    rows = [
        ("CPU", f"{m.get('cpu', '?')}, {m.get('cores', '?')} cores, {gib:.0f} GiB"),
        ("OS", f"{m.get('os', '?')} ({m.get('kernel', '?')}, {m.get('arch', '?')})"),
        ("frequency governor", m.get("governor")),
        ("rustc", m.get("rustc")),
        ("RUSTFLAGS", m.get("rustflags") or "(none)"),
        (
            "field backend features",
            ", ".join(f for f in m.get("target_features", []) if f in BACKEND_FEATURES)
            or "none: portable multiplies",
        ),
        (
            "commit",
            m.get("commit")
            and m["commit"] + (" (uncommitted edits)" if m.get("dirty") else ""),
        ),
        ("suites", " ".join(m.get("suites", []))),
        ("profile", m.get("profile")),
        ("criterion arguments", m.get("criterion_args")),
        ("run", f"{m.get('started', '?')} .. {m.get('finished', 'unfinished')}"),
    ]
    return pd.DataFrame(
        [(k, v) for k, v in rows if v], columns=["", m.get("name", "machine")]
    )


def provenance(m: dict | None) -> str:
    """One line for under every figure and the page: which run, machine,
    build and sources; empty without meta.json."""
    if not m:
        return ""
    commit = m.get("commit")
    parts = [
        # a label() repeats the parts below
        None if "id" in m else m.get("name"),
        m.get("cpu"),
        " ".join(m.get("rustc", "").split()[:2]),
        f"RUSTFLAGS={m.get('rustflags') or '(none)'}",
        commit and commit[:12] + (" + edits" if m.get("dirty") else ""),
        m.get("profile") and f"{m['profile']} profile",
        m.get("started", "")[:10],
    ]
    return " · ".join(p for p in parts if p)


def rho_bits(r: int, automorphisms: int) -> float:
    """log2 of rho's expected group operations in a group of prime order r,
    sqrt(pi r / (2 a)), with a the order of the automorphism group the walk
    quotients by: nominal, the generic bound and no structural attack."""
    return (math.log2(math.pi) + math.log2(r) - math.log2(2 * automorphisms)) / 2


def fixtures(t, meta) -> dict[str, dict]:
    """Each observed group-suite family's certified fixture, as bench-run
    recorded the executable's declaration, or None where it recorded none.
    Runs before 2026-10-03 recorded only the status: their fixtures have no
    r, so no rho."""
    records = (meta or {}).get("group_fixtures", {})
    if not isinstance(records, dict):
        records = {}
    group = t[t.group.str.startswith("group.") | (t.group == "h2c")]
    families = sorted({f.split("/")[0] for f in group.function if isinstance(f, str)})
    out = {}
    for f in families:
        rec = records.get(f)
        out[f] = {"status": rec} if isinstance(rec, str) else rec
    return out


def blocks(
    t,
    drawn,
    skipped,
    source,
    elem=None,
    elem_fig=None,
    meta=None,
    load=None,
    runs=None,
    selection_fig=None,
):
    b = [("h1", "Benchmark report")]
    intro = (
        f"{len(t)} benchmarks from run {source}. [[methodology.md#reading-a-run-report|Reading a run report]] explains the values, "
        "colours and tags."
    )
    b.append(("p", intro))
    b.append(("h2", "Machine"))
    if meta:
        b.append(("table", machine_table(meta)))
    else:
        # Criterion output alone does not identify the execution environment.
        unknown = (
            "Unknown: no meta.json accompanies these results. Machine, compiler "
            "flags, field backends and source revision are therefore unrecorded. "
            "The bench-run command records this metadata."
        )
        b.append(("p", unknown))
    if load is not None:
        b.append(
            (
                "p",
                (
                    "Load averages by phase while the run timed "
                    "([[methodology.md#measurement-methodology|Measurement methodology]])."
                ),
            )
        )
        b.append(("table", load_table(load)))
    if runs is not None and len(runs):
        b.extend(selection_blocks(runs, selection_fig, t))
    if elem is not None and len(elem):
        b.extend(elementary_blocks(elem, elem_fig))
    b.append(("h2", "Coverage"))
    have = (
        "Benchmarks per hot-path operation and family; an empty cell was not measured."
    )
    b += [("p", have), ("grid", coverage_grid(coverage(t)))]
    if skipped:
        b.append(("h2", "Skipped"))
        b.append(("p", "Unreadable results: " + "; ".join(skipped)))
    b.append(("h2", "All benchmarks"))
    every = "Every observation: a figure and a table per layer."
    b.append(("p", every))
    for layer, fig in drawn.items():
        d = t[t.layer == layer]
        b.append(("h3", layer))
        if layer in rules.LAYER_NOTES:
            b.append(("p", rules.LAYER_NOTES[layer]))
        if fig:
            b.append(("details", "Figure", [("img", fig)]))
        b.append(("details", f"{len(d)} benchmarks", [("table", display_table(d))]))
    return b


def selection_blocks(runs, fig=None, t=None):
    """The curve-selection section: what finding, proving and verifying
    each family's curve took per seed, the figure for the method that finds
    fastest, and every method's means; with the run's benchmarks `t`, the
    embedding-degree bound's share of accepting a curve."""
    chosen = fastest_method(runs)
    what = (
        "Curve selection as this run timed it, once per seed "
        f"({runs.seed.nunique()} seeds): find searches a seed's candidates for "
        "an acceptable curve; prove searches them again and certifies each "
        "rejection and the acceptance; verify checks that certificate. The "
        "figure shows each family's method that finds fastest, a bar at the "
        "seeds' mean in a box from the fastest seed to the slowest, and how "
        "the prover's time divides among the search's phases; the table lists "
        "every method ([[problem.md#parameter-selection-and-verification|Parameter selection and verification]])."
    )
    if t is not None:
        what += embedding_share(t)
    g = runs.groupby(["family", "method"])
    mean = g.mean(numeric_only=True)
    count = mean.prove_count_s
    table = pd.DataFrame(
        {
            "seeds": g.seed.count(),
            "candidates": mean.candidates.round(1),
            "point counts (prove)": mean.prove_counts.round(1),
            "find": [fmt_time(v * 1e9) for v in mean.find_s],
            "prove": [fmt_time(v * 1e9) for v in mean.prove_s],
            "verify": [fmt_time(v * 1e9) for v in mean.verify_s],
            "counting, of prove": [f"{100 * v:.0f}%" for v in count / mean.prove_s],
        }
    ).reset_index()
    shown = set(zip(chosen.family, chosen.method))
    table.insert(
        2,
        "shown",
        ["yes" if k in shown else "" for k in zip(table.family, table.method)],
    )
    b = [("h2", "Curve selection"), ("p", what)]
    if fig:
        b.append(("img", fig))
    b.append(("table", table))
    return b


def embedding_share(t: pd.DataFrame) -> str:
    """A sentence on the embedding-degree bound's time against the whole
    acceptance check, over the families the run timed both for, or ""."""

    def times(group):
        d = t[t.group == f"curvegen/{group}"]
        return pd.Series(d.value_ns.values, index=d.family.astype(str))

    e, a = times("embedding"), times("verify_accept")
    both = e.index.intersection(a.index)
    if both.empty:
        return ""
    e, a = e[both], a[both]
    share = e / a
    return (
        f" Accepting a curve ({fmt_time(a.min())} to {fmt_time(a.max())}, "
        "curvegen/verify_accept) includes the embedding-degree bound, a "
        "baby-step giant-step search for $q^k = 1 \\bmod r$ with "
        f"$k \\le 2^{{20}}$, which alone takes {fmt_time(e.min())} to "
        f"{fmt_time(e.max())} (curvegen/embedding), {100 * share.min():.0f}% to "
        f"{100 * share.max():.0f}% of it; rejected candidates do not reach it."
    )


def elementary_blocks(elem, elem_fig=None):
    what = (
        "The fastest measurement of each operation and mode per family, the "
        "hashes over every construction; rows grouped by field."
    )
    b = [("h2", "Group operations"), ("p", what)]
    if ((elem.layer == "field") == False).any():
        b.append(("grid", elementary_grid(elem, fields=False)))
    if "group" in (elem_fig or {}):
        b.append(("img", elem_fig["group"]))
    if (elem.layer == "field").any():
        b.append(("h2", "Field operations"))
        made = "Per field; batched inversion per element of a batch of 1024."
        b.append(("p", made))
        b.append(("grid", elementary_grid(elem, fields=True)))
        if "field" in (elem_fig or {}):
            b.append(("img", elem_fig["field"]))
    cells = pd.DataFrame(
        {
            "operation": elem.operation,
            "curve": elem.curve,
            "bench": elem.bench,
            "time": [fmt_time(v) for v in elem.value_ns],
            "95% CI": [
                f"{fmt_time(a)} .. {fmt_time(c)}"
                for a, c in zip(elem.value_lo_ns, elem.value_hi_ns)
            ],
        }
    )
    b.append(("details", "Which bench each cell is", [("table", cells)]))
    return b


def _md_cell(v) -> str:
    # renderers take hash_from_bytes<Sha512> for an HTML tag and drop it
    return html.escape(str(v), quote=False).replace("|", "\\|")


def _md_alt(s: str) -> str:
    return s.replace("[", "(").replace("]", ")")


DOC_LINK = re.compile(r"\[\[([^|\]]+)\|([^\]]+)\]\]")


def docs(text: str, book: bool) -> str:
    """Resolve [[chapter.md#anchor|title]]: a link in the book, whose run
    pages sit two levels below the chapters, else the title and the file."""
    if book:
        return DOC_LINK.sub(lambda m: f"[{m[2]}](../../{m[1]})", text)
    return DOC_LINK.sub(lambda m: f"{m[2]} (docs/{m[1].split('#')[0]})", text)


def to_markdown(bs, html_grids=False, book=False) -> str:
    """Markdown of the blocks; with html_grids, the comparison grids as the
    HTML report has them, shaded and with tooltips, for a renderer that
    passes HTML through and is given tables.STYLE; with book, laid out for
    BOOK_STYLE: each figure as its light and dark twins, which it shows by
    the theme, and each grid in a wrapper that may widen it."""
    out = []
    for kind, v, *rest in bs:
        if kind == "details":
            # GitHub renders markdown inside <details> given the blank lines
            inner = to_markdown(rest[0], html_grids, book).rstrip("\n")
            out.append(
                f"<details><summary>{html.escape(v)}</summary>\n\n{inner}\n\n</details>"
            )
        elif kind == "h1":
            out.append(f"# {v}")
        elif kind == "h2":
            out.append(f"## {v}")
        elif kind == "h3":
            out.append(f"### {v}")
        elif kind == "p":
            out.append(docs(html.escape(v, quote=False), book))
        elif kind == "footer":
            out.append("---\n\n" + html.escape(v, quote=False))
        elif kind == "img" and book:
            src, alt = html.escape(figures.svg(v)), html.escape(v.alt)
            out.append(
                f'<img class="fig-light" src="{src}" alt="{alt}">'
                f'<img class="fig-dark" src="{html.escape(figures.dark(figures.svg(v)))}" '
                f'alt="{alt}">'
            )
        elif kind == "img":
            out.append(f"![{_md_alt(v.alt)}]({figures.svg(v)})")
        elif kind == "grid" and html_grids and book:
            out.append(f'<div class="grid-wrap">{tables.to_html(v, fmt_time)}</div>')
        elif kind == "grid" and html_grids:
            out.append(tables.to_html(v, fmt_time))
        elif kind == "grid":
            out.append(tables.to_markdown(v, fmt_time, _md_cell))
        elif kind == "table":
            rows = [list(v.columns)] + v.values.tolist()
            lines = ["| " + " | ".join(_md_cell(c) for c in r) + " |" for r in rows]
            lines.insert(1, "|" + "---|" * len(v.columns))
            out.append("\n".join(lines))
    return "\n\n".join(out) + "\n"


# the text carries $...$ mathematics, as GitHub and the book render it
KATEX = (
    '<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16/dist/katex.min.css">'
    '<script defer src="https://cdn.jsdelivr.net/npm/katex@0.16/dist/katex.min.js"></script>'
    '<script defer src="https://cdn.jsdelivr.net/npm/katex@0.16/dist/contrib/auto-render.min.js"'
    """ onload="renderMathInElement(document.body, {delimiters: [{left: '$', right: '$', display: false}]})"></script>"""
)


def to_html(bs) -> str:
    return (
        '<!doctype html>\n<html><head><meta charset="utf-8"><title>Benchmark report</title>'
        f"<style>{STYLE}</style>{KATEX}</head><body>\n"
        + _html_body(bs)
        + "\n</body></html>\n"
    )


def _html_body(bs) -> str:
    body = []
    for kind, v, *rest in bs:
        if kind == "details":
            body.append(
                f"<details><summary>{html.escape(v)}</summary>\n"
                f"{_html_body(rest[0])}\n</details>"
            )
        elif kind in ("h1", "h2", "h3"):
            body.append(f"<{kind}>{html.escape(v)}</{kind}>")
        elif kind == "p":
            body.append(f"<p>{docs(html.escape(v), False)}</p>")
        elif kind == "footer":
            body.append(f"<footer>{html.escape(v)}</footer>")
        elif kind == "img":
            body.append(
                f'<p><img src="{html.escape(figures.svg(v))}" '
                f'alt="{html.escape(v.alt)}"></p>'
            )
        elif kind == "grid":
            body.append(tables.to_html(v, fmt_time))
        elif kind == "table":
            body.append(v.to_html(index=False, border=0))
    return "\n".join(body)


STYLE = (
    f"body{{font-family:system-ui,sans-serif;color:{TEXT};background:{SURFACE};"
    "max-width:72rem;margin:2rem auto;padding:0 1rem}"
    "img{max-width:100%}table{border-collapse:collapse;font-size:.8rem}"
    f"td,th{{padding:.15rem .6rem;text-align:left;border-bottom:1px solid {GRID}}}"
    "td:not(:first-child){font-variant-numeric:tabular-nums}"
    "summary{cursor:pointer;color:" + TEXT_2 + ";margin:.5rem 0}"
    "p{max-width:52rem;line-height:1.45}"
    f"footer{{margin-top:2rem;padding-top:.5rem;border-top:1px solid {GRID};"
    f"color:{TEXT_2};font-size:.75rem}}" + tables.STYLE
)


def analyse(
    root: Path, out: Path, stat="typical", formats=("png", "svg")
) -> tuple[pd.DataFrame, list]:
    """Write a run's figures and tables into out, and return its tidy table
    and the blocks of its report, closed by the footer."""
    df, skipped = load(root, stat)
    if df.empty:
        raise SystemExit(f"no criterion results under {root} (*/new/benchmark.json)")
    t = tidy(df)
    odd = unclassified(t)
    if len(odd):
        ids = "\n".join(
            f"  {r.full_id}  (layer {r.layer}, family {r.family}, operation {r.operation})"
            for r in odd.itertuples()
        )
        raise SystemExit(
            f"{len(odd)} benchmark ids match no rule in report/rules.py:\n{ids}"
        )
    out.mkdir(parents=True, exist_ok=True)
    meta = machine(root)
    footer = provenance(meta)
    figures.theme()
    drawn = {
        layer: plot_layer(t, layer, out, formats, footer)
        for layer in t.layer.cat.categories
        if (t.layer == layer).any()
    }
    elem = elementary(t)
    elem_fig = plot_elementary(elem, out, formats, footer) if len(elem) else {}
    runs = selection_runs(root)
    # runs whose curvegen.csv predates proving and its phases have no section
    runs = runs if "prove_count_s" in runs else None
    selection_fig = (
        plot_selection(fastest_method(runs), out, formats, footer)
        if runs is not None and len(runs)
        else None
    )
    t.to_csv(out / "tidy.csv", index=False)
    elem.to_csv(out / "elementary.csv", index=False)
    coverage(t).to_csv(out / "coverage.csv")
    if meta:
        (out / "meta.json").write_text(published(meta))
    bs = blocks(
        t,
        drawn,
        skipped,
        run_name(root, meta),
        elem,
        elem_fig,
        meta,
        machine_load(root),
        runs=runs,
        selection_fig=selection_fig,
    )
    bs.append(("footer", footer or run_name(root, meta)))
    for s in skipped:
        print(f"skipped {s}", file=sys.stderr)
    return t, bs


def report(
    root: Path, out: Path, stat="typical", formats=("png", "svg")
) -> pd.DataFrame:
    t, bs = analyse(root, out, stat=stat, formats=formats)
    (out / "report.md").write_text(to_markdown(bs))
    (out / "report.html").write_text(to_html(bs))
    print(
        f"{len(t)} benchmarks -> {out / 'report.md'}",
        file=sys.stderr,
    )
    return t


def ran(m: dict | None) -> str:
    """A sentence stating where, when and from what source a run timed."""
    if not m:
        return (
            "No meta.json accompanies this run: its machine, build and source "
            "revision are unrecorded."
        )
    run = f"Run {m['name']}" if m.get("name") else "The run"
    commit = m.get("commit", "an unrecorded commit")
    edits = " with uncommitted edits" if m.get("dirty") else ""
    return (
        f"{run} started {m.get('started', 'at an unrecorded time')} on "
        f"{m.get('cpu', 'an unrecorded CPU')} ({m.get('os', '?')}, "
        f"{m.get('arch', '?')}) and timed commit {commit}{edits}, built by "
        f"{m.get('rustc', 'an unrecorded rustc')} with RUSTFLAGS "
        f"{m.get('rustflags') or '(none)'}, in the {m.get('profile', '?')} profile."
    )


def chapter(root: Path, out: Path, title: str) -> str:
    """A run's report as a chapter of the book, under title and the
    sentence of ran(), its figures as SVG in out, each with a twin for the
    book's dark themes, and its grids shaded and with tooltips as the report
    states them to be."""
    _, bs = analyse(root, out, formats=("svg",))
    for name in _figures(bs):
        dark_twin(out / name)
    page = [("h1", title), ("p", ran(machine(root))), *bs[1:]]
    return f"<style>{tables.STYLE}{BOOK_STYLE}</style>\n\n" + to_markdown(
        page, html_grids=True, book=True
    )


# mdbook's dark themes name themselves on <html>, as in fungi-protocol/docs;
# print keeps the light rendering whatever the reader's theme
_DARK_THEMES = "html:is(.coal,.navy,.ayu)"
BOOK_STYLE = (
    # a grid wider than the text column is centred on it and may use the
    # width beside it, short of the sidebar and the page-turning arrows,
    # and scrolls beyond that
    ".grid-wrap{overflow-x:auto;width:max-content;"
    "max-width:max(100%,calc(100vw - var(--sidebar-width) - 200px));"
    "position:relative;left:50%;transform:translateX(-50%)}"
    "html.sidebar-hidden .grid-wrap{max-width:max(100%,calc(100vw - 200px))}"
    ".fig-dark{display:none}@media screen{"
    f"{_DARK_THEMES} .fig-dark{{display:inline}}"
    f"{_DARK_THEMES} .fig-light{{display:none}}"
    f"{_DARK_THEMES} table.grid :is(th,td){{border-color:var(--table-border-color)}}"
    f"{_DARK_THEMES} table.grid :is(th.sub,th.row .note,td.baseline)"
    "{color:color-mix(in srgb,currentColor 70%,transparent)}"
    f"{_DARK_THEMES} table.grid tr.group th"
    "{background:var(--table-header-bg);color:inherit}"
    f"{_DARK_THEMES} table.grid th.row .mark"
    "{background:color-mix(in srgb,#f1e4c8 20%,transparent);color:#e6c98f}}"
)


def _figures(bs):
    """The SVG file of each figure in the blocks."""
    for kind, v, *rest in bs:
        if kind == "details":
            yield from _figures(rest[0])
        elif kind == "img":
            yield figures.svg(v)


def export(run: Path, dest: Path) -> pd.DataFrame:
    """Write the compact form of a run into dest, from which the report
    renders without criterion's tree: each benchmark's RAW fields in
    benchmarks.csv, the published fields of meta.json, curvegen.csv
    when the run timed curve selection, and load.tsv when bench-run sampled
    the load. A run with unreadable results or
    without meta.json is refused, as its report would be incomplete or
    unattributed; so is a run of a source with uncommitted edits, which no
    published revision reproduces; and so is a dest that exists, which
    could keep a file of another run."""
    criterion = run / "criterion" if (run / "criterion").is_dir() else run
    df, skipped = load(criterion)
    if df.empty:
        raise SystemExit(f"no criterion results under {criterion}")
    if skipped:
        raise SystemExit("unreadable results:\n  " + "\n  ".join(skipped))
    meta = machine(criterion)
    if not meta:
        raise SystemExit(f"no meta.json in {run}: the run is unattributed")
    if not meta.get("commit") or meta.get("dirty") is not False:
        raise SystemExit(
            f"{run} was not run from a committed source: no revision reproduces it"
        )
    if dest.exists():
        raise SystemExit(f"{dest} exists: a run is exported once, into a new directory")
    dest.mkdir(parents=True)
    raw = df[RAW].astype({"elements": "Int64"})
    raw.to_csv(dest / COMPACT, index=False)
    (dest / "meta.json").write_text(published(meta))
    for name in ("curvegen.csv", LOAD):
        found = beside(criterion, name)
        if found:
            shutil.copyfile(found, dest / name)
    print(f"{len(df)} benchmarks -> {dest / COMPACT}", file=sys.stderr)
    return raw


INDEX_COLS = ["run", "started", "CPU", "rustc", "RUSTFLAGS", "commit", "profile"]


def index(runs: Path) -> pd.DataFrame:
    """List each report under runs, newest first, from the meta.json beside
    it, in runs/index.html and runs/index.md; the listing."""
    rows = []
    for page in sorted(runs.glob("**/report.html")):
        m = machine(page.parent) or {}
        commit = m.get("commit")
        rows.append(
            {
                "run": m.get("name") or page.parent.name,
                "href": page.relative_to(runs).as_posix(),
                "started": m.get("started", ""),
                "CPU": m.get("cpu", ""),
                "rustc": " ".join(m.get("rustc", "").split()[:2]),
                "RUSTFLAGS": m.get("rustflags", ""),
                "commit": commit
                and commit[:12] + (" + edits" if m.get("dirty") else ""),
                "profile": m.get("profile", ""),
            }
        )
    if not rows:
        raise SystemExit(f"no report.html under {runs}")
    listing = pd.DataFrame(rows).sort_values("started", ascending=False)
    esc = html.escape
    head = "".join(f"<th>{esc(c)}</th>" for c in INDEX_COLS)
    body = "".join(
        f'<tr><td><a href="{esc(r.href)}">{esc(r.run)}</a></td>'
        + "".join(f"<td>{esc(str(r[c] or ''))}</td>" for c in INDEX_COLS[1:])
        + "</tr>"
        for _, r in listing.iterrows()
    )
    (runs / "index.html").write_text(
        '<!doctype html>\n<html><head><meta charset="utf-8"><title>Benchmark runs'
        f"</title><style>{STYLE}</style></head><body>\n<h1>Benchmark runs</h1>\n"
        f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>\n"
        "</body></html>\n"
    )
    lines = ["# Benchmark runs", "", "| " + " | ".join(INDEX_COLS) + " |"]
    lines.append("|" + "---|" * len(INDEX_COLS))
    for _, r in listing.iterrows():
        md = r.href.removesuffix(".html") + ".md"
        cells = [f"[{_md_cell(r.run)}]({md})"]
        cells += [_md_cell(r[c] or "") for c in INDEX_COLS[1:]]
        lines.append("| " + " | ".join(cells) + " |")
    (runs / "index.md").write_text("\n".join(lines) + "\n")
    return listing


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "criterion",
        type=Path,
        nargs="?",
        help="criterion output, e.g. target/criterion",
    )
    p.add_argument("out", type=Path, nargs="?", help="output directory")
    p.add_argument(
        "--index", type=Path, metavar="DIR", help="list the reports under DIR instead"
    )
    p.add_argument(
        "--export",
        action="store_true",
        help="write the compact form of the run CRITERION into OUT instead",
    )
    p.add_argument(
        "--stat", default="typical", choices=["typical", "slope", "mean", "median"]
    )
    p.add_argument("--formats", default="png,svg")
    p.add_argument(
        "--check-ids",
        action="store_true",
        help="read benchmark ids from stdin, one per line (criterion's --list), "
        "and fail if there are none, or naming those no rule in rules.py classifies",
    )
    a = p.parse_args(argv)
    if a.check_ids:
        ids = [i for i in (ln.strip() for ln in sys.stdin) if i]
        if not ids:
            # a filter that selects nothing, or a --list format sed no
            # longer matches: either way nothing was checked
            raise SystemExit("no benchmark ids on stdin: the filter selects none")
        odd = check_ids(ids)
        if odd:
            raise SystemExit(
                f"{len(odd)} benchmark ids match no rule in report/rules.py:\n  "
                + "\n  ".join(odd)
            )
        return
    if a.index:
        index(a.index)
        return
    if a.out is None:
        p.error("give the criterion output and an output directory, or --index")
    if a.export:
        export(a.criterion, a.out)
        return
    report(a.criterion, a.out, stat=a.stat, formats=tuple(a.formats.split(",")))


if __name__ == "__main__":
    main()
