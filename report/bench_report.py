"""Criterion results -> tidy table, facet-grid plots and a markdown/HTML report.

  nix run .#bench-report -- target/criterion out/
  nix run .#bench-report -- --export bench-runs/ID results/ID
  nix run .#bench-report -- --index runs/
  BENCH --list | sed 's/: benchmark$//' | nix run .#bench-report -- --check-ids
  riblt --list | sed 's/: benchmark$//' | nix run .#bench-report -- --riblt-plan CRITERION
  group --list | sed 's/: benchmark$//' | nix run .#bench-report -- --group-plan CRITERION

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
from typing import NamedTuple

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
    plot_insert,
    plot_layer,
    plot_selection,
)
from tables import PLAIN, PLAIN_FACTOR, PLAIN_TIME, Col, Grid, Row


class Load(NamedTuple):
    """How many cells an item is inserted into: k."""

    key: str  # the cost tables' column
    label: str
    k: float


def mapping_degree(m: int) -> float:
    """The cells an item maps to among the first m, on average: the RIBLT
    mapping (src/riblt.rs) includes cell i with probability 1 / (1 + i/2),
    which sums to 2(H_{m+1} − 1)."""
    return 2 * sum(1 / j for j in range(2, m + 2))


def riblt_loads(ms) -> list[Load]:
    """A RIBLT's inserts at each of the coded-symbol counts ms."""
    return [
        Load(f"m={m}", f"$m = {m}$: $k = {mapping_degree(m):.2f}$", mapping_degree(m))
        for m in ms
    ]


def fixed_loads(ks) -> list[Load]:
    return [Load(f"k={k:g}", f"$k = {k:g}$", k) for k in ks]


# the table sizes benches/riblt.rs encodes into (its MS)
DEFAULT_MS = (5, 20, 150, 1350, 12150)
DEFAULT_LOADS = riblt_loads(DEFAULT_MS)


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


def h2c_token(t: pd.DataFrame, bench: str | None) -> str | None:
    """The h2c= parameter of benches/riblt.rs naming the construction a
    measured hash benchmark times, or None if RIBLT has no such parameter."""
    if bench is None:
        return None
    row = t[t.full_id == bench.split(", in batches of ")[0]]
    if row.empty:
        return None
    row = row.iloc[0]
    for token, (layer, op, pattern) in rules.H2C.items():
        if (
            row.layer == layer
            and row.operation == op
            and re.search(pattern, row.full_id)
        ):
            return token
    return None


def riblt_plan(
    t: pd.DataFrame, ids, load: str = "m=1350"
) -> dict[str, tuple[str, str]]:
    """How benches/riblt.rs is to time each family the listed ids offer
    under several hash constructions: (h2c, scope). The construction is the
    one whose modeled insertion (pipeline_cost) at `load` costs least among
    those RIBLT offers; the scope follows that cost's ratio to the least of
    any family, by rules.RIBLT_SCOPES. A family without a model is absent,
    and benches/riblt.rs spot-checks it under every construction."""
    offered: dict[str, set[str]] = {}
    for i in ids:
        group, _, rest = i.strip().partition("/")
        fam, _, param = rest.partition("/")
        m = re.search(r"(?:^|,)h2c=([^,/]+)", param)
        if group.startswith("riblt.") and m:
            offered.setdefault(fam, set()).add(m[1])
    model = pipeline_cost(t, [ld for ld in DEFAULT_LOADS if ld.key == load])
    best = {}
    for fam, tokens in sorted(offered.items()):
        family = rules.family(fam, "riblt.encode")[0]
        recipes = model[model.family.astype(str) == family] if len(model) else model
        costs = {}
        for _, r in recipes.iterrows():
            token = h2c_token(t, r["hash bench"])
            if token in tokens and r[load] < costs.get(token, math.inf):
                costs[token] = r[load]
        if costs:
            token = min(costs, key=costs.get)
            best[fam] = token, costs[token]
    if not best:
        return {}
    least = min(cost for _, cost in best.values())
    return {
        fam: (
            token,
            next(
                (scope for scope, ratio in rules.RIBLT_SCOPES if cost <= ratio * least),
                "spot",
            ),
        )
        for fam, (token, cost) in best.items()
    }


def group_plan(t: pd.DataFrame, ids, load: str = "m=1350") -> dict[str, str]:
    """Which families of the listed benches/group.rs ids its second pass
    (GROUP_PASS=wide) is to time, with their role: those whose modeled
    insertion (pipeline_cost) at `load` is within rules.RIBLT_SCOPES's
    buffer ratio of the least ("lead"), and, if none of them is over a
    field of odd characteristic, the cheapest family that is
    ("contrast")."""
    offered = {}
    for i in ids:
        group, _, rest = i.strip().partition("/")
        if group == "h2c" or group.startswith("group."):
            spelled = rest.partition("/")[0]
            offered[rules.family(spelled, group)[0]] = spelled
    model = pipeline_cost(t, [ld for ld in DEFAULT_LOADS if ld.key == load])
    if not len(model):
        return {}
    model = model[model.family.astype(str).isin(offered)]
    cost = model.groupby(model.family.astype(str))[load].min().dropna().sort_values()
    if cost.empty:
        return {}
    ratio = dict(rules.RIBLT_SCOPES)["buffer"]
    lead = list(cost.index[cost <= ratio * cost.iloc[0]])
    plan = {offered[f]: "lead" for f in lead}
    odd = [f for f in cost.index if rules.FAMILIES[f][0] != "gf2_127"]
    if odd and not set(odd) & set(lead):
        plan[offered[odd[0]]] = "contrast"
    return plan


def _one(rows: pd.DataFrame, what: str):
    """The one benchmark a cost term reads: with two, which one it reads
    would depend on their order."""
    if len(rows) != 1:
        raise ValueError(f"{what}: {len(rows)} benchmarks, {list(rows.full_id)}")
    return rows.iloc[0]


def sum_ci(terms) -> tuple[float, float, float]:
    """A weighted sum of independent estimates, [(benchmark or None,
    weight)], and its interval: each side's half-widths in quadrature, a
    weighted one scaled, being one estimate counted several times."""
    terms = [(r, w) for r, w in terms if r is not None]
    if not terms:
        return (float("nan"),) * 3
    v = sum(w * r.value_ns for r, w in terms)
    lo = math.hypot(*(w * (r.value_ns - r.value_lo_ns) for r, w in terms))
    hi = math.hypot(*(w * (r.value_hi_ns - r.value_ns) for r, w in terms))
    return v, v - lo, v + hi


def described(r) -> str | None:
    """A term's benchmark, and its batch size if it batches."""
    if r is None:
        return None
    if r["mode"] != "batch":
        return r.full_id
    n = r.batch_n if pd.notna(r.batch_n) else r.elements
    return f"{r.full_id}, in batches of {n:g}"


def comparison_recipe(row):
    """Known measured output contracts; an unrecognized map supplies no recipe.

    binary/map.rs returns Point, already binary::Accumulate::Addend.
    Elligator2 returns the family's affine point, and SSWU returns Weier
    affine: both feed the corresponding group.prepare operation.
    """
    if row.layer != rules.EXPERIMENTAL:
        return None
    name = str(row.function).split("/", 1)[-1].removesuffix(", batched")
    field = re.fullmatch(r"gf2_\d+(?:-gls)?(-u|-lambda|-w)?", str(row.family))
    if field and field[1] is None and name in ("pornin map x1", "pornin map x2"):
        return f"Pornin {name[-2:]}: hash to addend", False
    # the w codec's addends are λ-affine, as the λ family's are
    addend = {"-u": "(u, v)", "-lambda": "(x, λ)", "-w": "(x, λ)"}
    if field and field[1] and name == f"pornin map x1 to {addend[field[1]]}":
        return "Pornin x1: hash to addend", False
    if str(row.family).startswith(("edwards", "twisted")) and name == "elligator2 x1":
        return "Elligator 2 x1: hash + prepare", True
    if str(row.family).startswith("weier") and name == "sswu x1":
        return "SSWU x1: hash + prepare", True
    return None


def recipe_note(recipe):
    if recipe == "hash + prepare":
        return ""
    if recipe == "hash to addend":
        return "hashed straight to $(x, \\lambda)$: no prepare"
    if recipe == "hash is the addend":
        return "hash output added as is: no prepare"
    if recipe.endswith(": hash to addend"):
        return recipe + "; output already an addend"
    return recipe


def pipeline_cost(t: pd.DataFrame, loads=DEFAULT_LOADS) -> pd.DataFrame:
    """Estimate insertion cost from operations with compatible representations.

    The standard recipe uses h2c, group.prepare, and k group.add updates at
    throughput. A direct hash-to-addend measurement replaces both hashing
    and preparation; that recipe records its preparation cost as NaN. Each
    term is exactly one benchmark of its family, construction and mode.
    Comparison maps retain distinct constructions and use only known output
    contracts; an already prepared addend needs no separate conversion.
    The references in rules.HASH_IS_ADDEND hash one item at a time and add
    the hash output, timed by the compare suite.
    Every item also takes the salted map digest and walks the mapping's
    indices below m (riblt.mapping), the same for every family; with
    either unmeasured at some load, as in runs before they were, no load
    includes them.
    Every time X has its interval in "X lo" and "X hi", and each term its benchmark in
    "<term> bench".
    """
    per = t[(t.unit == "ns/elem") & ~t.partial_batch]
    rows = []
    mapped = mapping_terms(per, loads)

    def emit(fam, hashing, recipes, add):
        for recipe, hsh, prep in recipes:
            row = {"family": fam, "hashing": hashing, "recipe": recipe}
            sums = {
                "hash_ns": [(hsh, 1)],
                "prepare_ns": [(prep, 1)],
                "add_ns": [(add, 1)],
            }
            sums |= {f"{ld.key} map_ns": mapped.get(ld.key, []) for ld in loads}
            sums |= {
                ld.key: [(hsh, 1), (prep, 1), (add, ld.k), *mapped.get(ld.key, [])]
                for ld in loads
            }
            cis = {}
            for key, terms in sums.items():
                row[key], cis[f"{key} lo"], cis[f"{key} hi"] = sum_ci(terms)
            benches = {
                f"{term} bench": described(r)
                for term, r in (("hash", hsh), ("prepare", prep), ("add", add))
            } | {
                f"{ld.key} map bench": "; ".join(
                    described(r) for r, _ in mapped[ld.key]
                )
                for ld in loads
                if ld.key in mapped
            }
            rows.append(row | cis | benches)

    for fam in t.family.cat.categories:
        f = per[per.family == fam]
        if fam in rules.HASH_IS_ADDEND:
            add = f[(f.group == "add") & (f.operation == "add")]
            add = add[add["mode"] == "throughput"]
            h = f[(f.group == "hash_to_curve") & (f["mode"] == "per-element")]
            if not (add.empty or h.empty):
                h = _one(h, f"{fam} hash")
                emit(
                    fam,
                    "one at a time",
                    [("hash is the addend", h, None)],
                    _one(add, f"{fam} add"),
                )
            continue
        add = f[(f.group == "group.add") & (f["mode"] == "throughput")]
        if add.empty:
            continue
        add = _one(add, f"{fam} add")
        for hashing, mode in (("one at a time", "per-element"), ("batched", "batch")):
            h = f[(f.group == "h2c") & (f["mode"] == mode)]
            p = f[(f.group == "group.prepare") & (f["mode"] == mode)]
            fused = f[(f.operation == "hash to addend") & (f["mode"] == mode)]
            recipes = []
            if not (h.empty or p.empty):
                h, p = _one(h, f"{fam} {mode} hash"), _one(p, f"{fam} {mode} prepare")
                recipes.append(("hash + prepare", h, p))
            if not fused.empty:
                fused = _one(fused, f"{fam} {mode} hash to addend")
                recipes.append(("hash to addend", fused, None))
            maps = f[(f.layer == rules.EXPERIMENTAL) & (f["mode"] == mode)]
            constructions = {}
            for _, mapping in maps.iterrows():
                contract = comparison_recipe(mapping)
                if contract is None:
                    continue
                label, needs_prepare = contract
                constructions.setdefault((label, needs_prepare), []).append(mapping)
            for (label, needs_prepare), measurements in constructions.items():
                mapping = _one(pd.DataFrame(measurements), f"{fam} {mode} {label}")
                if needs_prepare:
                    preparation = f[(f.group == "group.prepare") & (f["mode"] == mode)]
                    if preparation.empty:
                        continue
                    prep = _one(preparation, f"{fam} {mode} prepare")
                else:
                    prep = None
                recipes.append((label, mapping, prep))
            emit(fam, hashing, recipes, add)
    return pd.DataFrame(rows)


def mapping_terms(per: pd.DataFrame, loads) -> dict:
    """Each load's mapping terms, {load key: [(benchmark, 1)]}: the map
    digest and the walk of rules.MAPPING below its m. Empty unless every
    load has both, so that no load's estimate counts a part another's
    leaves out."""
    m = per[per.group == "riblt.mapping"]
    digest = m[m.operation == "map digest"]
    terms = {}
    for ld in loads:
        walk = m[m.full_id == f"riblt.mapping/{rules.MAPPING}/item,{ld.key}"]
        if digest.empty or walk.empty:
            return {}
        terms[ld.key] = [
            (_one(digest, "map digest"), 1),
            (_one(walk, f"{ld.key} walk"), 1),
        ]
    return terms


# The point stated in docs/workload.md (Cost in repeated reconciliation): a
# namespace of n items, d of them new
# each round, reconciled through m = 1.35 d coded symbols
ROUND_N, ROUND_D = 10**5, 10**3
ROUND_M = round(1.35 * ROUND_D)
# the costs dominance compares; the regimes' estimates are sums of them
AXES = ["add", "hash + prepare", "encode", "decode"]


class Est(NamedTuple):
    value_ns: float
    value_lo_ns: float
    value_hi_ns: float


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


def decision(pipeline, elem, fixed, chosen) -> pd.DataFrame:
    """One row per family with a batched insertion estimate: log2 rho from
    its recorded fixture, the costs AXES names, the two regimes' estimates
    at ROUND_N, ROUND_D, curve selection, and the families that dominate it.
    A missing term leaves its estimates NaN rather than smaller."""
    m, k = ROUND_M, mapping_degree(ROUND_M)
    rho = {
        rules.family(name)[0]: rho_bits(int(rec["r"]), rec["automorphisms"])
        for name, rec in fixed.items()
        if rec and "r" in rec
    }
    sel = chosen.set_index("family")
    p = pipeline[pipeline.hashing == "batched"].copy()
    if p.empty:
        return pd.DataFrame()
    p["hp"] = p.hash_ns + p.prepare_ns.fillna(0)
    rows = []
    for fam, recs in p.groupby(p.family.astype(str), sort=False):
        r = recs.loc[recs.hp.idxmin()]

        def est(key, r=r):
            return Est(r[key], r[f"{key} lo"], r[f"{key} hi"])

        hp = [(est("hash_ns"), 1)]
        benches = {"add": r["add bench"], "hash + prepare": r["hash bench"]}
        if pd.notna(r.prepare_ns):
            hp.append((est("prepare_ns"), 1))
            benches["hash + prepare"] += "; " + r["prepare bench"]
        add = est("add_ns")
        cell = {}
        for op in ("encode", "decode"):
            e = elem[(elem.curve == fam) & (elem.operation == f"{op}, batched")]
            cell[op] = (
                Est(*e.iloc[0][["value_ns", "value_lo_ns", "value_hi_ns"]])
                if len(e)
                else None
            )
            benches[op] = e.iloc[0].bench if len(e) else None
        terms = {
            "add": [(add, 1)],
            "hash + prepare": hp,
            "encode": [(cell["encode"], 1)],
            "decode": [(cell["decode"], 1)],
            "retained": [*hp, (add, k)],
            "round": [
                (add, ROUND_N * k),
                *((e, w * ROUND_D) for e, w in hp),
                (cell["encode"], m),
                (cell["decode"], m),
            ],
        }
        c = rules.curve(fam)
        row = {
            "family": fam,
            "recipe": r.recipe,
            "rho": rho.get(fam, float("nan")),
            "counting tools": rules.COUNTING_TOOLS.get(c, ""),
        }
        for key, ts in terms.items():
            if any(e is None for e, _ in ts):
                row[key] = row[f"{key} lo"] = row[f"{key} hi"] = float("nan")
            else:
                row[key], row[f"{key} lo"], row[f"{key} hi"] = sum_ci(ts)
        row |= {f"{key} bench": b for key, b in benches.items()}
        for key in ("find", "verify"):
            row[key] = sel.at[c, f"{key}_s"] * 1e9 if c in sel.index else float("nan")
        row["seeds"] = sel.at[c, "seeds"] if c in sel.index else 0
        row["method"] = sel.at[c, "method"] if c in sel.index else None
        rows.append(row)
    d = pd.DataFrame(rows).set_index("family")
    # None: not compared
    doms = dominators(d)
    d["dominated by"] = [doms.get(f) for f in d.index]
    return d


def dominators(d: pd.DataFrame) -> dict[str, list[str]]:
    """For each family, those that dominate it: log2 rho at least its own,
    and on every one of AXES either the same measurement or faster beyond
    both confidence intervals, faster on at least one. Overlapping
    intervals are no evidence either way. A family missing rho, a cost or
    an interval is not compared."""
    need = ["rho", *AXES, *(f"{ax} {e}" for ax in AXES for e in ("lo", "hi"))]
    full = d[d[need].map(lambda v: math.isfinite(v)).all(axis=1)].index

    def faster(b, a, ax):
        return d.at[b, f"{ax} hi"] < d.at[a, f"{ax} lo"]

    def shared(b, a, ax):
        return d.at[b, f"{ax} bench"] == d.at[a, f"{ax} bench"]

    return {
        a: [
            b
            for b in full
            if b != a
            and d.at[b, "rho"] >= d.at[a, "rho"]
            and all(faster(b, a, ax) or shared(b, a, ax) for ax in AXES)
            and any(faster(b, a, ax) for ax in AXES)
        ]
        for a in full
    }


def riblt_cost(t: pd.DataFrame, loads=DEFAULT_LOADS) -> pd.DataFrame:
    """Estimate hash + k additions from each family's minimum measured costs.

    This lower-bound model permits incompatible point representations and
    omits conversion costs. It is not a measured implementation runtime.
    """
    per = t[(t.unit == "ns/elem") & ~t.partial_batch]
    hashes = per[
        ((per.layer == "hash to curve") & (per.operation == "hash to curve"))
        | (per.layer == rules.EXPERIMENTAL)
    ]
    adds = per[
        (per.layer == "group ops") & (per.operation == "add") & (per["mode"] != "batch")
    ]
    rows = []
    for fam in t.family.cat.categories:
        h, a = hashes[hashes.family == fam], adds[adds.family == fam]
        if h.empty or a.empty:
            continue
        add = a.loc[a.value_ns.idxmin()]
        variants = [
            ("one at a time", h[h["mode"] != "batch"]),
            ("batched", h[h["mode"] == "batch"]),
        ]
        for hashing, hh in variants:
            if hh.empty:
                continue
            hsh = hh.loc[hh.value_ns.idxmin()]
            row = {
                "family": fam,
                "hashing": hashing,
                "hash": hsh.label,
                "hash_ns": hsh.value_ns,
                "add": add.label,
                "add_ns": add.value_ns,
            }
            row.update({ld.key: hsh.value_ns + ld.k * add.value_ns for ld in loads})
            rows.append(row)
    return pd.DataFrame(rows)


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


def insert_grid(p: pd.DataFrame, loads, hashing: str) -> Grid:
    """Per-item cost of inserting into k cells, from pipeline_cost: its
    parts, then their sums for each k. Each measured compatible map or
    direct hash-to-addend construction gets its own recipe row."""
    p = p[p.hashing == hashing]
    mode = "batched" if hashing == "batched" else "alone"
    cols = [
        Col("hash_ns", "per item, once", f"hash ({mode})"),
        Col("prepare_ns", "per item, once", f"prepare ({mode})"),
        Col("add_ns", "per cell", "add (throughput)"),
    ] + [Col(ld.key, "per item, into k cells", ld.label) for ld in loads]
    have = set(p.family.astype(str))
    groups = grouped(rules.CURVE_GROUPS, have)
    for _, rows in groups:
        for row in list(rows):
            recs = p[p.family.astype(str) == row.label]
            recs = recs.sort_values("recipe", key=lambda s: s != "hash + prepare")
            at = rows.index(row)
            for n, (_, r) in enumerate(recs.iterrows()):
                target = row
                if r.recipe != "hash + prepare":
                    note = recipe_note(r.recipe)
                    marks = row.marks | rules.OTHER_HARNESS
                    if n:
                        target = Row(row.label, note, marks=marks)
                        rows.insert(at + n, target)
                    else:
                        row.note, row.marks = note, marks
                target.values = {c.key: r[c.key] for c in cols}
                target.cis = {c.key: (r[f"{c.key} lo"], r[f"{c.key} hi"]) for c in cols}
                target.tips = {
                    f"{term}_ns": r[f"{term} bench"]
                    for term in ("hash", "prepare", "add")
                    if isinstance(r[f"{term} bench"], str)
                } | {
                    ld.key: r[f"{ld.key} map bench"]
                    for ld in loads
                    if isinstance(r.get(f"{ld.key} map bench"), str)
                }
                if pd.isna(r.prepare_ns):
                    target.values["prepare_ns"] = "in the hash"
    return Grid(cols, groups, corner="curve")


def _tex_int(n: int) -> str:
    # in math a bare comma is punctuation and spaces the digits after it
    return f"{n:,}".replace(",", "{,}")


def dominance_label(dom: list[str] | None) -> str:
    return "not compared" if dom is None else ", ".join(dom) or "none"


def decision_grid(d: pd.DataFrame) -> Grid:
    """decision's rows, grouped by field as the other tables are."""
    k = mapping_degree(ROUND_M)
    cols = [
        Col("rho", "security", "$\\log_2$ rho", PLAIN),
        Col("add", "per cell", "add"),
        Col("encode", "per cell", "encode"),
        Col("decode", "per cell", "decode"),
        Col("hash + prepare", "per new item", "hash + prepare"),
        Col("retained", "per new item", f"sums retained, $k = {k:.2f}$"),
        Col(
            "round",
            "per party and round",
            f"sums rebuilt, $n = {_tex_int(ROUND_N)}$, $d = {_tex_int(ROUND_D)}$",
        ),
        # means over seeds, with no interval: shown, never marked best or tied
        Col("find", "curve selection", "find", PLAIN_TIME),
        Col("verify", "curve selection", "verify", PLAIN_TIME),
        Col("method", "curve selection", "measured method", PLAIN),
        Col("counting tools", "curve selection", "available tools", PLAIN),
        Col("dominated by", "dominated by", kind=PLAIN),
    ]
    groups = grouped(rules.CURVE_GROUPS, set(d.index))
    for _, rows in groups:
        for row in rows:
            r = d.loc[row.label]
            row.values = {c.key: r[c.key] for c in cols}
            row.values["dominated by"] = dominance_label(r["dominated by"])
            row.cis = {
                c.key: (r[f"{c.key} lo"], r[f"{c.key} hi"])
                for c in cols
                if f"{c.key} lo" in r
            }
            if r.recipe != "hash + prepare":
                row.note += "; " + recipe_note(r.recipe)
            row.tips = {
                ax: r[f"{ax} bench"] for ax in AXES if isinstance(r[f"{ax} bench"], str)
            }
            # the summed estimates name their terms' benchmarks too
            bench = {ax: row.tips.get(ax, "not measured") for ax in AXES}
            row.tips["retained"] = f"{bench['hash + prepare']}; k × ({bench['add']})"
            row.tips["round"] = (
                f"n k × ({bench['add']}); d × ({bench['hash + prepare']}); "
                f"m × ({bench['encode']}; {bench['decode']})"
            )
            if r.seeds:
                row.tips["find"] = f"{r.method}, mean of {r.seeds:g} seeds"
                row.tips["verify"] = row.tips["find"]
    return Grid(cols, groups, corner="curve")


def batching_grid(p: pd.DataFrame, e: pd.DataFrame) -> Grid:
    """Compare per-item batching savings for hashing, preparation, and codecs.

    Normalize each saving by the measured per-item inversion saving for
    that field. The ratio is a timing diagnostic, not an operation count.
    """
    inv = e[(e.layer == "field") & (e.facet == "field invert")]
    inv = {(r.curve, r.variant): r.value_ns for r in inv.itertuples()}
    p = p[p.recipe == "hash + prepare"]
    one = p[p.hashing == "one at a time"].set_index("family")
    bat = p[p.hashing == "batched"].set_index("family")
    codec = e[(e.layer != "field") & e.facet.isin(["encode", "decode"])]
    codec = {(r.curve, r.facet, r.variant): r.value_ns for r in codec.itertuples()}
    steps = [
        ("hash_ns", "hash to curve"),
        ("prepare_ns", "prepare"),
        ("encode", "encode"),
        ("decode", "decode"),
    ]

    def times(step, family):
        if step in ("encode", "decode"):
            a = codec.get((family, step, "one at a time"))
            b = codec.get((family, step, "batched"))
            return a, b
        if family in one.index and family in bat.index:
            return one.loc[family, step], bat.loc[family, step]
        return None, None

    have = set(one.index.astype(str)) & set(bat.index.astype(str))
    have |= {c for c, _, _ in codec}
    groups = grouped(rules.CURVE_GROUPS, have)
    cols = []
    for step, head in steps:
        if not any(
            None not in times(step, r.label) for _, rows in groups for r in rows
        ):
            continue
        cols += [
            Col((step, "alone"), head, "alone", PLAIN_TIME),
            Col((step, "batched"), head, "batched"),
            Col((step, "factor"), head, "faster by", PLAIN_FACTOR),
            Col((step, "inversions"), head, "inversions shared", PLAIN),
        ]
    for _, rows in groups:
        for row in rows:
            f = rules.FIELD_OF.get(row.label)
            i, ib = inv.get((f, "one at a time")), inv.get((f, "batched"))
            for step, _ in steps:
                a, b = times(step, row.label)
                if a is None or b is None:
                    continue
                row.values[(step, "alone")] = a
                row.values[(step, "batched")] = b
                row.values[(step, "factor")] = a / b
                if i and ib and i > ib:
                    n = (a - b) / (i - ib)
                    row.values[(step, "inversions")] = round(n, 1)
                    row.tips[(step, "inversions")] = (
                        f"saves {fmt_time(a - b)} per item; a {f} inversion alone "
                        f"is {fmt_time(i)}, batched {fmt_time(ib)} per item"
                    )
    return Grid(cols, groups, corner="curve")


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


def cost_table(c: pd.DataFrame) -> pd.DataFrame:
    c = c.copy()
    for col in c.columns:
        if col.endswith("_ns") or col.startswith(("k=", "m=")):
            # a NaN prepare: the hash's output is already the addend
            c[col] = ["in the hash" if pd.isna(v) else fmt_time(v) for v in c[col]]
    return c.rename(
        columns={"hash_ns": "hash", "prepare_ns": "prepare", "add_ns": "add"}
        if "prepare_ns" in c
        else {"hash_ns": "hash time", "add_ns": "add time"}
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
    costs,
    loads,
    drawn,
    skipped,
    source,
    elem=None,
    elem_fig=None,
    meta=None,
    pipeline=None,
    insert_fig=None,
    decided=None,
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
    if decided is not None and len(decided):
        b.extend(decision_blocks(decided))
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
    roles = (meta or {}).get("group_plan") or {}
    if roles:
        b.append(("h2", "Group benchmark plan"))
        b.append(
            (
                "p",
                (
                    "Which families the group suite's second pass timed. Its "
                    "first pass timed every family's core rows; from them and "
                    "the comparison maps, bench-run chose the families whose "
                    "modeled insertion at $m = 1350$ was within 50% of the "
                    "cheapest (lead) and, if none of them was over a field of "
                    "odd characteristic, the cheapest that was (contrast). The "
                    "second pass added the batch sweep, subtraction and equality "
                    "on pure cells "
                    "([[methodology.md#the-group-plan|The group plan]])."
                ),
            )
        )
        b.append(
            ("table", pd.DataFrame(list(roles.items()), columns=["family", "role"]))
        )
    plan = (meta or {}).get("riblt_plan") or {}
    if plan:
        b.append(("h2", "RIBLT benchmark plan"))
        b.append(
            (
                "p",
                (
                    "Which RIBLT benchmarks this run generated. Before the RIBLT "
                    "suite ran, bench-run chose from this run's group and "
                    "hash-to-curve results, for each family below, the hash its "
                    "RIBLT benchmarks use and their scope: every sweep (full), each "
                    "axis's endpoints and middle (buffer), or one point of each "
                    "benchmark (spot). The plan records what was measured, not a "
                    "result; the XOR checksums and ristretto255 are timed in full "
                    "regardless ([[methodology.md#the-riblt-plan|The RIBLT plan]])."
                ),
            )
        )
        b.append(
            (
                "table",
                pd.DataFrame(
                    [(f, v.get("h2c"), v.get("scope")) for f, v in plan.items()],
                    columns=["family", "hash", "scope"],
                ),
            )
        )
    if pipeline is not None and len(pipeline):
        b.extend(insert_blocks(pipeline, loads, insert_fig))
        if elem is not None and len(elem):
            b.extend(batching_blocks(pipeline, elem, batch_size_grid(t)))
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
    if len(costs):
        b.extend(lower_bound_blocks(costs, loads))
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


def decision_blocks(d):
    m = ROUND_M
    what = (
        "Each family's nominal $\\log_2$ rho, this run's costs (addition at "
        "throughput; hashing, preparation, encoding and decoding batched), the "
        "two regimes of [[workload.md#cost-in-repeated-reconciliation|Cost in "
        f"repeated reconciliation]] for a set of $n = 10^5$ items, a set difference "
        f"of $d = 10^3$ and a sketch of $m = {m}$ cells, in which a new item updates "
        f"$k = {mapping_degree(m):.2f}$ cells on average, and "
        "the mean time over seeds to find a curve and verify its certificate. "
        "Dominance is as defined in "
        "[[methodology.md#measurement-methodology|Measurement methodology]]."
    )
    share = (ROUND_N * mapping_degree(m) * d["add"] / d["round"]).dropna()
    if len(share):
        what += (
            f" Additions are {100 * share.min():.0f} to {100 * share.max():.0f}% "
            "of each rebuilt estimate."
        )
    if d["rho"].isna().all():
        what += " No fixture's $r$ was recorded, so no family is compared on security."
    return [("h2", "Families compared"), ("p", what), ("grid", decision_grid(d))]


def describe_loads(loads) -> str:
    return "; ".join(ld.label for ld in loads)


def insert_blocks(pipeline, loads, fig=None):
    what = (
        "Modeled insertion per item into a sketch of $m$ cells: a hash, a "
        "prepare, and an addition to each of the $k$ cells the item maps to on "
        f"average ({describe_loads(loads)}), by compatible recipes "
        "([[methodology.md#reading-a-run-report|Reading a run report]])."
    )
    refs = [f for f in rules.HASH_IS_ADDEND if f in set(pipeline.family.astype(str))]
    if refs:
        what += (
            f" The references ({', '.join(refs)}) have no batched hash: they hash "
            "one item at a time, add the hash output as is, and are in the figure "
            "and the second grid."
        )
    maps = [
        ld
        for ld in loads
        if f"{ld.key} map_ns" in pipeline and pipeline[f"{ld.key} map_ns"].notna().any()
    ]
    if maps:
        cost = ", ".join(
            f"{fmt_time(pipeline[f'{ld.key} map_ns'].iloc[0])} at $m = {ld.key[2:]}$"
            for ld in maps
        )
        what += (
            " Each total also counts, the same for every family, the salted "
            f"SHA-256 map digest and the {rules.MAPPING} walk through the item's "
            f"indices below $m$ (riblt.mapping): {cost}. The key XOR and count "
            "of each cell are not modelled; riblt.cells measures them with the "
            "rest of the cell updates."
        )
    else:
        what += (
            " The map digest and the mapping's walk, the same for every family, "
            "were not timed at every $m$ in this run and are not included."
        )
    b = [("h2", "Insertion"), ("p", what)]
    b.append(("grid", insert_grid(pipeline, loads, "batched")))
    if fig:
        b.append(("img", fig))
    b.append(
        (
            "details",
            "Hashing and preparing one item at a time",
            [("grid", insert_grid(pipeline, loads, "one at a time"))],
        )
    )
    return b


def batching_blocks(pipeline, elem, sweep=None):
    what = (
        "Batched against individual time per element. The inversions-shared "
        "column divides the saving by the difference between individual and "
        "batched field inversion; it also absorbs allocation and control flow."
    )
    b = [
        ("h2", "Batching"),
        ("p", what),
        ("grid", batching_grid(pipeline, elem)),
    ]
    if sweep is not None:
        by_size = "Time per element against batch size."
        b += [("p", by_size), ("grid", sweep)]
    return b


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


def lower_bound_blocks(costs, loads):
    what = (
        "One hash plus an addition to each of the $k$ cells an item maps to in "
        f"a sketch of $m$ cells ({describe_loads(loads)}), from each family's "
        "fastest hash and addition in any representation, without conversions; "
        "the only estimate for families without group-trait benchmarks."
    )
    return [
        (
            "details",
            "Lower bound, every family",
            [("p", what), ("table", cost_table(costs))],
        )
    ]


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
    root: Path, out: Path, loads=DEFAULT_LOADS, stat="typical", formats=("png", "svg")
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
    costs = riblt_cost(t, loads)
    pipeline = pipeline_cost(t, loads)
    decided = (
        decision(pipeline, elem, fixtures(t, meta), selection(root))
        if len(pipeline)
        else None
    )
    insert_fig = (
        plot_insert(pipeline, loads, out, formats, footer) if len(pipeline) else None
    )
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
    pipeline.to_csv(out / "riblt_cost.csv", index=False)
    costs.to_csv(out / "riblt_cost_lower_bound.csv", index=False)
    if decided is not None and len(decided):
        labels = decided["dominated by"].map(dominance_label)
        decided.assign(**{"dominated by": labels}).to_csv(out / "decision.csv")
    coverage(t).to_csv(out / "coverage.csv")
    if meta:
        (out / "meta.json").write_text(published(meta))
    bs = blocks(
        t,
        costs,
        loads,
        drawn,
        skipped,
        run_name(root, meta),
        elem,
        elem_fig,
        meta,
        pipeline,
        insert_fig,
        decided,
        machine_load(root),
        runs=runs,
        selection_fig=selection_fig,
    )
    bs.append(("footer", footer or run_name(root, meta)))
    for s in skipped:
        print(f"skipped {s}", file=sys.stderr)
    return t, bs


def report(
    root: Path, out: Path, loads=DEFAULT_LOADS, stat="typical", formats=("png", "svg")
) -> pd.DataFrame:
    t, bs = analyse(root, out, loads, stat, formats)
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
        "--m",
        default=",".join(map(str, DEFAULT_MS)),
        help="RIBLT coded-symbol counts, each costing k(m) adds per item",
    )
    p.add_argument("--k", help="adds per item instead, for exploration")
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
    p.add_argument(
        "--riblt-plan",
        action="store_true",
        help="read benchmark ids from stdin and print how benches/riblt.rs is to "
        "time each RIBLT family, from the group operations CRITERION measured: "
        "lines of family, h2c parameter and scope, tab-separated (RIBLT_PLAN)",
    )
    p.add_argument(
        "--group-plan",
        action="store_true",
        help="read benches/group.rs ids from stdin and print the families its "
        "second pass is to time, from the first pass CRITERION measured: lines "
        "of family and role, tab-separated (GROUP_PLAN)",
    )
    a = p.parse_args(argv)
    if a.riblt_plan or a.group_plan:
        flag = "--riblt-plan" if a.riblt_plan else "--group-plan"
        if a.criterion is None:
            p.error(f"{flag} reads the measurements under CRITERION")
        t, _ = load(a.criterion)
        if not len(t):
            return
        if a.riblt_plan:
            for fam, (token, scope) in riblt_plan(tidy(t), sys.stdin).items():
                print(f"{fam}\t{token}\t{scope}")
        else:
            for fam, role in group_plan(tidy(t), sys.stdin).items():
                print(f"{fam}\t{role}")
        return
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
    loads = (
        fixed_loads(float(k) for k in a.k.split(","))
        if a.k
        else riblt_loads(int(m) for m in a.m.split(","))
    )
    report(a.criterion, a.out, loads, a.stat, tuple(a.formats.split(",")))


if __name__ == "__main__":
    main()
