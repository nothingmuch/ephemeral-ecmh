"""Compare matched target configurations; ratios do not identify causal instructions."""

import math
import re
import shlex

import pandas as pd

FIELDS = {"gf2_109", "gf2_122", "gf2_127", "fp127"}
CURVES = {
    base + variant
    for base in ("gf2_109", "gf2_122-gls", "gf2_127")
    for variant in ("", "-lambda", "-u", "-w")
} | {"edwards127"}
# the references' additions, at throughput as the curves' are; a timing
# with one accumulator measures latency instead
REFERENCES = {"ristretto255": "point add throughput", "xor": "XOR add throughput"}
IDENTITY = ("cpu", "arch", "commit", "rustc", "profile")
ROW_IDENTITY = (
    "group",
    "function",
    "parameter",
    "stat",
    "family",
    "layer",
    "operation",
    "mode",
    "unit",
    "elements",
    "batch_n",
)


def _flags(meta):
    """Separate the declared CPU flag; other compiler flags must agree."""
    target = meta.get("target")
    if not target or not isinstance(meta.get("rustflags"), str):
        return None
    try:
        tokens = shlex.split(meta["rustflags"])
    except ValueError:
        return None
    cpus, rest = [], []
    i = 0
    while i < len(tokens):
        token = tokens[i]
        if token == "-C" and i + 1 < len(tokens):
            i += 1
            token += tokens[i]
        if token.startswith("-Ctarget-cpu="):
            cpus.append(token.split("=", 1)[1])
        else:
            rest.append(token)
        i += 1
    return tuple(rest) if cpus == [target] else None


def _compatible(native, generic):
    n, g = native.meta, generic.meta
    if n.get("dirty") is not False or g.get("dirty") is not False:
        return False
    if any(not n.get(k) or n.get(k) != g.get(k) for k in IDENTITY):
        return False
    cpu_target = re.sub(r"[^a-z0-9]+", "-", n["cpu"].lower()).strip("-")
    if n.get("target") not in ("native", cpu_target):
        return False
    # Older runs omit host IDs: use their recorded CPU/architecture, but never
    # pair two different explicit hosts or different recorded OS environments.
    if any(n.get(k) != g.get(k) for k in ("host", "os", "kernel")):
        return False
    nf, gf = _flags(n), _flags(g)
    return nf is not None and nf == gf


def _operation(row):
    if row.get("unit") != "ns/elem":
        return None
    if row.get("layer") == "field" and row.get("family") in FIELDS:
        if row.get("operation") == "mul" and row.get("mode") == "throughput":
            return "field mul throughput"
        if row.get("operation") == "invert":
            return "field invert"
    if row.get("family") in REFERENCES and (
        row.get("group"),
        row.get("operation"),
        row.get("mode"),
    ) == ("add", "add", "throughput"):
        return REFERENCES[row["family"]]
    if row.get("family") in CURVES:
        if row.get("group") == "group.add" and row.get("mode") == "throughput":
            return "point add throughput"
        if (
            row.get("group") == "group.prepare"
            and row.get("mode") == "batch"
            and row.get("batch_n") == 1024
        ):
            return "point prepare batch 1024"
    return None


def _rows(table):
    rows, duplicates = {}, set()
    for row in table.to_dict("records"):
        if any(pd.isna(row.get(k)) for k in ("group", "full_id", "stat", "elements")):
            continue
        operation = _operation(row)
        if not operation:
            continue
        key = tuple(None if pd.isna(row.get(k)) else row.get(k) for k in ROW_IDENTITY)
        if key in rows:
            duplicates.add(key)
        rows[key] = (operation, row)
    return {key: value for key, value in rows.items() if key not in duplicates}


def _valid(row):
    values = [row.get(k) for k in ("value_ns", "value_lo_ns", "value_hi_ns")]
    return (
        all(isinstance(v, (int, float)) and math.isfinite(v) and v > 0 for v in values)
        and values[1] <= values[2]
    )


def compare(runs):
    """Return generic/native time ratios for exact benchmark matches.

    `lo` and `hi` divide the marginal interval endpoints; they are bounds
    from the input intervals, not an independently estimated ratio CI.
    Unknown required provenance and ambiguous duplicate measurements are
    excluded. No minimum across implementations or batch sizes is selected.
    """
    runs = list(runs)
    native = [r for r in runs if r.meta.get("target") not in (None, "", "generic")]
    generic = [r for r in runs if r.meta.get("target") == "generic"]
    result = []
    for n in native:
        for g in generic:
            if not _compatible(n, g):
                continue
            nr, gr = _rows(n.table), _rows(g.table)
            for key, (operation, a) in nr.items():
                if key not in gr:
                    continue
                _, b = gr[key]
                if not (_valid(a) and _valid(b)):
                    continue
                result.append(
                    {
                        "run": f"{n.name} / {g.name}",
                        "family": a["family"],
                        "operation": operation,
                        "value": b["value_ns"] / a["value_ns"],
                        "lo": b["value_lo_ns"] / a["value_hi_ns"],
                        "hi": b["value_hi_ns"] / a["value_lo_ns"],
                        "native_id": a["full_id"],
                        "generic_id": b["full_id"],
                        "native_run": n.name,
                        "generic_run": g.name,
                        "native_source": n.meta["commit"],
                        "generic_source": g.meta["commit"],
                        "source": n.meta["commit"],
                        "native_target": n.meta["target"],
                        "generic_target": g.meta["target"],
                    }
                )
    return result
