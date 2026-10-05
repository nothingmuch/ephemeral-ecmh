"""Selection-stage costs in units of same-run measured insertion work.

``curvegen/certify`` already measures search and certificate construction;
``curvegen/find`` is an alternative search-only measurement, not another term.
PARI's ``curvegen/count`` times one accepted candidate, not curve selection.
Ranges below are extrema of observed seed point estimates, not confidence
intervals or break-even counts against an alternative implementation.
"""

import math
import re

import figures
import pandas as pd

FAMILIES = figures.REPRESENTATIVES
# Fixed groups: no curve is selected for them, so their selection cost is
# zero by construction rather than unmeasured.
BASELINES = (figures.BASELINE, "ristretto255")


def _positive(value):
    try:
        return math.isfinite(float(value)) and float(value) > 0
    except (TypeError, ValueError):
        return False


def _seed(value):
    """Criterion's fixture index, accepting numeric columns from saved CSVs."""
    if pd.isna(value):
        return None
    text = str(value)
    if re.fullmatch(r"\d+(?:\.0)?", text):
        return str(int(float(text)))
    return None


def default_map(table):
    """Rows under the default mapping and the cheapest hash. A `map=`
    parameter (benches/riblt.rs) isolates the mapping's share of a family's
    cost, and a `proj=` parameter keys its cells by IDs under a projection;
    neither is that family's measurement of the workload."""
    default = ~table.full_id.str.contains(r"(?:^|[/,])(?:map|proj)=", regex=True)
    # Under several hash constructions (an unplanned run, benches/riblt.rs
    # RIBLT_PLAN), a family's workload is timed under the cheapest
    hashed = r",?\bh2c=[^,/]+"
    varies = table.full_id.str.contains(hashed, regex=True)
    key = table.full_id.str.replace(hashed, "", regex=True)
    cost = table.value_ns.where(default & varies)
    return default & (~varies | cost.eq(cost.groupby(key).transform("min")))


def _insertion(table, family):
    rows = table[
        (table.group == "riblt.encode") & (table.family == family) & default_map(table)
    ]
    rows = rows[
        rows.full_id.map(
            lambda name: (
                dict(re.findall(r"\b(m|n)=(\d+)\b", name)) == {"m": "1350", "n": "3500"}
            )
        )
    ]
    # There is no justified choice between two measurements of this workload.
    if len(rows) != 1:
        return None
    row = rows.iloc[0]
    if row.unit != "ns/elem" or row.elements != 3500 or not _positive(row.value_ns):
        return None
    return row


def _record(run, family, method, stage, complete, rows, insert, seed_ids):
    values = [float(r.value_ns) for r in rows]
    cost = float(insert.value_ns)
    return {
        "run": run.name,
        "family": family,
        "method": method,
        "stage": stage,
        "complete": complete,
        "lo_items": min(values) / cost,
        "hi_items": max(values) / cost,
        "seeds": len(rows),
        "seed_ids": seed_ids,
        "raw_ids": [r.full_id for r in rows],
        "stage_ns": values,
        "insert_id": insert.full_id,
        "insert_ns_per_item": cost,
        "range_kind": "observed seed point estimates"
        if seed_ids
        else "one candidate fixture",
    }


def selection(runs) -> list[dict]:
    """Extract S6 records from runs exposing ``name``, ``meta`` and ``table``.

    Tables are ``bench_report.tidy`` results. Match only the exact family in
    the same run to measured ``riblt.encode`` at m=1350, n=3500. Prefer the
    measured certify stage over find; never fill missing certify seeds with
    find times. Keep only unambiguous, positive, finite stage groups. Missing
    measurements produce no record, and observed seed coverage is explicit.
    Count-only records have one fixture, identified by their raw benchmark;
    their ``seeds=1`` is not a claim of a measured search over that seed.
    """
    result = []
    for run in runs:
        table = run.table
        required = {
            "group",
            "family",
            "function",
            "parameter",
            "full_id",
            "unit",
            "elements",
            "value_ns",
        }
        if table.empty or not required.issubset(table.columns):
            continue
        for family in FAMILIES:
            insert = _insertion(table, family)
            if insert is None:
                continue
            rows = table[table.family == family]
            certify = rows[rows.group == "curvegen/certify"]
            stage_rows = (
                certify if len(certify) else rows[rows.group == "curvegen/find"]
            )
            # These are the complete seed-search measurements defined in agm.rs.
            if family == "gf2_127" and len(stage_rows):
                records = list(stage_rows.itertuples(index=False))
                seed_ids = [_seed(r.parameter) for r in records]
                valid = all(
                    r.unit == "ns/iter"
                    and _positive(r.value_ns)
                    and isinstance(r.function, str)
                    and re.fullmatch(r"gf2_127 agm\+sieve/index \d+", r.function)
                    for r in records
                )
                if (
                    valid
                    and None not in seed_ids
                    and len(set(seed_ids)) == len(seed_ids)
                ):
                    pairs = sorted(
                        zip(seed_ids, records), key=lambda pair: int(pair[0])
                    )
                    seed_ids, records = map(list, zip(*pairs))
                    complete = bool(len(certify))
                    result.append(
                        _record(
                            run,
                            family,
                            "agm+sieve",
                            "search + certificate"
                            if complete
                            else "search only (partial)",
                            complete,
                            records,
                            insert,
                            seed_ids,
                        )
                    )
            counts = rows[rows.group == "curvegen/count"]
            if len(counts) == 1:
                count = counts.iloc[0]
                if count.unit == "ns/iter" and _positive(count.value_ns):
                    result.append(
                        _record(
                            run,
                            family,
                            "PARI",
                            "one candidate point count (partial)",
                            False,
                            [count],
                            insert,
                            [],
                        )
                    )
    return result


def references(runs) -> list[dict]:
    """The BASELINES' measured insertion per item, matched as in selection():
    with no selection stage, it is their whole cost per item at any volume."""
    result = []
    for run in runs:
        table = run.table
        if table.empty or "full_id" not in table.columns:
            continue
        for family in BASELINES:
            insert = _insertion(table, family)
            if insert is not None:
                result.append(
                    {
                        "run": run.name,
                        "family": family,
                        "insert_id": insert.full_id,
                        "insert_ns_per_item": float(insert.value_ns),
                    }
                )
    return result
