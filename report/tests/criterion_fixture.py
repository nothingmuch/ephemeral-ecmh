"""Write a criterion 0.8 output tree from fixture rows.

  python criterion_fixture.py fixtures/2026-09-30.tsv DIR

Mirrors criterion's report.rs (BenchmarkId, make_filename_safe) and
estimate.rs (Estimates) serde output, so the report reads what criterion
writes. Rows without an element count get no slope, as with flat sampling.
"""

import csv
import json
import sys
from pathlib import Path


def filename_safe(s: str) -> str:
    for c in '?"/\\*<>:|^':
        s = s.replace(c, "_")
    return s.encode()[:64].decode(errors="ignore")  # 64 bytes, on a char boundary


def _estimate(point, lo, hi):
    return {
        "confidence_interval": {
            "confidence_level": 0.95,
            "lower_bound": lo,
            "upper_bound": hi,
        },
        "point_estimate": point,
        "standard_error": (hi - lo) / 4,
    }


def write_bench(
    root: Path,
    group,
    function=None,
    parameter=None,
    lo=1.0,
    point=1.0,
    hi=1.0,
    throughput=None,
    slope=True,
):
    parts = [p for p in (group, function, parameter) if p]
    d = root.joinpath(*map(filename_safe, parts))
    for sub in ("new", "base"):  # base: criterion keeps the previous run too
        (d / sub).mkdir(parents=True, exist_ok=True)
        bench = {
            "group_id": group,
            "function_id": function,
            "value_str": parameter,
            "throughput": throughput,
            "full_id": "/".join(parts),
            "directory_name": "/".join(map(filename_safe, parts)),
            "title": "/".join(parts),
        }
        (d / sub / "benchmark.json").write_text(json.dumps(bench))
        # the base run is 2x slower: a report that reads it is visibly wrong
        k = 1 if sub == "new" else 2
        est = _estimate(point * k, lo * k, hi * k)
        spread = _estimate((hi - lo) * k, 0.0, (hi - lo) * 2 * k)
        estimates = {
            # with a slope, the mean is off it, so a test can tell which was read
            "mean": _estimate(point * 1.01 * k, lo * 1.01 * k, hi * 1.01 * k)
            if slope
            else est,
            "median": est,
            "median_abs_dev": spread,
            "slope": est if slope else None,
            "std_dev": spread,
        }
        (d / sub / "estimates.json").write_text(json.dumps(estimates))
    return d


def read_rows(tsv: Path) -> list[dict]:
    """The rows of a fixture, after the `#` lines that say where it comes from."""
    lines = [s for s in tsv.read_text().splitlines() if not s.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def write_tree(rows, root: Path) -> Path:
    for r in rows:
        elements = int(r["elements"]) if r["elements"] else None
        write_bench(
            root,
            r["group"],
            r["function"] or None,
            r["parameter"] or None,
            float(r["lower_ns"]),
            float(r["point_ns"]),
            float(r["upper_ns"]),
            {"Elements": elements} if elements else None,
            slope=elements is not None,
        )
    return root


if __name__ == "__main__":
    write_tree(read_rows(Path(sys.argv[1])), Path(sys.argv[2]))
