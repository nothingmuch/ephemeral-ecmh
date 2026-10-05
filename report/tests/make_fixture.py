"""Build a fixture TSV from summarized bench results and criterion's raw log.

  python make_fixture.py OUT.tsv PARSED.txt [PARSED.txt ...] [--log RAW.txt]

PARSED lines are "<full id>  <value> <unit>[ (total)]", per element unless
"(total)". The raw log (cargo bench's stdout) supplies what the summary lost:
the confidence interval and the element count (time / throughput). Ids
missing from the log get the bench suite's N elements and a +-0.5% interval.
The first file that has an id wins.
"""

import argparse
import csv
import re
from pathlib import Path

UNITS = {"ps": 1e-3, "ns": 1.0, "µs": 1e3, "us": 1e3, "ms": 1e6, "s": 1e9}
ELEM = {"elem/s": 1.0, "Kelem/s": 1e3, "Melem/s": 1e6, "Gelem/s": 1e9}
# criterion's full id does not say which "/" splits function from parameter:
# only bench_with_input groups have a parameter, the last component.
PARAM_GROUPS = {"verify_certificate"}
N = 1024
FIELDS = [
    "group",
    "function",
    "parameter",
    "lower_ns",
    "point_ns",
    "upper_ns",
    "elements",
    "summary_ns",
]


def parse_summary(path: Path) -> dict[str, tuple[float, bool]]:
    out = {}
    for line in path.read_text().splitlines():
        m = re.match(r"^(\S.*?)\s{2,}([\d.]+) (\S+)( \(total\))?$", line)
        if m:
            out.setdefault(m[1], (float(m[2]) * UNITS[m[3]], bool(m[4])))
    return out


def parse_log(path: Path, ids) -> dict[str, dict]:
    out, cur = {}, None
    for line in path.read_text().splitlines():
        for i in ids:
            if line.rstrip() == i or re.match(re.escape(i) + r"\s+time:", line):
                cur = out.setdefault(i, {})
        if cur is None:
            continue
        m = re.search(r"(time|thrpt):\s+\[(\S+) (\S+) (\S+) (\S+) (\S+) (\S+)\]", line)
        if m and m[1] not in cur:
            vals = [float(m[k]) for k in (2, 4, 6)]
            units = (m[3], m[5], m[7])
            if m[1] == "time":
                cur["time"] = [v * UNITS[u] for v, u in zip(vals, units)]
            elif all(u in ELEM for u in units):
                cur["thrpt"] = [v * ELEM[u] for v, u in zip(vals, units)]
        if "change:" in line:
            cur = None  # the next "time:" is a relative change
    return out


def rows(summaries, log):
    seen = {}
    for s in summaries:
        for i, v in s.items():
            seen.setdefault(i, v)
    for i, (value, total) in seen.items():
        group, rest = i.split("/", 1)
        function, parameter = rest, ""
        if group in PARAM_GROUPS:
            function, parameter = rest.rsplit("/", 1)
        r = log.get(i, {})
        if total:
            elements = ""
        elif "thrpt" in r:
            elements = round(r["time"][1] * 1e-9 * r["thrpt"][1])
        else:
            elements = N
        lo, pt, hi = r.get("time") or [
            value * (elements or 1) * f for f in (0.995, 1, 1.005)
        ]
        yield dict(
            zip(FIELDS, [group, function, parameter, lo, pt, hi, elements, value])
        )


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("out", type=Path)
    p.add_argument("summary", type=Path, nargs="+")
    p.add_argument("--log", type=Path)
    a = p.parse_args()
    summaries = [parse_summary(s) for s in a.summary]
    ids = {i for s in summaries for i in s}
    log = parse_log(a.log, ids) if a.log else {}
    with a.out.open("w", newline="") as f:
        w = csv.DictWriter(f, FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows(summaries, log):
            w.writerow(
                {k: f"{v:.6g}" if isinstance(v, float) else v for k, v in r.items()}
            )


if __name__ == "__main__":
    main()
