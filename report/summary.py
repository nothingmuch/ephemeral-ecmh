"""Cross-run views of recorded observations, with explicit coverage and provenance.

Selection never pools runs or turns missing measurements into zero. The CSV
beside the figures keeps every displayed estimate and its benchmark identity.
"""

import html
import itertools
import re
import textwrap
import warnings
from dataclasses import dataclass
from pathlib import Path

import bench_report as br
import figures
import matplotlib.pyplot as plt
import pandas as pd
import rules
import summary_pairs
import summary_selection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import FuncFormatter, NullFormatter

# every view follows the same families, so that one can be read across them
FAMILIES = (figures.BASELINE, *figures.REPRESENTATIVES)
BASELINES = summary_selection.BASELINES


@dataclass
class Run:
    name: str
    meta: dict
    table: pd.DataFrame

    @classmethod
    def read(cls, path: Path):
        meta = br.machine(path)
        if not meta or not meta.get("commit") or meta.get("dirty") is not False:
            raise ValueError(f"{path.name}: summary requires recorded, clean source")
        raw, skipped = br.load(path)
        if skipped or raw.empty:
            raise ValueError(f"{path.name}: incomplete or empty measurements")
        t = br.tidy(raw)
        if len(br.unclassified(t)):
            raise ValueError(f"{path.name}: unclassified benchmark identities")
        return cls(meta.get("name") or path.name, meta, t)


def parameter(table, key, value):
    return table.full_id.str.contains(
        rf"(?:^|[/,]){re.escape(key)}={value}(?:,|$)", regex=True
    )


def subset(table, representatives=figures.REPRESENTATIVES):
    """Three lowest observations plus boundary overlaps and fixed references.

    Overlap retains coverage; it is not a statistical equivalence test.
    """
    d = table[~table.family.astype(str).str.startswith("xor")].sort_values("value_ns")
    if d.empty:
        return []
    first = d.head(figures.TOP)
    edge = first.iloc[-1]
    overlap = d[
        (d.value_lo_ns <= edge.value_hi_ns) & (d.value_hi_ns >= edge.value_lo_ns)
    ]
    keep = list(
        dict.fromkeys(
            [*figures.shown(d.set_index("family")), *overlap.family.astype(str)]
        )
    )
    keep += [
        f for f in representatives if f in set(d.family.astype(str)) and f not in keep
    ]
    return keep


def record(run, r, **extra):
    return dict(
        run=run.name,
        family=str(r.family),
        value=float(r.value_ns),
        lo=float(r.value_lo_ns),
        hi=float(r.value_hi_ns),
        benchmark=r.full_id,
        confidence=float(r.confidence),
        **extra,
    )


def _unique(d, description):
    if d.family.duplicated().any():
        raise ValueError(f"ambiguous {description}: multiple measurements per family")
    return d


def insertion(runs):
    by_run = []
    families = []
    for run in runs:
        t = run.table
        d = t[
            (t.group == "riblt.encode")
            & (t.unit == "ns/elem")
            & parameter(t, "m", 1350)
            & parameter(t, "n", 3500)
            & summary_selection.default_map(t)
        ]
        d = _unique(d, f"{run.name} insertion")
        keep = subset(d)
        families += [f for f in keep if f not in families]
        by_run.append((run, d))
    data = []
    for run, d in by_run:
        for _, r in d[d.family.isin([*families, "xor"])].iterrows():
            data.append(
                record(
                    run, r, available=len(d), selected=len(d[d.family.isin(families)])
                )
            )
    return data, families


SCALING = [
    ("encoding n = 3500 items, per item", "sketch cells m"),
    ("peeling, per difference", "set difference d (sketch of m = 2d cells)"),
    ("hashing to the curve, per item", "batch size (1: one at a time)"),
]
MAPS = [
    ("try-and-increment", "try-and-increment", "o"),
    ("pornin", "Pornin's map", "s"),
    ("elligator2", "Elligator 2", "^"),
    ("sswu", "simplified SWU", "v"),
]


def _map(function):
    return next((name for key, name, _ in MAPS if key in str(function)), None)


def _number(parameter, key):
    m = re.search(rf"(?:^|,){key}=(\d+)", str(parameter))
    return float(m.group(1)) if m else None


def scaling(runs, families):
    """Each run's measurements of families against a size: the sketch's cells,
    the differences peeled, or the hash's batch."""
    data = []
    for run in runs:
        t = run.table[run.table.family.isin(families) & (run.table.unit == "ns/elem")]
        encode = t[(t.group == "riblt.encode") & summary_selection.default_map(t)]
        peel = t[
            (t.group == "riblt.peel")
            & t.parameter.astype(str).str.contains("prefilter=true,batch=true")
            & summary_selection.default_map(t)
        ]
        h2c = t[t.group == "h2c"].assign(map=MAPS[0][1])
        maps = t[(t.layer == "comparison maps") & t.operation.str.startswith("one map")]
        h2c = pd.concat(
            [
                h2c,
                maps.assign(map=[_map(f) for f in maps.function]),
                # a family's only hash: drawn as the family's own mark
                _unbatched_hashes(t, set(h2c.family)).assign(map=None),
            ]
        )
        for (series, _), d, x in [
            (SCALING[0], encode, [_number(p, "m") for p in encode.parameter]),
            (SCALING[1], peel, [_number(p, "d") for p in peel.parameter]),
            (
                SCALING[2],
                h2c,
                [
                    1 if m != "batch" else n if pd.notna(n) else e
                    for n, e, m in zip(h2c.batch_n, h2c.elements, h2c["mode"])
                ],
            ),
        ]:
            d = d.assign(x=x)
            if "map" not in d:
                d = d.assign(map=None)
            if d.duplicated(["family", "map", "x"]).any():
                raise ValueError(f"ambiguous {run.name} {series}: two rows at a size")
            data += [
                record(run, r, series=series, x=float(r.x), map=r.map)
                for r in d.itertuples()
            ]
    return data


def insertion_findings(data):
    """Describe the observed leaders without treating interval overlap as equality."""
    lines = []
    for name in dict.fromkeys(r["run"] for r in data):
        rows = sorted(
            (r for r in data if r["run"] == name and not r["family"].startswith("xor")),
            key=lambda r: r["value"],
        )
        top = rows[: figures.TOP]
        if not top:
            continue
        lead = "; ".join(
            f"{_label(r['family']).replace(chr(10), ' ')} {r['value']:.1f} ns ({r['lo']:.1f}–{r['hi']:.1f})"
            for r in top
        )
        details = []
        if len(top) >= 2:
            a, b = top[:2]
            overlap = a["lo"] <= b["hi"] and b["lo"] <= a["hi"]
            details.append(
                "The first two recorded intervals "
                + ("overlap." if overlap else "are separated.")
            )
        if len(top) == figures.TOP:
            a, b = top[0], top[-1]
            details.append(
                f"The third estimate takes {100 * (b['value'] / a['value'] - 1):.1f}% more time than the first."
            )
            extras = [
                r["family"]
                for r in rows[figures.TOP :]
                if r["lo"] <= b["hi"] and r["hi"] >= b["lo"]
            ]
            if extras:
                details.append(
                    "Also overlapping the third interval: "
                    + ", ".join(_label(f).replace("\n", " ") for f in extras)
                    + "."
                )
        lines.append(f"**{name}**: {lead}. " + " ".join(details))
    return lines


# The hashes of the families outside the h2c suite, each its family's only
# construction, timed one item at a time: function to construction.
UNBATCHED_HASHES = {
    "sha256 (XOR baseline)": "SHA-256 digest, no curve",
    "ristretto255/hash_from_bytes<Sha512>": "upstream map",
    "secp256k1/ellswift decode": "upstream map",
}


def _unbatched_hashes(t, families=None):
    d = t[
        (t.group == "hash_to_curve")
        & (t["mode"] == "per-element")
        & (t.unit == "ns/elem")
        & t.function.isin(list(UNBATCHED_HASHES))
    ]
    return d if families is None else d[~d.family.isin(families)]


def hash_rows(run, mode):
    t = run.table
    d = t[(t.group == "h2c") & (t["mode"] == mode) & (t.unit == "ns/elem")]
    if mode == "batch":
        d = d[d.batch_n == 1024]
    if mode == "per-element":
        d = pd.concat([d, _unbatched_hashes(t, set(d.family))])
    return _unique(d, f"{run.name} hash {mode}")


def _hash_name(r):
    if r.group == "h2c":
        return MAPS[0][1]
    return UNBATCHED_HASHES.get(r.function) or _map(r.function) or "upstream map"


def batching(runs):
    """Each family's hash one at a time and in batches, per construction:
    try-and-increment, and the comparison layer's single maps."""
    data = []
    for run in runs:
        t = run.table
        maps = t[
            (t.layer == rules.EXPERIMENTAL)
            & (t.operation == "one map")
            & (t.unit == "ns/elem")
            & ((t["mode"] != "batch") | (t.elements == 1024))
        ]
        found = {}
        for mode in ("per-element", "batch"):
            rows = pd.concat([hash_rows(run, mode), maps[maps["mode"] == mode]])
            for _, r in rows.iterrows():
                found[r.family, _hash_name(r), mode] = r
        constructions = list(dict.fromkeys(name for _, name, _ in found))
        for family, name in [(f, n) for f in FAMILIES for n in constructions]:
            a = found.get((family, name, "per-element"))
            b = found.get((family, name, "batch"))
            if a is None and b is None:
                continue
            pair = a is not None and b is not None
            data.append(
                {
                    "run": run.name,
                    "family": family,
                    "map": name,
                    "individual": float(a.value_ns) if a is not None else None,
                    "individual_lo": float(a.value_lo_ns) if a is not None else None,
                    "individual_hi": float(a.value_hi_ns) if a is not None else None,
                    "batch": float(b.value_ns) if b is not None else None,
                    "batch_lo": float(b.value_lo_ns) if b is not None else None,
                    "batch_hi": float(b.value_hi_ns) if b is not None else None,
                    "value": float(a.value_ns / b.value_ns) if pair else None,
                    "lo": float(a.value_lo_ns / b.value_hi_ns) if pair else None,
                    "hi": float(a.value_hi_ns / b.value_lo_ns) if pair else None,
                    "saved_percent": float(100 * (1 - b.value_ns / a.value_ns))
                    if pair
                    else None,
                    "saved_lo_percent": float(100 * (1 - b.value_hi_ns / a.value_lo_ns))
                    if pair
                    else None,
                    "saved_hi_percent": float(100 * (1 - b.value_lo_ns / a.value_hi_ns))
                    if pair
                    else None,
                    "individual_id": a.full_id if a is not None else None,
                    "batch_id": b.full_id if b is not None else None,
                }
            )
    return data


def map_choices(runs):
    modes = {"per-element": "one at a time", "batch": "batches of 1024"}
    data = []
    for run in runs:
        rows = []
        for mode in modes:
            for _, r in hash_rows(run, mode).iterrows():
                if r.family in FAMILIES:
                    rows.append((r, _hash_name(r), None))
        t = run.table
        d = t[
            (t.layer == rules.EXPERIMENTAL)
            & t.family.isin(FAMILIES)
            & (t.unit == "ns/elem")
        ]
        for _, r in d.iterrows():
            # These legacy comparison benches process 1024 elements in each batch.
            if r["mode"] == "batch" and r.elements != 1024:
                continue
            name = _map(r.function)
            if r.operation == "two maps summed":
                name += ", two summed"
            rows.append((r, name, r.operation))
        order = [name for _, name, _ in MAPS]
        rows.sort(
            key=lambda x: (
                FAMILIES.index(x[0].family),
                list(modes).index(x[0]["mode"]),
                order.index(x[1].split(",")[0])
                if x[1].split(",")[0] in order
                else len(order),
                x[1],
            )
        )
        for r, name, semantics in rows:
            extra = {} if semantics is None else {"semantics": semantics}
            data.append(
                record(
                    run,
                    r,
                    label=f"{r.family}: {name}",
                    mode=modes[r["mode"]],
                    **extra,
                )
            )
    return data


def components(runs):
    data = []
    families = []
    all_runs = []
    for run in runs:
        t = run.table
        # All summary batch comparisons use 1024-item batches, including legacy fused benches.
        batched = t["mode"] == "batch"
        correct = (t.batch_n == 1024) | (t.batch_n.isna() & (t.elements == 1024))
        fixed = t[~batched | correct].copy()
        fixed["partial_batch"] = False
        p = br.pipeline_cost(fixed)
        if p.empty:
            continue
        # the baselines have no batched hash to choose
        p = p[(p.hashing == "batched") | p.family.astype(str).isin(BASELINES)]
        if p.empty:
            continue
        # Compare complete compatible recipes at the same fixed workload.
        p = p.loc[p.groupby("family", observed=True)["m=1350"].idxmin()].copy()
        d = p.rename(
            columns={
                "m=1350": "value_ns",
                "m=1350 lo": "value_lo_ns",
                "m=1350 hi": "value_hi_ns",
            }
        )
        keep = [*subset(d), *(f for f in BASELINES if f in set(d.family.astype(str)))]
        families += [f for f in dict.fromkeys(keep) if f not in families]
        all_runs.append((run, p))
    for run, p in all_runs:
        t = run.table
        measured = _unique(
            t[
                (t.group == "riblt.encode")
                & (t.unit == "ns/elem")
                & parameter(t, "m", 1350)
                & parameter(t, "n", 3500)
                & summary_selection.default_map(t)
            ],
            f"{run.name} insertion",
        ).set_index("family")
        for _, r in p[p.family.isin(families)].iterrows():
            data.append(
                {
                    "run": run.name,
                    "family": str(r.family),
                    "value": float(r["m=1350"]),
                    "lo": float(r["m=1350 lo"]),
                    "hi": float(r["m=1350 hi"]),
                    "hash_prepare": float(
                        r.hash_ns + (r.prepare_ns if pd.notna(r.prepare_ns) else 0)
                    ),
                    "adds": float(br.mapping_degree(1350) * r.add_ns),
                    "measured": float(measured.loc[r.family].value_ns)
                    if r.family in measured.index
                    else None,
                    "measured_id": measured.loc[r.family].full_id
                    if r.family in measured.index
                    else None,
                    "recipe": r.recipe,
                    "hashing": r.hashing,
                    "hash_id": r["hash bench"],
                    "prepare_id": r["prepare bench"],
                    "add_id": r["add bench"],
                }
            )
    families = sorted(families, key=lambda f: f != figures.BASELINE)
    return data, families


def _color(family, light=False):
    base = rules.FAMILIES.get(family, (family, None))[0]
    return figures.color(base, light)


def _label(text):
    text = re.sub(rf"^{figures.BASELINE}(?=$|[:,])", r"\g<0> (baseline)", text)
    text = (
        text.replace("gf2_", "binary ").replace("-lambda", " λ").replace("-gls", " GLS")
    )
    return "\n".join(textwrap.wrap(text, 29, break_long_words=False))


def _overlaps(fig):
    """Check the text boxes that readers see, not the hidden shared tick labels."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = list(fig.texts)
    for ax in fig.axes:
        if ax.get_legend() is not None:
            texts += ax.get_legend().get_texts()
        low, high = sorted(ax.get_xlim())
        ticks = [t for t in ax.get_xticklabels() if low <= t.get_position()[0] <= high]
        texts += [
            ax.title,
            ax.xaxis.label,
            ax.yaxis.label,
            *ax.texts,
            *ticks,
            *ax.get_yticklabels(),
        ]
    boxes = [
        (t.get_text(), t.get_window_extent(renderer))
        for t in texts
        if t.get_visible() and t.get_text().strip()
    ]
    return [
        (a, c)
        for i, (a, b) in enumerate(boxes)
        for c, d in boxes[i + 1 :]
        if b.overlaps(d)
    ]


def _save(fig, out, name):
    overlaps = _overlaps(fig)
    if overlaps:
        warnings.warn(f"{name}: overlapping figure text: {overlaps[:3]}", stacklevel=2)
    path = out / f"{name}.svg"
    figures.save(
        fig,
        path.with_suffix(""),
        ("svg",),
        ALT[name],
    )
    figures.dark_twin(path)
    return path.name


def _panels(names, labels, xlabel, ratio=False):
    figures.theme()
    # Keep the complete width at the page width, including family labels.
    height = 1.4 + 0.45 * sum(max(1, len(_label(s).splitlines())) for s in labels)
    fig, axes = plt.subplots(
        1,
        len(names),
        figsize=(figures.FIG_W, max(2.5, height)),
        squeeze=False,
        layout="constrained",
    )
    fig.get_layout_engine().set(h_pad=0.15, hspace=0.06, wspace=0.08)
    for i, (ax, name) in enumerate(zip(axes.flat, names)):
        ax.set_title("\n".join(textwrap.wrap(name, 18)), fontsize=figures.FONT, pad=12)
        ax.set_yticks(
            range(len(labels)),
            [_label(s) for s in labels] if i == 0 else [""] * len(labels),
        )
        ax.set_ylim(len(labels) - 0.5, -0.5)
        ax.set_xscale("log")
        ax.grid(axis="x", alpha=0.18)
        ax.set_axisbelow(True)
        ax.tick_params(labelsize=figures.SMALL)
        ax.xaxis.set_major_locator(figures.Decades(base=10, numticks=3))
        ax.xaxis.set_minor_formatter(NullFormatter())
        if not ratio:
            ax.xaxis.set_major_formatter(FuncFormatter(figures.fmt_tick))
        ax.set_xlabel(xlabel, fontsize=figures.SMALL)
        for side in ["top", "right"]:
            ax.spines[side].set_visible(False)
    return fig, list(axes.flat)


def _box(ax, y, lo, hi, value, family, filled=True, height=0.22):
    """A box from the interval's ends with a line at the estimate: a marker
    at the estimate would cover the ends of most of these intervals. Filled
    boxes take the family's light shade, which marks batching elsewhere."""
    color = _color(family)
    ax.add_patch(
        Rectangle(
            (lo, y - height / 2),
            hi - lo,
            height,
            facecolor=_color(family, light=True) if filled else "none",
            edgecolor=color,
            lw=1,
            zorder=2,
        )
    )
    ax.vlines(value, y - height / 2, y + height / 2, color=color, lw=1.5, zorder=3)


def _box_handle(filled):
    return Patch(facecolor=figures.GRID if filled else "none", edgecolor=figures.TEXT_2)


def _guide(ax, y, a, b):
    """A line from an item's estimate one at a time to its estimate in a
    batch, which pairs the two boxes across the gap batching opens."""
    ax.plot([a, b], [y, y], color=figures.TEXT_2, alpha=0.5, lw=0.8)


def _baseline_handle():
    return Line2D([], [], color=figures.TEXT_2, ls="--", lw=0.8)


MODES = {
    "one at a time": {"filled": False, "dy": -0.15},
    "batches of 1024": {"dy": 0.15},
}


def dots(
    data, out, name, families=None, xlabel="time per item · log scale", ratio=False
):
    names = list(dict.fromkeys(r["run"] for r in data))
    labels = list(families or [])
    labels += [
        label
        for label in dict.fromkeys(r.get("label", r["family"]) for r in data)
        if label not in labels and (not families or label == figures.BASELINE)
    ]
    labels.sort(key=lambda label: not label.startswith(figures.BASELINE))
    if not names or not labels:
        return None
    fig, axes = _panels(names, labels, xlabel, ratio)
    vals = [v for r in data for v in [r["lo"], r["hi"], r["value"]]]
    baseline = not ratio and any(r["family"] == figures.BASELINE for r in data)
    for ax, run in zip(axes, names):
        if ratio:
            ax.axvline(1, color=figures.TEXT_2, ls="--", lw=0.8)
        rows = [r for r in data if r["run"] == run]
        for r in rows:
            if baseline and r["family"] == figures.BASELINE:
                ax.axvline(r["value"], color=figures.TEXT_2, ls="--", lw=0.8)
            label = r.get("label", r["family"])
            mode = MODES.get(r.get("mode"), {})
            _box(
                ax,
                labels.index(label) + mode.get("dy", 0),
                r["lo"],
                r["hi"],
                r["value"],
                r["family"],
                mode.get("filled", True),
            )
        for label in labels:
            pair = {r.get("mode"): r for r in rows if r.get("label") == label}
            if set(MODES) <= set(pair):
                _guide(ax, labels.index(label), *(pair[m]["value"] for m in MODES))
        ax.set_xlim(min(vals) * 0.7, max(vals) * 1.4)
    modes = [m for m in MODES if any(r.get("mode") == m for r in data)]
    handles = [_box_handle(MODES[m].get("filled", True)) for m in modes]
    texts = list(modes)
    if baseline:
        handles.append(_baseline_handle())
        texts.append(f"{figures.BASELINE} (baseline) estimate, across the families")
    if handles:
        fig.legend(
            handles,
            texts,
            loc="outside lower center",
            ncol=len(handles),
            fontsize=figures.SMALL,
            frameon=False,
        )
    return _save(fig, out, name)


def paired(data, out):
    names = list(dict.fromkeys(r["run"] for r in data))
    keys = list(dict.fromkeys((r["family"], r["map"]) for r in data))
    if not names or not keys:
        return None
    alone = {
        k
        for k in keys
        if all(r["batch"] is None for r in data if (r["family"], r["map"]) == k)
    }
    fig, axes = _panels(
        names,
        [
            f"{f}: {m}" + (", one at a time only" if (f, m) in alone else "")
            for f, m in keys
        ],
        "hash ns/item · log scale",
    )
    vals = [
        r[k]
        for r in data
        for k in ["individual_lo", "individual_hi", "batch_lo", "batch_hi"]
        if r[k] is not None
    ]
    for ax, name in zip(axes, names):
        for r in [r for r in data if r["run"] == name]:
            y = keys.index((r["family"], r["map"]))
            if r["individual"] is not None and r["batch"] is not None:
                _guide(ax, y, r["individual"], r["batch"])
            for key, offset in [("individual", -0.15), ("batch", 0.15)]:
                if r[key] is None:
                    continue
                _box(
                    ax,
                    y + offset,
                    r[key + "_lo"],
                    r[key + "_hi"],
                    r[key],
                    r["family"],
                    key == "batch",
                )
        ax.set_xlim(min(vals) * 0.7, max(vals) * 1.4)
    fig.legend(
        [_box_handle(False), _box_handle(True)],
        ["one at a time", "batches of 1024"],
        loc="outside lower center",
        ncol=2,
        fontsize=figures.SMALL,
        frameon=False,
    )
    return _save(fig, out, "batching")


def stacked_log(data):
    """Whether the modeled totals span more than figures.LOG_RATIO."""
    return max(r["hi"] for r in data) > figures.LOG_RATIO * min(r["lo"] for r in data)


def stacks(data, families, out):
    if not data:
        return None
    names = list(dict.fromkeys(r["run"] for r in data))
    unbatched = {r["family"] for r in data if r.get("hashing") != "batched"}
    log = stacked_log(data)
    fig, axes = _panels(
        names,
        [f"{f}, one at a time" if f in unbatched else f for f in families],
        f"modeled ns/item · {'log' if log else 'linear'} scale",
    )
    for ax, name in zip(axes, names):
        if not log:
            ax.set_xscale("linear")
            ax.xaxis.set_major_locator(plt.MaxNLocator(3))
        for r in [r for r in data if r["run"] == name]:
            y = families.index(r["family"])
            ax.barh(
                y,
                r["hash_prepare"],
                height=0.6,
                color=_color(r["family"]),
                edgecolor=figures.TEXT_2,
                lw=0.6,
            )
            ax.barh(
                y,
                r["adds"],
                left=r["hash_prepare"],
                height=0.6,
                color=_color(r["family"], light=True),
                edgecolor=figures.TEXT_2,
                lw=0.6,
            )
            figures.whiskers(ax, y, r["lo"], r["hi"])
            if r.get("measured") is not None:
                ax.plot(
                    r["measured"],
                    y,
                    "|",
                    color=figures.TEXT,
                    markersize=12,
                    markeredgewidth=1.5,
                )
        top = max(max(r["hi"], r.get("measured") or 0) for r in data)
        if log:
            ax.set_xlim(min(r["hash_prepare"] for r in data) / 1.6, top * 1.4)
        else:
            ax.set_xlim(0, top * 1.12)
    # the shades in the neutral hue: each family's bars take its own
    shade = {"edgecolor": figures.TEXT_2, "lw": 0.6}
    handles = [
        Patch(facecolor=figures.color(rules.OTHER), **shade),
        Patch(facecolor=figures.color(rules.OTHER, light=True), **shade),
        Line2D([], [], color=figures.TEXT, lw=1, marker="|", markersize=5),
    ]
    labels = [
        "hash and prepare",
        "\N{MATHEMATICAL ITALIC SMALL K} additions",
        "propagated range",
    ]
    if any(r.get("measured") is not None for r in data):
        handles.append(
            Line2D(
                [],
                [],
                ls="",
                marker="|",
                markersize=12,
                markeredgewidth=1.5,
                color=figures.TEXT,
            )
        )
        labels.append("measured insertion")
    fig.legend(
        handles,
        labels,
        loc="outside lower center",
        ncol=len(handles),
        fontsize=figures.SMALL,
        frameon=False,
    )
    return _save(fig, out, "components")


def _reference_label(family):
    what = "no curve to select" if family == figures.BASELINE else "a fixed group"
    return _label(f"{family}: insertion alone, {what}")


def _stage_styles(refs=False):
    """The line styles the selection figures share, in the neutral hue."""
    handles = [
        Line2D([], [], color=figures.TEXT_2, ls="-"),
        Line2D([], [], color=figures.TEXT_2, ls="--"),
        (
            Patch(facecolor=figures.TEXT_2, alpha=0.16, lw=0),
            Line2D([], [], color=figures.TEXT_2),
        ),
    ]
    labels = [
        "solid: search + certificate",
        "dashed: partial stage",
        "band: fastest seed (thick) to slowest (thin)",
    ]
    if refs:
        handles.append(Line2D([], [], color=figures.TEXT_2, ls=":", lw=1.4))
        labels.append("dotted: baseline, no selection")
    return handles, labels


def _selection_curves(data, refs, out, name, cost, ylabel, guides=()):
    """A panel per run of cost(record, items, x) against the items x
    inserted per salt, items the seeds' fastest and slowest stage in units
    of insertion; refs are flat."""
    names = list(dict.fromkeys(r["run"] for r in [*data, *refs]))
    if not data or not names:
        return None
    figures.theme()
    fig, axes = plt.subplots(
        len(names),
        1,
        # the references' legend entries need the height
        figsize=(figures.FIG_W, (4 if refs else 3.3) * len(names)),
        squeeze=False,
        layout="constrained",
    )
    fig.get_layout_engine().set(h_pad=0.15, hspace=0.06, wspace=0.08)
    xs = [10 ** (2 + i / 12) for i in range(61)]
    for ax, run in zip(axes.flat, names):
        handles, labels = [], []
        for r in [r for r in data if r["run"] == run]:
            low = [cost(r, r["lo_items"], x) for x in xs]
            high = [cost(r, r["hi_items"], x) for x in xs]
            seeds = f"{r['seeds']} {'seeds' if r['seed_ids'] else 'fixture'}"
            color = _color(r["family"])
            ls = "-" if r.get("complete", False) else "--"
            (line,) = ax.plot(xs, low, color=color, ls=ls)
            if r["hi_items"] != r["lo_items"]:
                band = ax.fill_between(xs, low, high, color=color, alpha=0.16, lw=0)
                ax.plot(xs, high, color=color, lw=0.7, ls=ls)
                handles.append((band, line))
                seeds += ", band from fastest to slowest"
            else:
                handles.append(line)
            labels.append(_label(f"{r['family']}: {r['stage']}; {seeds}"))
        for r in [r for r in refs if r["run"] == run]:
            handles.append(
                ax.axhline(
                    r["insert_ns_per_item"], color=_color(r["family"]), ls=":", lw=1.4
                )
            )
            labels.append(_reference_label(r["family"]))
        for y, text in guides:
            ax.axhline(y, color=figures.TEXT_2, ls=":", lw=0.8)
            ax.text(
                xs[-1],
                y,
                text,
                ha="right",
                va="bottom",
                fontsize=figures.SMALL,
                color=figures.TEXT_2,
            )
        ax.set_xscale("log")
        if refs:
            figures.axis(ax, True, "y")
            low = min(r["insert_ns_per_item"] for r in [*data, *refs])
            ax.set_ylim(bottom=low / 2)
        else:
            ax.set_yscale("log")
        ax.set_title(run, fontsize=figures.FONT)
        ax.set_xlabel("items inserted per salt, N · log scale", fontsize=figures.SMALL)
        ax.set_ylabel(ylabel, fontsize=figures.SMALL)
        ax.tick_params(labelsize=figures.SMALL)
        ax.legend(
            handles,
            labels,
            loc="upper left",
            bbox_to_anchor=(1.01, 1),
            fontsize=figures.SMALL,
            frameon=False,
        )
    handles, labels = _stage_styles(bool(refs))
    fig.legend(
        handles,
        labels,
        loc="outside lower center",
        ncol=len(handles),
        fontsize=figures.SMALL,
        frameon=False,
    )
    return _save(fig, out, name)


def amortization(data, out):
    return _selection_curves(
        data,
        [],
        out,
        "selection-amortization",
        lambda r, items, x: 100 * items / x,
        "selection / insertion work (%) · log scale",
        [(100, "equal work"), (10, "one tenth")],
    )


def per_item(data, refs, out):
    """Insertion and selection amortised over N items, in ns per item."""
    return _selection_curves(
        data,
        refs,
        out,
        "selection-per-item",
        lambda r, items, x: r["insert_ns_per_item"] * (1 + items / x),
        "\N{MATHEMATICAL ITALIC SMALL C} + \N{MATHEMATICAL ITALIC CAPITAL S}"
        "/\N{MATHEMATICAL ITALIC CAPITAL N}, time per item · log scale",
    )


def selection_times(data, refs, out):
    """Each stage's time per salt: its seeds' estimates and their range."""
    if not data:
        return None
    names = list(dict.fromkeys(r["run"] for r in [*data, *refs]))
    stages = list(dict.fromkeys(f"{r['family']}: {r['stage']}" for r in data))
    shown = [f for f in BASELINES if any(r["family"] == f for r in refs)]
    labels = [*shown, *stages]
    fig, axes = _panels(names, labels, "selection time per salt · log scale")
    values = [v for r in data for v in r["stage_ns"]]
    height = 0.3
    for ax, name in zip(axes, names):
        for r in [r for r in data if r["run"] == name]:
            y = labels.index(f"{r['family']}: {r['stage']}")
            color = _color(r["family"])
            ls = "-" if r.get("complete", False) else "--"
            lo, hi = min(r["stage_ns"]), max(r["stage_ns"])
            if hi > lo:
                ax.add_patch(
                    Rectangle(
                        (lo, y - height / 2),
                        hi - lo,
                        height,
                        facecolor=color,
                        alpha=0.16,
                        lw=0,
                        zorder=2,
                    )
                )
                ax.add_patch(
                    Rectangle(
                        (lo, y - height / 2),
                        hi - lo,
                        height,
                        facecolor="none",
                        edgecolor=color,
                        ls=ls,
                        lw=1,
                        zorder=2,
                    )
                )
            ax.vlines(
                r["stage_ns"],
                y - height / 2,
                y + height / 2,
                color=color,
                ls=ls,
                lw=1.5,
                zorder=3,
            )
        for r in [r for r in refs if r["run"] == name and r["family"] in shown]:
            ax.text(
                0.02,
                labels.index(r["family"]),
                "none: no curve to select"
                if r["family"] == figures.BASELINE
                else "none: a fixed group",
                transform=ax.get_yaxis_transform(),
                va="center",
                fontsize=figures.SMALL,
                color=figures.TEXT_2,
            )
        ax.set_xlim(min(values) * 0.5, max(values) * 2)
    handles = [
        Line2D([], [], color=figures.TEXT_2, ls="none", marker="|", ms=10, mew=1.5),
        Patch(facecolor=figures.TEXT_2, alpha=0.16, edgecolor="none"),
        Line2D([], [], color=figures.TEXT_2, ls="--"),
    ]
    fig.legend(
        handles,
        [
            "a seed's or fixture's estimate",
            "fastest to slowest seed",
            "dashed: partial stage",
        ],
        loc="outside lower center",
        ncol=len(handles),
        fontsize=figures.SMALL,
        frameon=False,
    )
    return _save(fig, out, "selection-time")


RUN_STYLES = ["-", "--", ":"]
# Tableau 10, the 2016 palette
CATEGORICAL = [
    "#4e79a7",
    "#f28e2b",
    "#e15759",
    "#76b7b2",
    "#59a14f",
    "#edc948",
    "#b07aa1",
    "#ff9da7",
    "#9c755f",
    "#bab0ac",
]
# Tableau 20's light tints of the hues above, grey excepted
TINTS = [
    "#a0cbe8",
    "#ffbe7d",
    "#ff9d9a",
    "#86bcb6",
    "#8cd17d",
    "#f1ce63",
    "#d4a6c8",
    "#fabfd2",
    "#d7b5a6",
]


def _colors(families):
    """A colour per family: the first family of a base keeps its base's,
    as in the other figures, and the rest take the categorical colours no
    base here holds, then their tints. Lines of one hue would not be told
    apart where they cluster; the selection of several runs can exceed
    both palettes, and colours then repeat rather than fail."""
    colors = {}
    for f in families:
        c = figures.color(rules.FAMILIES.get(f, (f, None))[0])
        if c not in colors.values():
            colors[f] = c
    spare = itertools.chain(
        (c for c in CATEGORICAL + TINTS if c not in colors.values()),
        itertools.cycle(CATEGORICAL),
    )
    return {f: colors.get(f) or next(spare) for f in families}


def curves(data, families, out):
    names = list(dict.fromkeys(r["run"] for r in data))
    panels = [s for s in SCALING if any(r["series"] == s[0] for r in data)]
    if not panels:
        return None
    figures.theme()
    fig = plt.figure(figsize=(figures.FIG_W, 4.4), layout="constrained")
    axes = fig.subplots(1, len(panels), squeeze=False).flat
    families = [f for f in families if any(r["family"] == f for r in data)]
    colors = _colors(families)
    maps = [m for m in MAPS if any(r.get("map") == m[1] for r in data)]
    markers = {name: marker for _, name, marker in maps}
    for ax, (series, xlabel) in zip(axes, panels):
        for run, ls in zip(names, RUN_STYLES):
            for f, m in [(f, m) for f in families for m in [None, *markers]]:
                g = sorted(
                    (
                        r
                        for r in data
                        if (r["series"], r["run"], r["family"], r.get("map"))
                        == (series, run, f, m)
                    ),
                    key=lambda r: r["x"],
                )
                if not g:
                    continue
                c = colors[f]
                x = [r["x"] for r in g]
                ax.plot(
                    x,
                    [r["value"] for r in g],
                    ls=ls,
                    marker=markers.get(m, "o"),
                    ms=2.5 if m in (None, MAPS[0][1]) else 4,
                    mfc="none" if m not in (None, MAPS[0][1]) else c,
                    lw=1.3,
                    color=c,
                )
                ax.fill_between(
                    x,
                    [r["lo"] for r in g],
                    [r["hi"] for r in g],
                    color=c,
                    alpha=0.18,
                    lw=0,
                )
        ax.set_title(
            "\n".join(textwrap.wrap(series, 30, break_on_hyphens=False)),
            loc="left",
            fontsize=figures.FONT,
        )
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
        ax.xaxis.set_minor_formatter(NullFormatter())
        figures.axis(ax, True, "y")
        ax.grid(axis="x", which="major", color=figures.GRID, linewidth=0.8)
        ax.set_xlabel(f"{xlabel}, log scale")
    handles = [
        Line2D(
            [],
            [],
            color=colors[f],
            marker="o",
            ms=2.5,
        )
        for f in families
    ]
    labels = [f"{f} (insecure baseline)" if f == "xor" else f for f in families]
    if len(names) > 1:
        handles += [
            Line2D([], [], color=figures.TEXT_2, ls=ls)
            for _, ls in zip(names, RUN_STYLES)
        ]
        labels += names[: len(RUN_STYLES)]
    if len(maps) > 1:
        handles += [
            Line2D(
                [],
                [],
                color=figures.TEXT_2,
                ls="none",
                marker=marker,
                ms=2.5 if name == MAPS[0][1] else 4,
                mfc=figures.TEXT_2 if name == MAPS[0][1] else "none",
            )
            for _, name, marker in maps
        ]
        labels += [name for _, name, _ in maps]
    fig.legend(handles, labels, loc="outside lower center", ncol=4)
    return _save(fig, out, "scaling")


ALT = {
    "scaling": "Line graphs on log axes of measured time per item or difference against the number of cells, the number of differences and the hash batch size, a line per family, map and run, intervals as bands.",
    "insertion": "Box plot of measured RIBLT encoding time per item by family and run, on a log scale, with XOR as a dashed reference.",
    "batching": "Paired hollow and filled boxes compare individual and batched hashing times for each family and run on a log scale.",
    "native-generic": "Box plot of generic-to-native time ratios for matched operations, with a dashed parity line.",
    "map-choice": "Box plot of measured hash-to-curve time per item for each family and map construction, hollow boxes one item at a time and filled boxes in batches of 1024, on a log scale.",
    "components": "Stacked bars split modeled insertion work into hashing plus preparation and lighter addition work; whiskers show propagated ranges and ticks mark measured insertion.",
    "selection-amortization": "Log-log curves show selection-stage work relative to measured insertion work as item volume grows, with observed seed ranges and dashed partial stages.",
    "selection-time": "Box plot on a log scale of each curve selection stage's time per salt by family and run, a tick at each seed's estimate and a band across the seeds, with rows stating that the baselines select no curve.",
    "selection-per-item": "Log-log curves show insertion plus amortised selection time per item as item volume grows, with observed seed ranges, dashed partial stages, and the baselines' insertion as flat dotted lines.",
}


def _figure(svg, caption):
    if svg is None:
        return ""
    alt = html.escape(ALT[svg.removesuffix(".svg")], quote=True)
    return (
        f'<img class="fig-light" src="{svg}" alt="{alt}">'
        f'<img class="fig-dark" src="{svg.removesuffix(".svg")}-dark.svg" alt="{alt}">\n\n{caption}\n'
    )


def chapter(runs, out):
    if not runs:
        return None
    out.mkdir(parents=True, exist_ok=True)
    text = [
        "# Summary across runs",
        (
            "Each run's measurements stay separate; boxes span recorded 95% intervals, with a line at the estimate, "
            "read as in [Reading a run report](../methodology.md#reading-a-run-report). "
            "[The displayed measurements](evidence.csv) keep their benchmark identifiers."
        ),
    ]
    text += [
        f"- **{r.name}**: {br.ran({k: v for k, v in r.meta.items() if k != 'name'})}"
        for r in runs
    ]
    evidence = []
    coverage = {r.name: [] for r in runs}
    views_rendered = []

    def section(title, data, svg, caption):
        views_rendered.append(title)
        text.extend(
            [
                f"## {title}",
                _figure(svg, caption)
                if svg
                else (
                    "This comparison needs a native and a generic run of the same commit on the same recorded machine and compiler."
                    if title == "Native and generic builds"
                    else "Not measured: this view has no eligible observations."
                ),
            ]
        )
        have = {r["run"] for r in data}
        have.update(r[k] for r in data for k in ["native_run", "generic_run"] if k in r)
        missing = [r.name for r in runs if r.name not in have]
        for name, views in coverage.items():
            if name not in missing:
                views.append(title)
        evidence.extend({"view": title, **r} for r in data)

    ins, families = insertion(runs)
    section(
        "Insertion",
        ins,
        dots(ins, out, "insertion", families),
        "RIBLT encoding per item: $n = 3500$ items into a sketch of $m = 1350$ cells, "
        "the $1.35d$ coded symbols that decode $d = 1000$ differences, "
        "under each family's cheapest measured hash; "
        "each run's three fastest curves, those whose intervals overlap the third, and the representative families, "
        "ristretto255 among them. "
        "The first row is the insecure XOR baseline, SHA-256 digests XORed into 32-byte cells; "
        "the dashed line carries its estimate across the families.",
    )
    text.extend(insertion_findings(ins))
    if ins:
        text.append(
            "Coverage: "
            + "; ".join(
                f"{name}: {next(r['selected'] for r in ins if r['run'] == name)} displayed curves of {next(r['available'] for r in ins if r['run'] == name)} measured rows (including baselines)"
                for name in dict.fromkeys(r["run"] for r in ins)
            )
            + "."
        )
    scale = scaling(runs, families + ["xor"])
    section(
        "Scaling",
        scale,
        curves(scale, families + ["xor"], out),
        "Encoding per item against the sketch's cells $m$, $n = 3500$ items; "
        "peeling per difference against the set difference $d$, in a sketch of $m = 2d$ cells "
        "with the prefilter and batched hashes; hashing per item against the batch size, "
        "try-and-increment filled and the comparison maps open; "
        "the hashes of ristretto255 and of the XOR baseline, a SHA-256 digest, have no batched form "
        "and are a single mark at batch size 1. Bands are recorded intervals"
        + ("; line style tells the runs apart." if len(runs) > 1 else "."),
    )
    bat = batching(runs)
    section(
        "Batching",
        bat,
        paired(bat, out),
        "Hashing per item by construction, one item at a time and in batches of 1024; "
        "a grey line joins a construction's two estimates. "
        "Batching shares the field inversion; the prime-field hashes are bound by a square root per item, which it does not share. "
        "ristretto255's hash and the XOR baseline's SHA-256 digest have no batched form: "
        "each is a single hollow box, labelled one at a time only. "
        "The table gives the time saved, $100 \\times (1 - \\text{batched} / \\text{individual})$.",
    )
    if bat:
        text.append(
            "| Run | Family | Hash | Batch time change |\n|---|---|---|---|\n"
            + "\n".join(
                f"| {r['run']} | {_label(r['family']).replace(chr(10), ' ')} | {r['map']} | "
                + (
                    f"{abs(r['saved_percent']):.1f}% {'saved' if r['saved_percent'] >= 0 else 'more time'}"
                    if r["saved_percent"] is not None
                    else (
                        "individual only; the digest is timed one at a time"
                        if r["family"] == figures.BASELINE
                        and r["individual"] is not None
                        else "individual only; comparison harness has no batched entrypoint"
                        if r["family"] in ("ristretto255", "secp256k1")
                        and r["individual"] is not None
                        else "only one mode measured"
                    )
                )
                + " |"
                for r in bat
            )
        )
    pairs = summary_pairs.compare(runs)
    for r in pairs:
        r["label"] = f"{r['family']}: {r['operation']}"
    section(
        "Native and generic builds",
        pairs,
        dots(
            pairs,
            out,
            "native-generic",
            xlabel="generic / native · log scale",
            ratio=True,
        ),
        "Generic over native time for benchmarks matched on one CPU, source and compiler; "
        "the dashed line is parity, boxes the propagated ranges. "
        "ristretto255's point addition and the XOR baseline's addition of 32-byte digests, "
        "both at throughput, are matched as references where both runs timed them.",
    )
    for name in dict.fromkeys(r["generic_run"] for r in pairs):
        models = sorted(
            {
                r["family"]
                for r in pairs
                if r["generic_run"] == name and r["operation"].startswith("point ")
            }
        )
        if models:
            text.append(
                f"Matched group models in {name}: "
                + ", ".join(_label(f).replace("\n", " ") for f in models)
                + ". Other models need corresponding generic measurements."
            )
    maps = map_choices(runs)
    section(
        "Hash to curve",
        maps,
        dots(maps, out, "map-choice"),
        "Hashing per item by construction for the representative families, one item at a time "
        "and in batches of 1024; a grey line joins a construction's two estimates. "
        "Constructions differ in the representation they return, "
        "as the benchmark identifiers in the data state. "
        "The XOR baseline's row is its SHA-256 digest, which maps to no curve, "
        "and the dashed line carries its estimate across the families; "
        "ristretto255's and secp256k1's hashes are their libraries' maps, timed one at a time only.",
    )
    comp, families = components(runs)
    absent = [
        name
        for name in dict.fromkeys(r["run"] for r in comp)
        if not any(r["run"] == name and r["family"] == figures.BASELINE for r in comp)
    ]
    section(
        "Insertion cost model",
        comp,
        stacks(comp, families, out),
        "Modeled insertion per item under each family's cheapest compatible recipe, "
        "a fused hash to the addend counting as hash and prepare, and additions at throughput. "
        "A sketch of $m = 1350$ cells, the $1.35d$ coded symbols that decode $d = 1000$ differences, "
        f"maps an item to $k(m) = {br.mapping_degree(1350):.2f}$ cells on average "
        "([Workload](../workload.md#cost-in-repeated-reconciliation)). "
        "The measured insertion includes what the model omits. "
        "ristretto255 and the XOR baseline have no batched hash: they hash one item at a time "
        "and add the hash output as is, with nothing to prepare."
        + (
            " The XOR baseline is absent from "
            + ", ".join(absent)
            + ": its addition was not timed at throughput, with eight accumulators, there."
            if absent
            else ""
        )
        + (
            " The time axis is logarithmic, as the totals span more than a factor of ten: "
            "a segment's boundaries read as cumulative times, its length not as its share."
            if comp and stacked_log(comp)
            else ""
        ),
    )
    sel = summary_selection.selection(runs)
    refs = [r for r in summary_selection.references(runs) if sel]
    seeds = (
        "The band spans the seeds, the fastest thick and the slowest thin; "
        "dashed curves are partial stages, lower bounds on $S$."
    )
    section(
        "Curve selection",
        sel,
        amortization(sel, out),
        "Selection work per salt $S$ as a percentage of the insertion work $Nc$ it is amortised over, "
        "$100 S / (N c)$, against the items $N$ inserted per salt, "
        "with $c$ the measured insertion per item "
        "($n = 3500$ items into a sketch of $m = 1350$ cells, as in Insertion). "
        f"{seeds} "
        "Normalised by each family's own $c$, the figure reads as the volume that amortises selection, "
        "but it credits a family that is slow to insert with a smaller ratio for the same $S$; "
        "the next two figures give $S$ and the cost per item in absolute time. "
        "The baselines, XOR and ristretto255, select no curve: their ratio is zero, "
        "which a log axis cannot show, and they appear in the absolute figures.",
    )
    section(
        "Curve selection time",
        [*sel, *refs],
        selection_times(sel, refs, out),
        "The time $S$ per salt of each selection stage, on a log scale: "
        "a tick at each seed's estimate, or at the single candidate fixture of a point count, "
        "and a band from the fastest seed to the slowest; dashed marks are partial stages, lower bounds on $S$. "
        "Unlike the ratio above, these times compare across families whatever their insertion cost. "
        "XOR has no curve and ristretto255 is a fixed group: neither selects a curve, as their rows state.",
    )
    section(
        "Cost per item with selection",
        [*sel, *refs],
        per_item(sel, refs, out),
        "Insertion and amortised selection per item, $c + S / N$, against the items $N$ inserted per salt, "
        "with $c$ the measured insertion per item "
        "($n = 3500$ items into a sketch of $m = 1350$ cells) and $S$ the selection time per salt above. "
        "As $N$ grows each curve falls to its family's $c$. "
        f"{seeds} "
        "The dotted flat lines are the insertion of the baselines, XOR and ristretto255, "
        "which select no curve, so that $c$ is their whole cost per item.",
    )
    gaps = [
        f"{name} appears in {', '.join(views) if views else 'no eligible summary views'}."
        for name, views in coverage.items()
        if len(views) < len(views_rendered)
    ]
    if gaps:
        text.insert(2, "Coverage: " + " ".join(gaps))
    pd.DataFrame(evidence).to_csv(out / "evidence.csv", index=False)
    return f"<style>{br.BOOK_STYLE}</style>\n\n" + "\n\n".join(text) + "\n"
