"""The report's figures, drawn at the width at which the book shows them and in
a palette whose neutrals a dark twin swaps.

A comparison across families shows the TOP fastest, the REPRESENTATIVES and
the baseline. Every other family is still in the grid or table beside the
figure, which keeps every measurement, and the figure says how many it shows.
A parameter swept by a suite (cells m, differences d, batch size n) is an
axis, not a row per value.
"""

import math
import re
from pathlib import Path
from typing import NamedTuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import rules
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter, PercentFormatter

TEXT = "#0b0b0b"
TEXT_2 = "#52514e"
GRID = "#e4e3df"
SURFACE = "#fcfcfb"
# the neutrals of a figure's twin for the book's dark themes; the categorical
# colours hold up on either surface and are kept
DARK = {TEXT: "#e8e6e3", TEXT_2: "#a8a6a1", GRID: "#3a3d41", SURFACE: "#1c1e21"}

# The book's text column is 900px wide, the width of an SVG this many inches
# wide at CSS's 96px to the inch, so type is shown at its nominal size.
FIG_W = 9.0
FONT = 9
SMALL = 8
TITLE = 10
HEAD = 12
# a figure whose times span more than this ratio is on a log scale throughout
LOG_RATIO = 10

TOP = 3
# the README's candidates, then the references riblt-ecmh and Bitcoin use
REPRESENTATIVES = (
    "gf2_109",
    "gf2_127",
    "gf2_122-gls",
    "edwards127",
    "weier127",
    "ristretto255",
    "secp256k1",
)
BASELINE = "xor"


class Figure(NamedTuple):
    """A figure's files and the sentence that describes it, its alt text
    and caption."""

    names: list[str]
    alt: str


def theme():
    """Set every neutral to one of the palette's, which the dark twin swaps."""
    plt.rcParams.update(
        {
            "svg.hashsalt": "bench-report",
            "svg.fonttype": "none",
            "font.size": FONT,
            "axes.edgecolor": TEXT_2,
            "axes.facecolor": SURFACE,
            "axes.labelcolor": TEXT,
            "axes.labelsize": SMALL,
            "axes.titlesize": TITLE,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "figure.facecolor": SURFACE,
            "legend.frameon": False,
            "legend.fontsize": SMALL,
            "text.color": TEXT,
            "xtick.color": TEXT_2,
            "xtick.labelsize": SMALL,
            "ytick.color": TEXT_2,
            "ytick.labelsize": FONT,
        }
    )


def save(fig, stem: Path, formats, alt: str, footer: str = "") -> Figure:
    """Write fig as stem.<format> for each format, footer under it, and
    close it."""
    if footer:
        # above the title, clear of the legends the layout puts below
        fig.text(
            1, 1, footer, fontsize=SMALL - 1, color=TEXT_2, ha="right", va="bottom"
        )
    names = []
    for fmt in formats:
        name = f"{stem.name}.{fmt}"
        meta = {"Date": None} if fmt == "svg" else {}
        fig.savefig(
            stem.parent / name,
            dpi=150,
            bbox_inches="tight",
            pad_inches=0.1,
            facecolor=SURFACE,
            metadata=meta,
        )
        names.append(name)
    plt.close(fig)
    return Figure(names, alt)


def dark(name: str) -> str:
    """The file name of a figure's twin for the book's dark themes."""
    return name.removesuffix(".svg") + "-dark.svg"


def dark_twin(svg: Path) -> Path:
    """Write svg's twin beside it: the figure with its neutrals swapped for
    DARK's."""
    text = svg.read_text()
    for light, shade in DARK.items():
        text = text.replace(light, shade)
    twin = svg.with_name(dark(svg.name))
    twin.write_text(text)
    return twin


def svg(fig: Figure) -> str:
    """The file of fig that the book and the reports show."""
    return next((n for n in fig.names if n.endswith(".svg")), fig.names[0])


def fmt_time(ns: float) -> str:
    for scale, unit in ((1e9, "s"), (1e6, "ms"), (1e3, "µs"), (1, "ns")):
        if abs(ns) >= scale or unit == "ns":
            v = ns / scale
            return (
                f"{v:.2f} {unit}"
                if abs(v) < 10
                else f"{v:.1f} {unit}"
                if abs(v) < 100
                else f"{v:.0f} {unit}"
            )
    raise AssertionError


def fmt_tick(ns: float, _=None) -> str:
    for scale, unit in ((1e9, "s"), (1e6, "ms"), (1e3, "µs"), (1, "ns")):
        if ns >= scale or unit == "ns":
            return f"{ns / scale:g} {unit}"
    raise AssertionError


def batch_change(alone: float, batched: float) -> str:
    """What batching saves on the time alone, or costs, as a percentage."""
    if batched <= alone:
        return f"{1 - batched / alone:.1%} saved"
    return f"{batched / alone - 1:.1%} more"


def wrap_label(s: str, width: int = 50) -> str:
    """s in lines of at most width where it can be: broken after a space, a
    slash or a comma, never inside an identifier or a parameter."""
    lines, line = [], ""
    for part in re.findall(r"[^ /,]*[ /,]?", s):
        if line and len(line) + len(part.rstrip()) > width:
            lines.append(line.rstrip())
            line = ""
        line += part
    return "\n".join(lines + [line.rstrip()])


def color(base: str, light: bool = False) -> str:
    shades = rules.LIGHT if light else rules.COLORS
    return shades.get(base, shades[rules.OTHER])


def shown(best: pd.DataFrame, top: int = TOP) -> list:
    """Of best, a row per family indexed by it with value_ns, value_lo_ns and
    value_hi_ns, the families a comparison shows, in best's order: the top
    fastest that are not baselines, up to as many more whose intervals reach
    the slowest of those (a tie at the cut is not cut), the representatives
    and the baseline."""
    rivals = best[~best.index.str.startswith(BASELINE)].sort_values("value_ns")
    first = rivals.head(top)
    rest = rivals.iloc[top:]
    tied = rest[rest.value_lo_ns <= first.value_hi_ns.max()].head(top)
    keep = {*first.index, *tied.index, *REPRESENTATIVES, BASELINE}
    return [f for f in best.index if f in keep]


def per_base(best: pd.DataFrame, bases: dict) -> list:
    """Of best, as for shown(), the fastest family of each base family in
    bases, family to base, in the order of rules.FAMILY_ORDER: lines of one
    hue would not be told apart, and the grids list the rest."""
    first = {}
    for f in best.sort_values("value_ns").index:
        first.setdefault(bases[f], f)
    order = {b: i for i, b in enumerate(rules.FAMILY_ORDER)}
    return sorted(first.values(), key=lambda f: order.get(bases[f], len(order)))


def of(n: int, total: int) -> str:
    return f"{n} of {total} shown" if n < total else ""


def whiskers(ax, y, lo, hi):
    """Each interval [lo, hi] at its y, drawn from its own endpoints: a
    bootstrap interval can exclude the slope's point estimate, and is drawn
    as measured all the same."""
    ax.hlines(y, lo, hi, color=TEXT, linewidth=1, zorder=3)
    for x in (lo, hi):
        ax.scatter(x, y, marker="|", s=25, color=TEXT, linewidths=1, zorder=3)


def interval_boxes(ax, y, lo, hi, value, fills, *, height=0.3):
    """Measured confidence bounds as boxes, with the estimate as a line.

    Keep the estimate independent of the box: a bootstrap interval can
    exclude the fitted slope. Never widen a narrow interval to fit a marker.
    """
    for yi, low, high, estimate, fill in zip(y, lo, hi, value, fills):
        ax.add_patch(
            Rectangle(
                (low, yi - height / 2),
                high - low,
                height,
                facecolor=fill,
                edgecolor=TEXT,
                linewidth=0.8,
                zorder=3,
            )
        )
        ax.vlines(
            estimate,
            yi - height / 2,
            yi + height / 2,
            color=TEXT,
            linewidth=1.2,
            zorder=4,
        )


def interval_caption(fig):
    fig.supxlabel(
        "Boxes span 95% confidence intervals of Criterion’s estimates; vertical lines mark the estimates.\n"
        "Hollow boxes: individual operations. Filled boxes: batched. Very narrow intervals may be subpixel.",
        fontsize=SMALL,
        color=TEXT_2,
    )


class Decades(LogLocator):
    """Ticks at the powers of ten, and at 2 and 5 times them on an axis
    that spans too few powers to be read by them alone."""

    def tick_values(self, vmin, vmax):
        few = 0 < vmin < vmax and math.log10(vmax / vmin) < 2.5
        self.set_params(subs=(1.0, 2.0, 5.0) if few else (1.0,))
        return super().tick_values(vmin, vmax)


def axis(ax, log: bool, grid_axis: str = "x"):
    """Time ticks on the value axis, its scale, and the grid behind."""
    which = ax.xaxis if grid_axis == "x" else ax.yaxis
    if log:
        (ax.set_xscale if grid_axis == "x" else ax.set_yscale)("log")
        which.set_major_locator(Decades(base=10))
        which.set_minor_formatter(NullFormatter())
    which.set_major_formatter(FuncFormatter(fmt_tick))
    ax.grid(axis=grid_axis, which="major", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def spans(d: pd.DataFrame) -> bool:
    """Whether d's times, the baseline aside, span more than LOG_RATIO."""
    d = d[d.base != BASELINE] if (d.base != BASELINE).any() else d
    return d.value_hi_ns.max() > LOG_RATIO * d.value_lo_ns.min()


def value_axis(ax, d: pd.DataFrame, log: bool):
    axis(ax, log)
    top = d.value_hi_ns.max()
    if log:
        ax.set_xlim(d.value_lo_ns.min() / 1.6, top * 3.5)
    else:
        ax.set_xlim(0, top * 1.3)


def time_labels(ax, y, d: pd.DataFrame, log: bool):
    top = d.value_hi_ns.max()
    for yi, v, hi in zip(y, d.value_ns, d.value_hi_ns):
        x = max(v, hi) * 1.12 if log else max(v, hi) + 0.02 * top
        ax.text(x, yi, fmt_time(v), va="center", fontsize=SMALL, color=TEXT_2)


def shaded_barh(ax, y, d: pd.DataFrame, light, height):
    """Bars in their family's hue, light where light holds; a light bar is
    edged in the dark shade, which keeps its extent visible on a pale
    surface."""
    ax.barh(
        y,
        d.value_ns,
        height=height,
        color=[color(b, li) for b, li in zip(d.base, light)],
        edgecolor=[color(b) if li else SURFACE for b, li in zip(d.base, light)],
        linewidth=0.8,
    )


# Bars: one row per family or method, an operation per facet.


def pair_bars(
    ax, d: pd.DataFrame, rows: list, variants: list, log: bool, row="curve", slots=None
):
    """A facet's horizontal bars: a row per entry of rows, top down, by
    d[row], and per row a bar per variant of variants that it has, the first
    dark and above, the second light and below; slots rows tall, so that
    facets side by side draw bars alike. A batched variant is labelled
    beyond the axis with what it saves on the first."""
    pair = len(variants) > 1
    height = 0.38 if pair else 0.7
    offset = {v: (0.2 - 0.4 * i if pair else 0) for i, v in enumerate(variants)}
    d = d[d[row].isin(rows) & d.variant.isin(variants)]
    y = [-rows.index(r) + offset[v] for r, v in zip(d[row], d.variant)]
    second = [pair and v == variants[1] for v in d.variant]
    shaded_barh(ax, y, d, second, height)
    interval_boxes(
        ax,
        y,
        d.value_lo_ns,
        d.value_hi_ns,
        d.value_ns,
        [
            color(base, True) if variant == "batched" else SURFACE
            for base, variant in zip(d.base, d.variant)
        ],
        height=0.25 if pair else 0.4,
    )
    value_axis(ax, d, log)
    time_labels(ax, y, d, log)
    first = {
        r: v for r, var, v in zip(d[row], d.variant, d.value_ns) if var == variants[0]
    }
    for yi, r, s, v in zip(y, d[row], second, d.value_ns):
        if s and variants[1] == "batched" and r in first:
            # beyond the axis, not after the time, so that a time and a
            # saving are not run together
            ax.text(
                1.02,
                yi,
                batch_change(first[r], v),
                va="center",
                ha="left",
                fontsize=SMALL,
                style="italic",
                color=TEXT_2,
                transform=ax.get_yaxis_transform(),
            )
    ax.set_yticks([-i for i in range(len(rows))], [wrap_label(r, 40) for r in rows])
    ax.set_ylim(-(slots or len(rows)) + 0.4, 0.6)
    ax.tick_params(axis="y", length=0)


def titled(ax, title: str, note: str = ""):
    """title over the axes, and note in smaller type between them."""
    ax.set_title(title, loc="left", color=TEXT, pad=16 if note else 6)
    if note:
        ax.annotate(
            note,
            (0, 1),
            xycoords="axes fraction",
            xytext=(0, 4),
            textcoords="offset points",
            fontsize=SMALL,
            color=TEXT_2,
            va="bottom",
        )


def variant_key(variants: list) -> str:
    return f"dark {variants[0]}, light {variants[1]}" if len(variants) > 1 else ""


def plot_elementary(
    e: pd.DataFrame, out: Path, formats, footer=""
) -> dict[str, Figure]:
    """Two figures, field operations and the rest: a facet per operation,
    a row per field or curve that the comparison shows."""
    names = {}
    for kind, (stem, title) in {
        "field": ("elementary-field", "Field operations, time per operation"),
        "group": (
            "elementary-group",
            "Group operations and hash to curve, time per operation",
        ),
    }.items():
        d = e[(e.layer == "field") == (kind == "field")]
        if len(d):
            names[kind] = _plot_elementary(d, out / stem, title, formats, footer)
    return names


def _plot_elementary(e: pd.DataFrame, stem: Path, title: str, formats, footer):
    # Every facet on a log scale: the operations differ by orders of
    # magnitude, and a linear facet among log ones reads as a different
    # spread.
    ncol = 2
    plans = []
    for facet in dict.fromkeys(e.facet):
        d = e[e.facet == facet]
        variants = list(dict.fromkeys(d.variant))
        lead = d[d.variant == variants[0]].set_index("curve")
        plans.append((facet, d, variants, shown(lead), len(lead)))
    grid = [plans[i : i + ncol] for i in range(0, len(plans), ncol)]
    heights = [
        0.9 + max(len(p[3]) * (0.5 if len(p[2]) > 1 else 0.3) for p in row)
        for row in grid
    ]
    fig = plt.figure(figsize=(FIG_W, sum(heights) + 0.6), layout="constrained")
    fig.get_layout_engine().set(h_pad=0.15, hspace=0.06, wspace=0.08)
    axes = fig.subplots(len(grid), ncol, height_ratios=heights, squeeze=False)
    for r, row in enumerate(grid):
        for c in range(ncol):
            ax = axes[r][c]
            if c >= len(row):
                ax.set_visible(False)
                continue
            facet, d, variants, rows, total = row[c]
            slots = max(len(p[3]) for p in row)
            pair_bars(ax, d, rows, variants, log=True, slots=slots)
            note = ", ".join(
                s for s in (variant_key(variants), of(len(rows), total)) if s
            )
            titled(ax, facet.removeprefix("field "), note)
    fig.suptitle(title, x=0, ha="left", fontsize=HEAD, color=TEXT)
    interval_caption(fig)
    alt = (
        f"{title}, on log scales: a panel per operation with a bar per field "
        "or curve and its 95% confidence-interval box with an estimate line, two variants of an operation as a "
        "dark and a light bar of the family's colour; each panel shows the "
        "three fastest, the representative families and the XOR baseline, "
        "and the grid above lists the rest."
    )
    return save(fig, stem, formats, alt, footer)


# Layers: the reference figure of each suite layer.


def plot_layer(
    t: pd.DataFrame, layer: str, out: Path, formats, footer=""
) -> Figure | None:
    """The layer's figure, or None where the elementary figures already show
    its operations (fields, group operations)."""
    d = t[t.layer == layer]
    stem = out / slug(layer)
    if layer in ("field", "group ops"):
        return None
    sweep = (t.group == "h2c").any()
    if layer == "hash to curve" and sweep:
        return _hash_sweep(d, stem, formats, footer)
    if layer == rules.EXPERIMENTAL and sweep:
        return _maps(t, stem, formats, footer)
    if layer == "curve generation":
        return _seeds(d, stem, formats, footer)
    return _layer_bars(d, layer, stem, formats, footer)


def slug(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", s.lower()).strip("-")


def _layer_bars(d: pd.DataFrame, layer: str, stem: Path, formats, footer):
    ops = list(dict.fromkeys(d.operation))
    log = spans(d)
    heights = [0.7 + 0.3 * (d.operation == op).sum() for op in ops]
    fig = plt.figure(figsize=(FIG_W, sum(heights) + 0.5), layout="constrained")
    axes = fig.subplots(len(ops), 1, height_ratios=heights, squeeze=False)
    for ax, op in zip(axes[:, 0], ops):
        o = d[d.operation == op]
        y = [-i for i in range(len(o))]
        # light only where a family has both, so that the shade contrasts
        # two measurements rather than recording how one was taken
        both = {f for f, m in o.groupby("family")["mode"] if m.nunique() > 1}
        light = [m == "batch" and f in both for f, m in zip(o.family, o["mode"])]
        shaded_barh(ax, y, o, light, 0.7)
        interval_boxes(
            ax,
            y,
            o.value_lo_ns,
            o.value_hi_ns,
            o.value_ns,
            [
                color(base, True) if mode == "batch" else SURFACE
                for base, mode in zip(o.base, o["mode"])
            ],
            height=0.4,
        )
        value_axis(ax, d, log)
        time_labels(ax, y, o, log)
        ax.set_yticks(y, [wrap_label(s, 40) for s in o.label])
        ax.set_ylim(-len(o) + 0.4, 0.6)
        ax.tick_params(axis="y", length=0)
        if op != layer:
            ax.set_title(op, loc="left", color=TEXT)
    fig.suptitle(layer, x=0, ha="left", fontsize=HEAD, color=TEXT)
    interval_caption(fig)
    alt = (
        f"{layer}: a bar per benchmark and its 95% confidence-interval box with an estimate line, a panel per "
        "operation" + (", on a log scale" if log else "") + "."
    )
    return save(fig, stem, formats, alt, footer)


def markers(d: pd.DataFrame, families: list) -> dict:
    """A marker per family of families that d measures, distinct among those
    sharing a hue."""
    bases = dict(zip(d.family, d.base))
    seen: dict = {}
    out = {}
    for f in families:
        if f not in bases:
            continue
        b = bases[f]
        out[f] = MARKERS[seen.get(b, 0) % len(MARKERS)]
        seen[b] = seen.get(b, 0) + 1
    return out


MARKERS = "os^Dv"


def lines(ax, d: pd.DataFrame, x: str, families: list, alone=None, mark=None):
    """A line per family of families over x, its interval a band. With alone,
    the name of a boolean column, the rows where it holds are a second line,
    dark and dashed, and the rest light and solid. mark, from markers() over
    every panel's rows, keeps a family's marker when a panel lacks others."""
    mark = mark or markers(d, families)
    for fam in families:
        f = d[d.family == fam]
        for one in (True, False) if alone else (None,):
            g = (f[f[alone] == one] if alone else f).sort_values(x)
            if g.empty:
                continue
            base = g.base.iloc[0]
            light = alone is not None and not one
            c = color(base, light)
            style = ":" if base == BASELINE else "--" if one else "-"
            ax.plot(
                g[x],
                g.value_ns,
                style,
                color=c,
                marker=mark[fam],
                markersize=3.5,
                markeredgecolor=color(base),
                linewidth=2 if light else 1.4,
            )
            ax.fill_between(
                g[x], g.value_lo_ns, g.value_hi_ns, color=c, alpha=0.2, linewidth=0
            )
    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_formatter(NullFormatter())
    axis(ax, True, "y")
    ax.grid(axis="x", which="major", color=GRID, linewidth=0.8)


def base_of(d: pd.DataFrame, family: str) -> str:
    return d[d.family == family].base.iloc[0]


def family_legend(fig, d: pd.DataFrame, families: list, extra=(), light=False):
    """A family's dark line, and with light its light line beside it."""
    mark = markers(d, families)

    def handle(f, shade):
        return Line2D(
            [],
            [],
            color=color(base_of(d, f), shade),
            linestyle=":" if base_of(d, f) == BASELINE else "-",
            linewidth=2 if shade else 1.4,
            marker=mark[f],
            markersize=3.5,
            markeredgecolor=color(base_of(d, f)),
        )

    handles = [
        (handle(f, False), handle(f, True)) if light else handle(f, False)
        for f in families
    ]
    labels = [f"{f} (baseline)" if base_of(d, f) == BASELINE else f for f in families]
    for h, label in extra:
        handles.append(h)
        labels.append(label)
    fig.legend(
        handles,
        labels,
        loc="outside lower center",
        ncol=4,
        handler_map={tuple: HandlerTuple(ndivide=None, pad=0.6)},
        handlelength=4 if light else 2,
    )


def leaders(d: pd.DataFrame) -> pd.DataFrame:
    """Each family's fastest row of d, the measure shown() ranks by."""
    return d.loc[d.groupby("family").value_ns.idxmin()].set_index("family")[
        ["value_ns", "value_lo_ns", "value_hi_ns"]
    ]


AXES = {
    "m": "sketch cells m, log scale",
    "d": "set difference d (sketch of 2d cells or more), log scale",
}


def _hash_sweep(d: pd.DataFrame, stem: Path, formats, footer):
    h = d[(d.operation == "hash to curve") & d.group.isin(["h2c", "hash_to_curve"])]
    # the coordinates and codecs of a curve hash alike: the map is the curve's
    h = h[[rules.curve(f) == f for f in h.family]]
    h = h[(h.group == "h2c") | (h.family == BASELINE)]
    h = h.assign(n=[n if m == "batch" else 1 for n, m in zip(h.batch_n, h["mode"])])
    # the representatives, and the fastest of each base family beside them
    fastest = per_base(leaders(h), dict(zip(h.family, h.base)))
    order = {b: i for i, b in enumerate(rules.FAMILY_ORDER)}
    families = sorted(
        {f for f in [*REPRESENTATIVES, *fastest] if (h.family == f).any()},
        key=lambda f: (order.get(base_of(h, f), len(order)), f),
    )
    drawn = [f for f in families if f != BASELINE]
    fig, ax = plt.subplots(figsize=(FIG_W, 4.6), layout="constrained")
    lines(ax, h, "n", drawn)
    extra = []
    if BASELINE in families:
        # measured one at a time only: a level to compare against
        x = h[h.family == BASELINE].value_ns.min()
        ax.axhline(x, color=color(BASELINE), linestyle=":", linewidth=1.4)
        ax.annotate(
            f"{BASELINE}, one at a time: {fmt_time(x)}",
            (1, x),
            xycoords=ax.get_yaxis_transform(),
            xytext=(0, 3),
            textcoords="offset points",
            ha="right",
            fontsize=SMALL,
            color=TEXT_2,
        )
        level = Line2D([], [], color=color(BASELINE), linestyle=":", linewidth=1.4)
        extra.append((level, f"{BASELINE} (baseline), one at a time"))
    ax.set_xlabel("batch size n (1: one at a time), log scale")
    ax.set_title("try-and-increment, time per item", loc="left", color=TEXT)
    family_legend(fig, h, drawn, extra)
    shows = of(len(families), h.family.nunique()) or "all shown"
    fig.suptitle(
        f"Hash to curve against batch size ({shows})",
        x=0,
        ha="left",
        fontsize=HEAD,
        color=TEXT,
    )
    alt = (
        "Hash to curve by try-and-increment: time per item against batch size "
        "on log axes, one at a time at n = 1; the representative curves and "
        "the fastest of each base family, a marker per curve within a hue, "
        "and the XOR baseline as a level, intervals as bands. "
        "The other maps are in the comparison maps figure, hash to addend "
        "and the steps of a hash in the table."
    )
    return save(fig, stem, formats, alt, footer)


def _maps(t: pd.DataFrame, stem: Path, formats, footer):
    """The comparison maps beside try-and-increment on the same curves and
    the fastest try-and-increment curves: a row per curve and method, one
    at a time dark and batched light."""
    ex = t[t.layer == rules.EXPERIMENTAL]
    h2c = t[(t.layer == "hash to curve") & (t.group == "h2c")]
    h2c = h2c[[rules.curve(f) == f for f in h2c.family]]
    largest = h2c[h2c["mode"] == "batch"].batch_n.max()
    h2c = h2c[(h2c["mode"] != "batch") | (h2c.batch_n == largest)]
    native = t[(t.layer == "hash to curve") & t.label.str.contains("try-and-increment")]
    native = native[native.operation == "hash to curve"]
    best = leaders(h2c)
    fastest = list(best.sort_values("value_ns").index[:TOP])
    curves = [*fastest, *(f for f in dict.fromkeys(ex.family) if f not in fastest)]
    rows = []
    for src, method in ((h2c, "try-and-increment"), (native, None), (ex, None)):
        for r in src[src.family.isin(curves)].itertuples():
            name = method or r.label.split("/", 1)[1].removesuffix(", batched")
            rows.append(
                {
                    "row": f"{r.family}: {name}",
                    "family": r.family,
                    "base": r.base,
                    "variant": "batched" if r.mode == "batch" else "one at a time",
                    "value_ns": r.value_ns,
                    "value_lo_ns": r.value_lo_ns,
                    "value_hi_ns": r.value_hi_ns,
                }
            )
    m = pd.DataFrame(rows)

    def rank(r):
        fam, method = r.split(": ", 1)
        return curves.index(fam), method != "try-and-increment", method

    order = sorted(dict.fromkeys(m.row), key=rank)
    fig, ax = plt.subplots(
        figsize=(FIG_W, 0.9 + 0.5 * len(order)), layout="constrained"
    )
    variants = ["one at a time", "batched"]
    pair_bars(ax, m, order, variants, log=True, row="row")
    titled(
        ax,
        "time per item",
        f"{variant_key(variants)}, try-and-increment in batches of {largest:g}",
    )
    fig.suptitle(
        "Comparison maps beside try-and-increment, time per item",
        x=0,
        ha="left",
        fontsize=HEAD,
        color=TEXT,
    )
    alt = (
        "Comparison maps (Pornin's binary map, Elligator 2, SSWU) on a log "
        "scale beside try-and-increment on the same curves and on the three "
        "curves where it is fastest: a row per curve and method, one at a time "
        "and batched as a dark and a light bar of the family's colour."
    )
    return save(fig, stem, formats, alt, footer)


def _seeds(d: pd.DataFrame, stem: Path, formats, footer):
    ops = list(dict.fromkeys(d.operation))
    log = spans(d.assign(value_lo_ns=d.value_ns, value_hi_ns=d.value_ns))
    ncol = 2
    grid = [ops[i : i + ncol] for i in range(0, len(ops), ncol)]
    heights = [
        0.8 + 0.28 * max(d[d.operation == o].family.nunique() for o in row)
        for row in grid
    ]
    fig = plt.figure(figsize=(FIG_W, sum(heights) + 0.6), layout="constrained")
    axes = fig.subplots(len(grid), ncol, height_ratios=heights, squeeze=False)
    for r, row in enumerate(grid):
        for c in range(ncol):
            ax = axes[r][c]
            if c >= len(row):
                ax.set_visible(False)
                continue
            o = d[d.operation == row[c]]
            fams = list(dict.fromkeys(o.family))
            y = [-fams.index(f) for f in o.family]
            ax.scatter(o.value_ns, y, s=14, color=[color(b) for b in o.base], zorder=3)
            axis(ax, log)
            ax.set_yticks([-i for i in range(len(fams))], fams)
            ax.set_ylim(-len(fams) + 0.4, 0.6)
            ax.tick_params(axis="y", length=0)
            ax.set_title(wrap_label(row[c], 40), loc="left", color=TEXT)
    fig.suptitle(
        "Curve generation, total time per seed",
        x=0,
        ha="left",
        fontsize=HEAD,
        color=TEXT,
    )
    alt = "Curve generation: a dot per seed for each family, a panel per stage" + (
        ", on log scales." if log else "."
    )
    return save(fig, stem, formats, alt, footer)


# Curve selection, end to end (curvegen.csv).

# the phases curvegen_times records, in search order, in the Tableau 10
# hues that no family in the figure has; find has no factor or witness phase
PHASES = [
    ("candidate", "#9c755f"),
    ("quick", "#76b7b2"),
    ("sieve", "#ff9da7"),
    ("count", "#e15759"),
    ("factor", "#edc948"),
    ("witness", "#b07aa1"),
    ("accept", "#bab0ac"),
]
UNTIMED = SURFACE


def plot_selection(runs: pd.DataFrame, out: Path, formats, footer="") -> Figure:
    """Each family's curve selection by one method, a row per seed in runs:
    finding a curve and proving it, the seeds' mean as a bar and their range
    as a box; verifying the certificate; and the prover's time by phase."""
    s = runs.assign(base=[rules.family(f)[1] for f in runs.family])
    fams = list(s.groupby("family").prove_s.mean().sort_values().index)
    long = pd.concat(
        [
            s.assign(variant=v, value=s[f"{v}_s"] * 1e9)
            for v in ("find", "prove", "verify")
        ]
    )
    d = (
        long.groupby(["family", "variant", "base"])
        .value.agg(value_ns="mean", value_lo_ns="min", value_hi_ns="max")
        .reset_index()
    )
    fig = plt.figure(figsize=(FIG_W, 1.2 + 0.42 * len(fams)), layout="constrained")
    search, check, share = fig.subplots(
        1, 3, sharey=True, gridspec_kw={"width_ratios": [2.2, 1.2, 1.6]}
    )
    pair_bars(search, d, fams, ["find", "prove"], True, row="family")
    titled(search, "find (dark) and prove (light)", "per seed: mean bar, range box")
    pair_bars(check, d, fams, ["verify"], True, row="family")
    # a narrow panel: decades only
    check.xaxis.set_major_locator(LogLocator(base=10, subs=(1.0,)))
    titled(check, "verify the certificate")
    mean = s.groupby("family").mean(numeric_only=True)
    for i, f in enumerate(fams):
        shares = [
            (mean.at[f, f"prove_{p}_s"] / mean.at[f, "prove_s"], c)
            for p, c in PHASES
            if f"prove_{p}_s" in mean
        ]
        shares.append((1 - sum(w for w, _ in shares), UNTIMED))
        left = 0.0
        for w, c in shares:
            if w > 0:
                share.barh(
                    -i,
                    w,
                    left=left,
                    height=0.7,
                    color=c,
                    edgecolor=TEXT_2 if c == UNTIMED else SURFACE,
                    linewidth=0.6,
                )
                left += w
    share.set_xlim(0, 1)
    share.xaxis.set_major_formatter(PercentFormatter(1))
    share.tick_params(axis="y", length=0)
    titled(share, "the prover's time by phase")
    present = [
        (p, c) for p, c in PHASES if (s.get(f"prove_{p}_s", pd.Series([0])) > 0).any()
    ]
    fig.legend(
        [Patch(facecolor=c) for _, c in present]
        + [Patch(facecolor=UNTIMED, edgecolor=TEXT_2)],
        [p for p, _ in present] + ["outside the phases"],
        loc="outside lower center",
        ncol=len(present) + 1,
    )
    methods = ", ".join(sorted(set(s.method)))
    fig.suptitle(
        f"Curve selection, end to end, over {s.seed.nunique()} seeds ({methods})",
        x=0,
        ha="left",
        fontsize=HEAD,
        color=TEXT,
    )
    alt = (
        "Curve selection per family on log time axes: the mean time to find a "
        "curve (dark) and to find and certify one (light), with a box from the "
        "fastest to the slowest seed; the mean time to verify the certificate; "
        "and, as shares of the mean proving time, the time spent deriving "
        "candidates, in quick rejection, sieving, point counting, factoring, "
        "rejection witnesses and acceptance."
    )
    return save(fig, out / "selection", formats, alt, footer)


# The per-item cost of insertion.

# Tableau hues that no family is drawn in, so that a part is not read as one
