"""Comparison tables: rows grouped under headings, columns under shared
headers, and in each comparable column the best value marked and the rest
shaded by how far they are from it. A best that other values tie within
their confidence intervals is no win: those cells are marked as a tie
instead. Rendered to HTML and to markdown.
"""

import html
import math
from dataclasses import dataclass, field
from itertools import pairwise

# kinds of column: a time (lower is better), a factor (higher is better), or
# a time or a number shown for reading, never compared
TIME, FACTOR, PLAIN_TIME, PLAIN_FACTOR, PLAIN = (
    "time",
    "factor",
    "plain-time",
    "plain-factor",
    "plain",
)
# a column's shading saturates at this ratio to its best
SHADE_RATIO = 16
# every compared cell is coloured by its ratio to the column's best, on a log
# scale along matplotlib's inferno, sampled at 1.0, 0.9, ..., 0.3: pale yellow
# at the best, through orange and red, to purple at SHADE_RATIO. inferno's
# lightness is monotone, so the order survives greyscale and colour blindness
SHADE_STOPS = [
    (i / 7, c)
    for i, c in enumerate(
        [
            (252, 255, 164),
            (246, 215, 70),
            (252, 165, 10),
            (243, 120, 25),
            (221, 81, 58),
            (188, 55, 84),
            (147, 38, 103),
            (106, 23, 110),
        ]
    )
]
# under every grid: the marks are conventions the reader meets on each table,
# and a tooltip needs a pointer; the cells repeat the same words for the ear
KEY = (
    "<b>bold</b>: the best of its column; <i>italic, dotted</i>: tied with "
    "the best within the confidence intervals; shading: ratio to the best, "
    f"pale yellow at 1× to purple at {SHADE_RATIO}× or more; grey: baseline, "
    "not compared."
)


@dataclass(frozen=True)
class Col:
    key: object
    head: str
    sub: str = ""
    kind: str = TIME
    note: str = ""


@dataclass
class Row:
    label: str
    note: str = ""
    values: dict = field(default_factory=dict)
    tips: dict = field(default_factory=dict)
    # False: shown, but neither best nor shaded (the XOR baseline)
    compare: bool = True
    # what sets the row apart from the others it's compared with, {mark: why}
    marks: dict = field(default_factory=dict)
    # a time's confidence interval (lo, hi), by column key
    cis: dict = field(default_factory=dict)


@dataclass
class Grid:
    cols: list
    groups: list  # [(heading or "", [Row])]
    corner: str = ""

    def rows(self):
        return [r for _, rows in self.groups for r in rows]

    def shown(self) -> list:
        """The columns with a value in some row: a column blank throughout,
        like log2 rho for a run that recorded no r, is left out rather than
        shown empty."""
        rows = self.rows()
        return [c for c in self.cols if any(_filled(r.values.get(c.key)) for r in rows)]

    def best(self, col: Col):
        """The column's best value among the compared rows, or None."""
        if col.kind not in (TIME, FACTOR):
            return None
        vs = [v for r in self.rows() if r.compare and _num(v := r.values.get(col.key))]
        if not vs:
            return None
        return min(vs) if col.kind == TIME else max(vs)

    def leaders(self, col: Col) -> set:
        """The ids of the compared rows the best can't be told from: those
        within rounding of it, and for times, those whose interval overlaps
        the best's. One leader is the best; more are a tie."""
        b = self.best(col)
        if b is None:
            return set()
        rows = [r for r in self.rows() if r.compare and _num(r.values.get(col.key))]
        top = [r for r in rows if ratio(col, r.values[col.key], b) <= 1.0005]
        ids = {id(r) for r in top}
        if col.kind == TIME:
            hi = max(r.cis.get(col.key, (b, b))[1] for r in top)
            ids |= {id(r) for r in rows if r.cis.get(col.key, (math.inf,))[0] <= hi}
        return ids


def _num(v) -> bool:
    return isinstance(v, (int, float)) and not math.isnan(v)


def _filled(v) -> bool:
    return _num(v) or (v is not None and not isinstance(v, float) and str(v) != "")


def fmt(col: Col, v, fmt_time) -> str:
    if not _num(v):
        return "" if v is None or (isinstance(v, float) and math.isnan(v)) else str(v)
    if col.kind in (TIME, PLAIN_TIME):
        return fmt_time(v)
    if col.kind in (FACTOR, PLAIN_FACTOR):
        return f"{v:.1f}×" if v < 100 else f"{v:.0f}×"
    v = 0.0 if abs(v) < 0.05 else v  # no "-0"
    return f"{v:.1f}" if v != int(v) else f"{v:.0f}"


def ratio(col: Col, v, best) -> float:
    """How many times worse than the best: 1 for the best itself."""
    if not (_num(v) and _num(best)) or best <= 0 or v <= 0:
        return float("nan")
    return v / best if col.kind == TIME else best / v


def _shade(r: float) -> str:
    """The cell's style: its colour, and black or white text, whichever
    contrasts more with it, so that the cell reads alike on any page."""
    t = min(1.0, math.log(r) / math.log(SHADE_RATIO)) if r > 1 else 0.0
    for (t0, c0), (t1, c1) in pairwise(SHADE_STOPS):
        if t <= t1:
            f = (t - t0) / (t1 - t0)
            c = [round(a + (b - a) * f) for a, b in zip(c0, c1)]
            # WCAG's crossover: black and white contrast equally at this luminance
            text = "#000" if _luminance(c) > 0.179 else "#fff"
            return "background:rgb({},{},{});color:{}".format(*c, text)
    raise AssertionError


def _luminance(rgb) -> float:
    """WCAG relative luminance of an sRGB colour."""
    lin = [
        v / 12.92 if v <= 0.04045 else ((v + 0.055) / 1.055) ** 2.4
        for v in (x / 255 for x in rgb)
    ]
    return 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]


def _headers(cols):
    """The top header's runs: (head, span, has subs)."""
    runs = []
    for c in cols:
        if runs and runs[-1][0] == c.head and c.sub:
            runs[-1][1] += 1
        else:
            runs.append([c.head, 1, bool(c.sub)])
    return runs


def to_html(g: Grid, fmt_time) -> str:
    esc = html.escape
    cols = g.shown()
    two = any(c.sub for c in cols)
    out = ['<table class="grid"><thead><tr>']
    out.append(f'<th class="corner"{" rowspan=2" if two else ""}>{esc(g.corner)}</th>')
    i = 0
    for head, span, subs in _headers(cols):
        note = cols[i].note
        title = f' title="{esc(note)}"' if note else ""
        if subs:
            out.append(f'<th class="head" colspan="{span}"{title}>{esc(head)}</th>')
        else:
            out.append(
                f'<th class="head"{" rowspan=2" if two else ""}{title}>{esc(head)}</th>'
            )
        i += span
    out.append("</tr>")
    if two:
        out.append("<tr>")
        out += [f'<th class="sub">{esc(c.sub)}</th>' for c in cols if c.sub]
        out.append("</tr>")
    out.append("</thead>")
    best = {c.key: g.best(c) for c in cols}
    lead = {c.key: g.leaders(c) for c in cols}
    for heading, rows in g.groups:
        out.append("<tbody>")
        if heading:
            out.append(
                f'<tr class="group"><th colspan="{len(cols) + 1}">{esc(heading)}</th></tr>'
            )
        for r in rows:
            note = f'<span class="note">{esc(r.note)}</span>' if r.note else ""
            tags = "".join(
                f'<span class="mark" title="{esc(why)}">{esc(m)}</span>'
                for m, why in r.marks.items()
            )
            out.append(f'<tr><th class="row">{esc(r.label)}{tags}{note}</th>')
            for c in cols:
                v = r.values.get(c.key)
                text = fmt(c, v, fmt_time)
                tip = r.tips.get(c.key, "")
                cls, style, said = [c.kind], "", ""
                b = best[c.key]
                if r.compare and _num(v) and b is not None:
                    q = ratio(c, v, b)
                    style = f' style="{_shade(q)}"'
                    if id(r) in lead[c.key] and len(lead[c.key]) == 1:
                        cls.append("best")
                        said = "best"
                    elif id(r) in lead[c.key]:
                        cls.append("tie")
                        said = "tied for the best, within noise"
                        if q > 1.0005:
                            said += f", {q:.2f}×"
                    elif q > 1.0005:
                        said = f"{q:.2f}× the best"
                    if c.key in r.cis:
                        lo, hi = r.cis[c.key]
                        tip = f"{fmt_time(lo)} to {fmt_time(hi)}" + (
                            f"; {tip}" if tip else ""
                        )
                    tip = f"{said}; {tip}" if tip else said
                elif not r.compare:
                    cls.append("baseline")
                title = f' title="{esc(tip)}"' if tip else ""
                # the mark's meaning as text too: tooltips need a pointer, and
                # bold or shading say nothing to a screen reader
                aloud = f'<span class="sr"> ({esc(said)})</span>' if said else ""
                out.append(
                    f'<td class="{" ".join(cls)}"{style}{title}>{esc(text)}{aloud}</td>'
                )
            out.append("</tr>")
        out.append("</tbody>")
    out.append("</table>")
    out.append(f'<p class="key">{KEY}</p>')
    return "".join(out)


def to_markdown(g: Grid, fmt_time, cell) -> str:
    """One header row (head, sub), headings as bold rows, the best in bold
    and a tie for it in italics."""
    cols = g.shown()
    heads = [g.corner] + [f"{c.head}, {c.sub}" if c.sub else c.head for c in cols]
    lines = ["| " + " | ".join(cell(h) for h in heads) + " |"]
    lines.append("|" + "---|" * len(heads))
    lead = {c.key: g.leaders(c) for c in cols}
    for heading, rows in g.groups:
        if heading:
            lines.append(f"| **{cell(heading)}** |" + " |" * len(cols))
        for r in rows:
            label = cell(r.label) + "".join(f" [{cell(m)}]" for m in r.marks)
            label += f" ({cell(r.note)})" if r.note else ""
            cells = [label]
            for c in cols:
                v = r.values.get(c.key)
                text = cell(fmt(c, v, fmt_time))
                if id(r) in lead[c.key]:
                    text = f"*{text}*" if len(lead[c.key]) > 1 else f"**{text}**"
                cells.append(text)
            lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


# in em, not rem: the book's root is 10px, the standalone report's 16px
STYLE = (
    "table.grid{border-collapse:separate;border-spacing:0;font-size:.8em;"
    "margin:.625em 0 1.25em}"
    "table.grid th,table.grid td{padding:.25em .69em;border-bottom:1px solid #ecebe7}"
    "table.grid thead th{font-weight:600;text-align:center;border-bottom:1px solid #b9b8b3}"
    "table.grid th.head{border-left:1px solid #ecebe7}"
    "table.grid th.sub{font-weight:400;color:#52514e;font-size:.94em}"
    "table.grid th.row{text-align:left;font-weight:500;white-space:nowrap}"
    "table.grid th.row .note{display:block;font-weight:400;color:#6b6a66;"
    "font-size:.875em;white-space:normal;max-width:27.5em}"
    "table.grid tr.group th{text-align:left;background:#f1f0ec;color:#2b2b29;"
    "font-weight:600;padding-top:.44em}"
    "table.grid td{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}"
    "table.grid td.best{font-weight:700}"
    "table.grid td.tie{font-style:italic;text-decoration:underline dotted}"
    "table.grid td.baseline{color:#6b6a66}"
    "table.grid th.row .mark{margin-left:.5em;padding:0 .375em;border-radius:.25em;"
    "background:#f1e4c8;color:#6b4e16;font-size:.8125em;font-weight:500}"
    # visually hidden, read aloud: the usual clip recipe, not display:none,
    # which assistive technology skips
    ".sr{position:absolute;width:1px;height:1px;overflow:hidden;"
    "clip:rect(0 0 0 0);clip-path:inset(50%);white-space:nowrap}"
    "p.key{font-size:.8em;margin:-.5em 0 1.5em}"
)
