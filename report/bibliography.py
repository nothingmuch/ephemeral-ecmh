"""references.bib, its entries formatted as markdown, and the citations of
the chapters.

A citation is a bracket group of keys and locators, `[key]` or
`[key, key, Section 4; key, Section 5]`: a key is a lowercase word, a
locator starts with a capital letter, and a group with no key is not a
citation. Code spans and fenced blocks hold no citations. The annotated
bibliography, docs/literature.md, heads each entry with the backticked keys
it annotates.
"""

import re
import unicodedata
from collections.abc import Callable, Iterator

CODE = re.compile(r"(```.*?```|`[^`\n]*`)", re.DOTALL)
GROUP = re.compile(r"\[([^\]\[]+)\](?![(\[])")
KEY = re.compile(r"[a-z][a-z0-9-]*")
LOCATOR = re.compile(r"[A-Z]")
ANNOTATION = re.compile(r"^- ((?:`[^`]+`[,;]?\s+)*`[^`]+`)\.", re.MULTILINE)

Entry = tuple[str, dict[str, str]]


def citation_parts(group: str) -> list[str] | None:
    """The keys and locators of a bracket group's contents, or None if the
    group is not a citation."""
    parts = [p.strip() for p in re.split(r"[,;]", group)]
    if not all(KEY.fullmatch(p) or LOCATOR.match(p) for p in parts):
        return None
    if not any(KEY.fullmatch(p) for p in parts):
        return None
    return parts


def citations(text: str) -> Iterator[list[str]]:
    """The parts of each citation in markdown text."""
    for i, chunk in enumerate(CODE.split(text)):
        if i % 2 == 0:
            for m in GROUP.finditer(chunk):
                if (parts := citation_parts(m.group(1))) is not None:
                    yield parts


def link_citations(text: str, href: Callable[[str], str | None]) -> str:
    """Markdown text whose citations link each key for which `href` gives a
    target. The separators and locators of a citation are kept."""

    def sub(m: re.Match) -> str:
        if citation_parts(m.group(1)) is None:
            return m.group(0)

        def part(p: re.Match) -> str:
            url = href(p.group(0)) if KEY.fullmatch(p.group(0)) else None
            return f"[{p.group(0)}]({url})" if url else p.group(0)

        return "[" + re.sub(r"[^,;\s][^,;]*?(?=\s*(?:[,;]|$))", part, m.group(1)) + "]"

    chunks = CODE.split(text)
    return "".join(GROUP.sub(sub, c) if i % 2 == 0 else c for i, c in enumerate(chunks))


def annotated_keys(head: str) -> list[str]:
    return re.findall(r"`([^`]+)`", head)


# -- parsing ----------------------------------------------------------------


def _braced(src: str, i: int) -> tuple[str, int]:
    """The contents of the brace group opening at src[i], and the index past
    its close."""
    depth = 0
    for j in range(i, len(src)):
        if src[j] == "{":
            depth += 1
        elif src[j] == "}":
            depth -= 1
            if depth == 0:
                return src[i + 1 : j], j + 1
    raise ValueError(f"unbalanced braces at offset {i}")


ENTRY_START = re.compile(r"@(\w+)\s*\{\s*([^,\s]+)\s*,")
SEPARATOR = re.compile(r"(?:[\s,]|%[^\n]*)*")
FIELD = re.compile(r"\s*(\w+)\s*=\s*")
BARE = re.compile(r"[\w-]+")


def parse(bib: str) -> dict[str, Entry]:
    """Entries by key, in file order: the entry type and the fields, with
    lowercase names and values as written."""
    out: dict[str, Entry] = {}
    pos = 0
    while m := ENTRY_START.search(bib, pos):
        kind, key = m.group(1).lower(), m.group(2)
        if key in out:
            raise ValueError(f"{key}: duplicate key")
        fields: dict[str, str] = {}
        i = m.end()
        while True:
            i = SEPARATOR.match(bib, i).end()
            if bib[i] == "}":
                break
            f = FIELD.match(bib, i)
            if not f:
                raise ValueError(f"{key}: cannot parse field at offset {i}")
            i = f.end()
            if bib[i] == "{":
                value, i = _braced(bib, i)
            elif bib[i] == '"':
                j = bib.index('"', i + 1)
                value, i = bib[i + 1 : j], j + 1
            else:
                v = BARE.match(bib, i)
                if not v:
                    raise ValueError(f"{key}: cannot parse value at offset {i}")
                value, i = v.group(0), v.end()
            fields[f.group(1).lower()] = value
        out[key] = (kind, fields)
        pos = i + 1
    return out


# -- LaTeX to markdown ------------------------------------------------------

ACCENT = {
    "'": "\u0301",
    "`": "\u0300",
    '"': "\u0308",
    "^": "\u0302",
    "~": "\u0303",
    "=": "\u0304",
    ".": "\u0307",
    "c": "\u0327",
}
LETTER = {"i": "i", "aa": "å", "AA": "Å", "ae": "æ", "AE": "Æ", "o": "ø", "O": "Ø"}
LETTER |= {"ss": "ß", "oe": "œ", "OE": "Œ", "l": "ł", "L": "Ł"}
COMMAND = re.compile(r"\\(url|cite)\{([^}]*)\}")
MARKDOWN = re.compile(r"([\\`*_\[\]<>|])")
DATED = ("journal", "booktitle", "howpublished", "number", "note")
MONTHS = [
    "January",
    "February",
    "March",
    "April",
    "May",
    "June",
    "July",
    "August",
    "September",
    "October",
    "November",
    "December",
]


def _accents(s: str) -> str:
    s = re.sub(
        r"\\([a-zA-Z]+)\b(?:\{\})?",
        lambda m: LETTER.get(m.group(1), m.group(0)),
        s,
    )
    s = re.sub(
        r"\\([\'`\"^~=.c])\s*(?:\{([^{}\\]?)\}|([A-Za-z]))",
        lambda m: (m.group(2) or m.group(3) or "") + ACCENT[m.group(1)],
        s,
    )
    return unicodedata.normalize("NFC", s)


def _plain(s: str) -> str:
    """A text field without markup: accents and escapes resolved, math left
    as written without its delimiters, protective braces dropped."""
    s = re.sub(r"\\mathbb\{(\w)\}", r"\1", s)
    # Escaped symbols are set aside so that brace and $ removal spare them.
    s = re.sub(r"\\([&#_%${}])", lambda m: "\0" + str(ord(m.group(1))) + "\0", s)
    s = _accents(s)
    while (t := re.sub(r"(?<!_)\{([^{}]*)\}", r"\1", s)) != s:
        s = t
    s = s.replace("$", "")
    s = s.replace("---", "\u2014").replace("--", "\u2013")
    s = s.replace("``", "\u201c").replace("''", "\u201d")
    s = re.sub(r"\s+", " ", s).strip()
    s = re.sub(r"\0(\d+)\0", lambda m: chr(int(m.group(1))), s)
    if "\\" in s:
        raise ValueError(f"unsupported LaTeX in {s!r}")
    return s


def _escape(s: str) -> str:
    return MARKDOWN.sub(r"\\\1", s)


def _text(s: str, href: Callable[[str], str]) -> str:
    """A text field as markdown: `\\url` becomes an autolink and `\\cite` a
    citation."""
    out = []
    for i, part in enumerate(COMMAND.split(s)):
        if i % 3 == 0:
            out.append(
                re.sub(
                    r"\S(.*\S)?",
                    lambda m: _escape(_plain(m.group(0))),
                    part,
                    flags=re.DOTALL,
                )
            )
        elif i % 3 == 1:
            command = part
        elif command == "url":
            out.append(f"<{part}>")
        else:
            out.append(f"[[{part}]({href(part)})]")
    return re.sub(r"\s+", " ", "".join(out)).strip()


# -- names ------------------------------------------------------------------


def _split_and(s: str) -> list[str]:
    out, depth, start = [], 0, 0
    for m in re.finditer(r"[{}]|\s+and\s+", s):
        if m.group(0) == "{":
            depth += 1
        elif m.group(0) == "}":
            depth -= 1
        elif depth == 0:
            out.append(s[start : m.start()])
            start = m.end()
    out.append(s[start:])
    return [a.strip() for a in out]


def _initials(given: str) -> str:
    """`Daniel J.` as `D. J.`, `Jean-Pierre` as `J.-P.`."""
    out = []
    for word in given.split():
        parts = [p for p in word.split("-") if p]
        out.append("-".join(p[0] + "." for p in parts))
    return " ".join(out)


def name(author: str) -> str:
    """One author as initials and surname; a braced name is a corporate
    author and kept whole."""
    if author.startswith("{") and _braced(author, 0)[1] == len(author):
        return _plain(author)
    if author == "others":
        return "et al."
    parts = [p.strip() for p in _plain(author).split(",")]
    if len(parts) == 3:
        surname, suffix, given = parts
        surname = f"{surname}, {suffix}"
    elif len(parts) == 2:
        surname, given = parts
    else:
        words = parts[0].split()
        # BibTeX's `First von Last`: the surname starts at the first
        # lowercase word, else it is the last word.
        cut = next(
            (i for i, w in enumerate(words[:-1]) if i and w[0].islower()),
            len(words) - 1,
        )
        given, surname = " ".join(words[:cut]), " ".join(words[cut:])
    return f"{_initials(given)} {surname}".strip()


def authors(field: str) -> str:
    return ", ".join(name(a) for a in _split_and(field))


# -- entries ----------------------------------------------------------------


def _pages(p: str) -> str:
    p = _plain(p)
    return f"pp. {p}" if re.search("[–,]", p) else f"p. {p}"


def link(fields: dict[str, str]) -> str | None:
    if doi := fields.get("doi"):
        return f"https://doi.org/{doi}"
    return fields.get("url") or fields.get("eprint")


def render(entry: Entry, href: Callable[[str], str]) -> str:
    """One entry as a markdown sentence sequence: authors, title, venue and
    date, link. `href` gives the target of a key cited in a note."""
    kind, f = entry

    def t(name: str) -> str | None:
        return _text(f[name], href) if name in f else None

    venue: list[str | None] = []
    if kind == "article":
        journal = t("journal")
        if "volume" in f:
            journal += " " + _plain(f["volume"])
            if "number" in f:
                journal += f"({_plain(f['number'])})"
        venue = [journal]
    elif kind in ("inproceedings", "incollection"):
        container = t("booktitle")
        if kind == "incollection":
            container = "In: " + container
            if "editor" in f:
                eds = authors(f["editor"])
                container += f", ed. {eds}"
        venue = [container]
    elif kind == "techreport":
        venue = [t("institution"), t("number")]
    elif kind == "book":
        venue = []
    else:
        venue = [t("howpublished")]
    if "series" in f:
        series = t("series")
        venue.append(f"{series} {_plain(f['volume'])}" if "volume" in f else series)
    if kind in ("book", "incollection"):
        venue.append(t("publisher"))
    if "pages" in f:
        venue.append(_pages(f["pages"]))
    venue = [v for v in venue if v]
    year = f.get("year")
    named = " ".join(_plain(COMMAND.sub(r"\2", f[k])) for k in DATED if k in f)
    if year and not re.search(rf"\b{year}\b", named):
        month = f.get("month")
        month = next((m for m in MONTHS if month and m.lower().startswith(month)), None)
        venue.append(f"{month} {year}" if month else year)
    out = []
    if "author" in f:
        out.append(authors(f["author"]))
    out.append(t("title"))
    if venue:
        out.append(", ".join(venue))
    if note := t("note"):
        out.append(note)
    out = [s if s.endswith((".", "?", "!")) else s + "." for s in out]
    if url := link(f):
        out.append(f"<{url}>")
    return " ".join(out)
