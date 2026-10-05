"""Check the citations of README.md and docs/ against references.bib.

Citations are recognized as bibliography.py describes; the frontmatter is
not scanned. Every cited key must be in the bibliography, and
docs/literature.md, from which the book renders each entry, must annotate
every key of the bibliography exactly once.

    python bibcheck.py ROOT
"""

import sys
from pathlib import Path

import bibliography
from bibliography import ANNOTATION, KEY
from book import strip_frontmatter


def citations(text: str, keys: set[str]) -> tuple[list[str], list[str]]:
    """The keys a chapter cites, and the keys of its citations that are not
    in the bibliography."""
    cited, bad = [], []
    for parts in bibliography.citations(strip_frontmatter(text)):
        for p in parts:
            if p in keys:
                cited.append(p)
            elif KEY.fullmatch(p):
                bad.append(p)
    return cited, bad


def annotations(text: str) -> list[str]:
    return [
        k
        for m in ANNOTATION.finditer(strip_frontmatter(text))
        for k in bibliography.annotated_keys(m.group(1))
    ]


def _once(where: str, found: list[str], expected: set[str]) -> list[str]:
    out = []
    for k in sorted(expected - set(found)):
        out.append(f"{where}: missing {k}")
    for k in sorted(set(found) - expected):
        out.append(f"{where}: unexpected {k}")
    for k in sorted({k for k in found if found.count(k) > 1}):
        out.append(f"{where}: {k} listed more than once")
    return out


def problems(root: Path) -> list[str]:
    keys = set(bibliography.parse((root / "references.bib").read_text()))
    docs = root / "docs"
    out = []
    for page in [root / "README.md", *sorted(docs.glob("*.md"))]:
        if page.name == "literature.md":
            continue
        _, bad = citations(page.read_text(), keys)
        for k in bad:
            out.append(f"{page.relative_to(root)}: [{k}] is not in references.bib")
    out += _once(
        "docs/literature.md", annotations((docs / "literature.md").read_text()), keys
    )
    return out


def main(argv=None) -> None:
    (root,) = argv if argv is not None else sys.argv[1:]
    errors = problems(Path(root))
    for e in errors:
        print(e, file=sys.stderr)
    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()
