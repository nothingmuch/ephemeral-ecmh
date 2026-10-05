"""Assemble the mdbook source from README.md and the docs/ bundle.

docs/index.md, the bundle's OKF index, states the chapter order: its
headings become the parts of SUMMARY.md and its entries the chapters,
without their descriptions. The README is the introduction, and the index,
with its descriptions, is the contents page to which links to docs/index.md
lead. Frontmatter is dropped, as mdbook would render it as text, and links
that leave docs/ point at the repository at the built revision. Each key
cited links to its entry in the annotated bibliography, docs/literature.md,
whose entries are headed by their citations as formatted from
references.bib.

    python book.py ROOT OUT --repository URL --rev REV
"""

import argparse
import posixpath
import re
from collections.abc import Container
from pathlib import Path

import bibliography

FRONTMATTER = re.compile(r"\A---\n.*?\n---\n+", re.DOTALL)
LINK = re.compile(r"\]\(([^)\s]+)\)")
SCHEME = re.compile(r"[A-Za-z][A-Za-z0-9+.-]*:|#")
ENTRY = re.compile(r"\* \[([^\]]+)\]\(([^)\s]+)\)(?: - .*)?")
CONTENTS = "contents.md"
LITERATURE = "literature.md"
HEAD = re.compile(bibliography.ANNOTATION.pattern + r"\s*", re.MULTILINE)


def strip_frontmatter(text: str) -> str:
    return FRONTMATTER.sub("", text, count=1)


def relink(text: str, base: str, blob: str) -> str:
    """Links of a file in directory `base` of the repository, rewritten for
    a book whose source is docs/: targets in docs/ become book-relative,
    and every other target a URL under `blob`."""

    def sub(m: re.Match) -> str:
        url = m.group(1)
        if SCHEME.match(url):
            return m.group(0)
        path, sep, frag = url.partition("#")
        target = posixpath.normpath(posixpath.join(base, path))
        if target == "docs/index.md":
            return f"]({CONTENTS}{sep}{frag})"
        if target.startswith("docs/"):
            return f"]({target.removeprefix('docs/')}{sep}{frag})"
        return f"]({blob}/{target}{sep}{frag})"

    return LINK.sub(sub, text)


def link_citations(text: str, keys: Container[str]) -> str:
    return bibliography.link_citations(
        text, lambda k: f"{LITERATURE}#{k}" if k in keys else None
    )


def literature(text: str, entries: dict[str, bibliography.Entry]) -> str:
    """The annotated bibliography with each entry's head of keys replaced
    by an anchor and the formatted citation per key, the annotation
    following as a paragraph of the same item."""

    def head(key: str) -> str:
        cite = bibliography.render(entries[key], lambda k: f"#{k}")
        return f'<a id="{key}"></a>\\[{key}\\] {cite}'

    def sub(m: re.Match) -> str:
        if text.startswith("- ", m.end()) or m.end() == len(text):
            raise ValueError(f"no annotation after {m.group(0)!r}")
        heads = [head(k) for k in bibliography.annotated_keys(m.group(1))]
        return "- " + "\\\n  ".join(heads) + "\n\n  "

    return HEAD.sub(sub, text)


def summary(index: str) -> str:
    """SUMMARY.md from an OKF index: `# Part` headings and
    `* [Title](file.md) - description` entries."""
    out = ["# Summary", "", "[Introduction](README.md)", f"[Contents]({CONTENTS})"]
    for line in strip_frontmatter(index).splitlines():
        if line.startswith("# "):
            out += ["", line, ""]
        elif m := ENTRY.fullmatch(line):
            out.append(f"- [{m.group(1)}]({m.group(2)})")
        elif line.strip():
            raise ValueError(f"unexpected line in index: {line!r}")
    return "\n".join(out) + "\n"


def contents(index: str) -> str:
    """The contents page: the index with its descriptions, its parts one
    heading level down."""
    body = strip_frontmatter(index)
    body = re.sub(r"^# ", "## ", body, flags=re.MULTILINE)
    return "# Contents\n\n" + body


# mdbook's 750px column leaves most of a wide screen empty and the grids and
# figures cramped; 900px still keeps the prose near a hundred characters
CSS = ":root{--content-max-width:900px}\n"


def book_toml(repository: str) -> str:
    return (
        "[book]\n"
        'title = "Ephemeral ECMH"\n'
        'language = "en"\n'
        'src = "src"\n'
        "\n"
        "[output.html]\n"
        f'git-repository-url = "{repository}"\n'
        'additional-css = ["book.css"]\n'
        # GitHub Pages serves a repository's site under /NAME/.
        f'site-url = "/{repository.rstrip("/").rsplit("/", 1)[-1]}/"\n'
        # $...$ is the math syntax GitHub renders, so the chapters read the
        # same there and in the book
        "\n"
        "[preprocessor.katex]\n"
        'after = ["links"]\n'
    )


def assemble(root: Path, out: Path, repository: str, rev: str) -> None:
    blob = f"{repository}/blob/{rev}"
    src = out / "src"
    src.mkdir(parents=True)
    (out / "book.toml").write_text(book_toml(repository))
    (out / "book.css").write_text(CSS)
    docs = root / "docs"
    index = (docs / "index.md").read_text()
    (src / "SUMMARY.md").write_text(summary(index))
    (src / CONTENTS).write_text(contents(index))
    entries = bibliography.parse((root / "references.bib").read_text())
    readme = (root / "README.md").read_text()
    (src / "README.md").write_text(link_citations(relink(readme, ".", blob), entries))
    for page in sorted(docs.glob("*.md")):
        # docs/README.md directs readers of the repository to the book; as a
        # chapter it would overwrite the introduction.
        if page.name in ("index.md", "README.md"):
            continue
        text = strip_frontmatter(page.read_text())
        text = link_citations(relink(text, "docs", blob), entries)
        if page.name == LITERATURE:
            text = literature(text, entries)
        (src / page.name).write_text(text)


def main(argv=None) -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("root", type=Path)
    p.add_argument("out", type=Path)
    p.add_argument("--repository", required=True)
    p.add_argument("--rev", required=True)
    a = p.parse_args(argv)
    assemble(a.root, a.out, a.repository, a.rev)


if __name__ == "__main__":
    main()
