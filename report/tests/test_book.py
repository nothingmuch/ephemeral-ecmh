import json
import re

import bench_report
import bibliography
import book
import figures
import pytest
from criterion_fixture import write_bench

BLOB = "https://example.org/r/blob/abc"
BIB = (
    "@misc{a-2001,\n  author = {Doe, Jane},\n  title = {A},\n  year = {2001},\n}\n"
    "@misc{b-2002,\n  title = {B},\n  note = {After \\cite{a-2001}},\n}\n"
)


def test_frontmatter_is_dropped_once():
    text = "---\ntype: Chapter\n---\n\n# T\n\n---\n\nbody\n"
    assert book.strip_frontmatter(text) == "# T\n\n---\n\nbody\n"
    assert book.strip_frontmatter("# T\n") == "# T\n"


@pytest.mark.parametrize(
    "base, link, expected",
    [
        ("docs", "literature.md", "literature.md"),
        ("docs", "literature.md#wagner", "literature.md#wagner"),
        ("docs", "../src/hash.rs", f"{BLOB}/src/hash.rs"),
        ("docs", "../sage/", f"{BLOB}/sage"),
        (".", "docs/ecc_security.md", "ecc_security.md"),
        (".", "references.bib", f"{BLOB}/references.bib"),
        ("docs", "https://github.com/x", "https://github.com/x"),
        ("docs", "#section", "#section"),
        (".", "docs/index.md", "contents.md"),
        ("docs", "index.md#study", "contents.md#study"),
    ],
)
def test_links_leave_the_book_only_for_the_repository(base, link, expected):
    assert book.relink(f"[a]({link})", base, BLOB) == f"[a]({expected})"


def test_citations_link_to_the_literature():
    text = "[a-2001, Section 2; c-2003] and `[a-2001]`.\n"
    assert book.link_citations(text, {"a-2001"}) == (
        "[[a-2001](literature.md#a-2001), Section 2; c-2003] and `[a-2001]`.\n"
    )


def test_literature_heads_become_anchored_citations():
    text = (
        "# L\n\n- `a-2001`. One,\n  wrapped.\n"
        "- `a-2001`, `b-2002`.\n  Two.\n\nNot `a-2001`. here.\n"
    )
    assert book.literature(text, bibliography.parse(BIB)) == (
        "# L\n\n"
        '- <a id="a-2001"></a>\\[a-2001\\] J. Doe. A. 2001.\n\n  One,\n  wrapped.\n'
        '- <a id="a-2001"></a>\\[a-2001\\] J. Doe. A. 2001.\\\n'
        '  <a id="b-2002"></a>\\[b-2002\\] B. After [[a-2001](#a-2001)].\n\n  Two.\n'
        "\nNot `a-2001`. here.\n"
    )


def test_literature_heads_need_an_annotation():
    with pytest.raises(ValueError):
        book.literature("- `a-2001`.\n- `b-2002`. B.\n", bibliography.parse(BIB))


def test_summary_follows_the_index_without_descriptions():
    index = (
        '---\nokf_version: "0.2"\n---\n\n# Study\n\n'
        "* [Background](background.md) - The checksum.\n\n"
        "# Reference\n\n* [Literature](literature.md) - Entries.\n"
    )
    assert book.summary(index) == (
        "# Summary\n\n[Introduction](README.md)\n[Contents](contents.md)\n\n"
        "# Study\n\n- [Background](background.md)\n\n"
        "# Reference\n\n- [Literature](literature.md)\n"
    )


def test_contents_keep_descriptions_under_one_title():
    index = '---\nokf_version: "0.2"\n---\n\n# Study\n\n* [B](b.md) - The b.\n'
    assert book.contents(index) == "# Contents\n\n## Study\n\n* [B](b.md) - The b.\n"


def test_summary_rejects_what_it_cannot_place():
    with pytest.raises(ValueError):
        book.summary("# Study\n\nprose\n")


def test_assembly_covers_every_chapter(tmp_path):
    root = tmp_path / "repo"
    (root / "docs").mkdir(parents=True)
    (root / "references.bib").write_text(BIB)
    (root / "README.md").write_text("# R\n\n[c](docs/c.md) [a-2001]\n")
    (root / "docs" / "index.md").write_text("# Part\n\n* [C](c.md) - d\n")
    (root / "docs" / "c.md").write_text("---\ntype: Chapter\n---\n\n[s](../src/x.rs)\n")
    (root / "docs" / "literature.md").write_text("- `a-2001`. See [b-2002].\n")
    (root / "docs" / "README.md").write_text("---\ntype: Readme\n---\n\nR\n")
    out = tmp_path / "book"
    book.assemble(root, out, "https://example.org/r", "abc")
    src = out / "src"
    assert sorted(p.name for p in src.iterdir()) == [
        "README.md",
        "SUMMARY.md",
        "c.md",
        "contents.md",
        "literature.md",
    ]
    assert (src / "README.md").read_text() == (
        "# R\n\n[c](c.md) [[a-2001](literature.md#a-2001)]\n"
    )
    assert (src / "literature.md").read_text() == (
        '- <a id="a-2001"></a>\\[a-2001\\] J. Doe. A. 2001.\n\n'
        "  See [[b-2002](literature.md#b-2002)].\n"
    )
    assert (src / "c.md").read_text() == f"[s]({BLOB}/src/x.rs)\n"
    toml = (out / "book.toml").read_text()
    assert 'git-repository-url = "https://example.org/r"' in toml
    assert 'site-url = "/r/"' in toml
    assert 'additional-css = ["book.css"]' in toml
    assert "--content-max-width" in (out / "book.css").read_text()


INDEX = (
    "# Study\n\n* [C](c.md) - d\n\n# Evidence\n\n* [G](g.md) - e\n\n"
    "# Reference\n\n* [L](literature.md) - f\n"
)


def test_entries_close_their_part():
    assert book.with_entries(INDEX, "# Evidence", ["* [R](r.md) - x"]) == (
        "# Study\n\n* [C](c.md) - d\n\n# Evidence\n\n* [G](g.md) - e\n"
        "* [R](r.md) - x\n\n# Reference\n\n* [L](literature.md) - f\n"
    )
    assert book.with_entries(INDEX, "# Evidence", []) == INDEX
    last = book.with_entries(INDEX, "# Reference", ["* [R](r.md) - x"])
    assert last.endswith("* [L](literature.md) - f\n* [R](r.md) - x\n")
    with pytest.raises(ValueError):
        book.with_entries(INDEX, "# Results", ["* [R](r.md) - x"])


META = {
    "name": "m4-native",
    "host": "box.local",
    "started": "2026-10-04T09:00:00Z",
    "cpu": "Apple M4",
    "os": "macOS 26.6.2",
    "arch": "arm64",
    "rustc": "rustc 1.98.1 (48a229cea 2026-09-01)",
    "rustflags": "-C target-cpu=apple-m4",
    "commit": "0123456789abcdef0123456789abcdef01234567",
    "dirty": False,
    "profile": "full",
    "group_fixtures": {
        "binary.127": {"r": str(2**126 + 1), "cofactor": 2, "automorphisms": 2}
    },
}


def repository(root, published=()):
    (root / "docs").mkdir(parents=True)
    (root / "references.bib").write_text(BIB)
    (root / "README.md").write_text("# R\n")
    (root / "docs" / "index.md").write_text(INDEX)
    for name in ("c", "g", "literature"):
        (root / "docs" / f"{name}.md").write_text(f"# {name}\n")
    for name in published:
        run = root / "bench-runs" / name
        n = 1024
        for g, v in zip(
            ["group.add", "h2c", "group.prepare", "group.encode", "group.decode"],
            [10, 100, 20, 30, 40],
        ):
            mode = "mode=throughput" if g == "group.add" else f"mode=batch,n={n}"
            write_bench(
                run / "criterion",
                g,
                f"binary.127/{mode}",
                None,
                0.9 * v * n,
                v * n,
                1.1 * v * n,
                {"Elements": n},
            )
        (run / "meta.json").write_text(json.dumps(META | {"name": name}))
        bench_report.export(run, root / "results" / name)
    return root


def test_without_published_runs_the_index_is_the_summary(tmp_path):
    root = repository(tmp_path / "repo")
    (root / "results").mkdir()
    book.assemble(root, tmp_path / "book", "https://example.org/r", "abc")
    src = tmp_path / "book" / "src"
    assert (src / "SUMMARY.md").read_text() == book.summary(INDEX)
    assert not (src / "results").exists()


def test_each_published_run_is_a_chapter_of_evidence(tmp_path):
    root = repository(tmp_path / "repo", ["x86-native", "m4-native"])
    book.assemble(root, tmp_path / "book", "https://example.org/r", "abc")
    src = tmp_path / "book" / "src"
    summary = (src / "SUMMARY.md").read_text()
    assert (
        "# Evidence\n\n- [G](g.md)\n"
        "- [Benchmark run m4-native](results/m4-native/report.md)\n"
        "- [Benchmark run x86-native](results/x86-native/report.md)\n\n# Reference"
    ) in summary
    contents = (src / "contents.md").read_text()
    assert (
        "* [Benchmark run m4-native](results/m4-native/report.md) - "
        "Apple M4 · rustc 1.98.1 · RUSTFLAGS=-C target-cpu=apple-m4 · "
        "0123456789ab · full profile · 2026-10-04\n"
    ) in contents
    page = src / "results" / "m4-native" / "report.md"
    style, text = page.read_text().split("\n\n", 1)
    assert style.startswith("<style>table.grid{")
    assert text.startswith(
        "# Benchmark run m4-native\n\nRun m4-native started 2026-10-04T09:00:00Z "
        "on Apple M4 (macOS 26.6.2, arm64) and timed commit "
        "0123456789abcdef0123456789abcdef01234567, built by rustc 1.98.1 "
        "(48a229cea 2026-09-01) with RUSTFLAGS -C target-cpu=apple-m4, in the "
        "full profile.\n\n"
    )
    assert "## Families compared" in text and "box.local" not in text
    # the grids shaded, as the report's introduction states
    assert '<table class="grid">' in text and "background:" in text
    light = re.findall(r'<img class="fig-light" src="([^"]+)"', text)
    dark = re.findall(r'<img class="fig-dark" src="([^"]+)"', text)
    assert light and all(f.endswith(".svg") for f in light)
    assert dark == [f.removesuffix(".svg") + "-dark.svg" for f in light]
    assert all((page.parent / f).is_file() for f in light + dark)
    # the twin is the figure with the neutrals swapped, and only those
    for f, d in zip(light, dark):
        svg, twin = ((page.parent / n).read_text() for n in (f, d))
        for old, new in figures.DARK.items():
            svg = svg.replace(old, new)
        assert twin == svg and bench_report.SURFACE not in twin
    assert not list(page.parent.glob("*.png")) + list(page.parent.glob("*.html"))


def test_a_run_name_must_be_a_path_segment(tmp_path):
    root = repository(tmp_path / "repo", ["m4 native"])
    with pytest.raises(ValueError, match="m4 native"):
        book.assemble(root, tmp_path / "book", "https://example.org/r", "abc")
