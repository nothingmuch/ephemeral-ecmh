import bibliography
import book
import pytest

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
