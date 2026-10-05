import bibcheck

KEYS = {"a-2001", "b-2002", "efd"}


def test_citations_are_key_groups_outside_code_and_links():
    text = (
        "---\ntags: [a-2001, x]\n---\n"
        "See [a-2001] and [a-2001, Section 4;\nb-2002, Section 5]\n"
        "and [b-2002, c-2003], [d-2004] and [e-2005, a-2001].\n"
        "Not [Adversary](x.md#a) nor [TODO] nor [see below]: `[a-2001]` nor\n"
        "```\n[efd]\n```\nbut [efd].\n"
    )
    assert bibcheck.citations(text, KEYS) == (
        ["a-2001", "a-2001", "b-2002", "b-2002", "a-2001", "efd"],
        ["c-2003", "d-2004", "e-2005"],
    )


def test_annotations_list_keys():
    assert bibcheck.annotations(
        "- `a-2001`. One.\n- `b-2002`, `efd`. Two,\n  wrapped.\n- `x`, not a bullet\n"
    ) == ["a-2001", "b-2002", "efd"]


def test_problems_name_every_gap(tmp_path):
    (tmp_path / "references.bib").write_text(
        "@misc{a-2001,\n  title = {A},\n}\n@article{b-2002,\n  title = {B},\n}\n"
    )
    docs = tmp_path / "docs"
    docs.mkdir()
    (tmp_path / "README.md").write_text("[a-2001] and [a-2001, c-2003].\n")
    (docs / "x.md").write_text("---\ntype: Chapter\n---\n\nbody\n")
    (docs / "literature.md").write_text("- `a-2001`. A.\n- `a-2001`. A.\n")
    assert bibcheck.problems(tmp_path) == [
        "README.md: [c-2003] is not in references.bib",
        "docs/literature.md: missing b-2002",
        "docs/literature.md: a-2001 listed more than once",
    ]
