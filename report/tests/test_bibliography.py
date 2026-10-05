import bibliography
import pytest

BIB = r"""
% a comment, with an @ sign
@inproceedings{blf-2008,
  author    = {Bernstein, Daniel J. and Lange, Tanja and Rezaeian Farashahi, Reza},
  title     = {Binary {Edwards} Curves},
  booktitle = {Cryptographic Hardware and Embedded Systems -- {CHES} 2008},
  series    = {LNCS},
  volume    = {5154},
  publisher = {Springer},
  year      = {2008},
  pages     = {244--265},
  doi       = {10.1007/978-3-540-85053-3_16},
  eprint    = {https://eprint.iacr.org/2008/171},
  file      = {pdfs/blf-2008.pdf},
}

@article{phan-wagner-2006,
  author  = {Phan, Raphael C.-W. and Wagner, David},
  title   = {Security Considerations for Incremental Hash Functions},
  journal = {Computers \& Security},
  volume  = {25},
  number  = {2},
  year    = {2006},
  pages   = {131--136},
  url     = {https://people.eecs.berkeley.edu/~daw/papers/inchash-cs06.pdf},
}

@incollection{hess-2005,
  author    = {Hess, Florian},
  title     = {Weil Descent Attacks},
  booktitle = {Advances in Elliptic Curve Cryptography},
  editor    = {Blake, Ian F. and Seroussi, Gadiel},
  series    = {London Mathematical Society Lecture Note Series},
  volume    = {317},
  publisher = {Cambridge University Press},
  year      = {2005},
  pages     = {151--180},
}

@techreport{fips-180-4,
  author      = {{National Institute of Standards and Technology}},
  title       = {Secure Hash Standard ({SHS})},
  number      = {FIPS PUB 180-4},
  institution = {NIST},
  year        = {2015},
}

@misc{rfc9496,
  author       = {Henry de Valence and Jack Grigg},
  title        = {The ristretto255 and decaf448 Groups},
  howpublished = {{RFC} 9496},
  year         = {2023},
  month        = dec,
  % a comment inside an entry
  doi          = {10.17487/RFC9496},
  url          = {https://www.rfc-editor.org/rfc/rfc9496},
}

@misc{safecurves,
  author       = {Bernstein, Daniel J. and Lange, Tanja},
  title        = {{SafeCurves}: choosing safe curves},
  howpublished = {\url{https://safecurves.cr.yp.to}},
  note         = {Fields: \url{https://safecurves.cr.yp.to/field.html}},
}

@unpublished{mestre-2000,
  author = {Mestre, Jean-Fran{\c{c}}ois},
  title  = {Lettre adress{\'e}e {\`a} {Gaudry}},
  note   = {Letter, December 2000; cited through \cite{blf-2008}},
  year   = {2000},
}
"""


def render(key: str) -> str:
    return bibliography.render(bibliography.parse(BIB)[key], lambda k: f"#{k}")


def test_parse_keeps_order_and_raw_values():
    entries = bibliography.parse(BIB)
    assert list(entries)[:2] == ["blf-2008", "phan-wagner-2006"]
    kind, fields = entries["rfc9496"]
    assert kind == "misc"
    assert fields["month"] == "dec"
    assert fields["howpublished"] == "{RFC} 9496"


def test_parse_rejects_a_duplicate_key():
    with pytest.raises(ValueError, match="duplicate"):
        bibliography.parse("@misc{a,\n title={A},\n}\n@misc{a,\n title={B},\n}\n")


@pytest.mark.parametrize(
    "key, expected",
    [
        (
            "blf-2008",
            (
                "D. J. Bernstein, T. Lange, R. Rezaeian Farashahi. Binary Edwards "
                "Curves. Cryptographic Hardware and Embedded Systems – CHES 2008, "
                "LNCS 5154, pp. 244–265. <https://doi.org/10.1007/978-3-540-85053-3_16>"
            ),
        ),
        (
            "phan-wagner-2006",
            (
                "R. C.-W. Phan, D. Wagner. Security Considerations for Incremental "
                "Hash Functions. Computers & Security 25(2), pp. 131–136, 2006. "
                "<https://people.eecs.berkeley.edu/~daw/papers/inchash-cs06.pdf>"
            ),
        ),
        (
            "hess-2005",
            (
                "F. Hess. Weil Descent Attacks. In: Advances in Elliptic Curve "
                "Cryptography, ed. I. F. Blake, G. Seroussi, London Mathematical "
                "Society Lecture Note Series 317, Cambridge University Press, "
                "pp. 151–180, 2005."
            ),
        ),
        (
            "fips-180-4",
            (
                "National Institute of Standards and Technology. Secure Hash "
                "Standard (SHS). NIST, FIPS PUB 180-4, 2015."
            ),
        ),
        (
            "rfc9496",
            (
                "H. de Valence, J. Grigg. The ristretto255 and decaf448 Groups. "
                "RFC 9496, December 2023. <https://doi.org/10.17487/RFC9496>"
            ),
        ),
        (
            "safecurves",
            (
                "D. J. Bernstein, T. Lange. SafeCurves: choosing safe curves. "
                "<https://safecurves.cr.yp.to>. "
                "Fields: <https://safecurves.cr.yp.to/field.html>."
            ),
        ),
        (
            "mestre-2000",
            (
                "J.-F. Mestre. Lettre adressée à Gaudry. Letter, December 2000; "
                "cited through [[blf-2008](#blf-2008)]."
            ),
        ),
    ],
)
def test_entries_render_by_type(key, expected):
    assert render(key) == expected


@pytest.mark.parametrize(
    "author, expected",
    [
        ("Bernstein, Daniel J.", "D. J. Bernstein"),
        ("de Valence, Henry", "H. de Valence"),
        ("Henry de Valence", "H. de Valence"),
        ('Ulrich Hab{\\"o}ck', "U. Haböck"),
        ("Aumasson, Jean-Philippe", "J.-P. Aumasson"),
        ("Brand{\\~a}o, Lu{\\'\\i}s T. A. N.", "L. T. A. N. Brandão"),
        ("{RISC Zero Team}", "RISC Zero Team"),
        ("others", "et al."),
    ],
)
def test_names_are_initials_and_surname(author, expected):
    assert bibliography.name(author) == expected


def test_text_is_escaped_for_markdown_and_math_unwrapped():
    entry = ("misc", {"title": r"{BEC\_Small} over $\mathbb{F}_{2^m}$ [draft]"})
    assert bibliography.render(entry, str) == r"BEC\_Small over F\_{2^m} \[draft\]."


def test_unknown_latex_is_refused():
    with pytest.raises(ValueError, match="unsupported LaTeX"):
        bibliography.render(("misc", {"title": r"\emph{x}"}), str)


def test_citations_are_key_groups_outside_code_and_links():
    text = (
        "See [a-2001, Section 4;\nb-2002] and `[c-2003]`, [x](y.md), [TODO]\n"
        "```\n[efd]\n```\nbut [efd].\n"
    )
    assert list(bibliography.citations(text)) == [
        ["a-2001", "Section 4", "b-2002"],
        ["efd"],
    ]


def test_linked_citations_keep_locators_and_separators():
    text = "See [a-2001, Section 4;\nb-2002], `[a-2001]`, [x](y.md) and [TODO].\n"
    href = {"a-2001": "l.md#a-2001", "b-2002": "l.md#b-2002"}.get
    assert bibliography.link_citations(text, href) == (
        "See [[a-2001](l.md#a-2001), Section 4;\n[b-2002](l.md#b-2002)], "
        "`[a-2001]`, [x](y.md) and [TODO].\n"
    )
