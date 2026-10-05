"""Check the relative links of a rendered mdbook.

Every relative href or src in an HTML page of the book must name a file
inside the book, and a fragment must name an id in the target page. Links with
a scheme are not followed, nor absolute paths, which depend on where the
book is served (mdbook's 404 page uses them, under `site-url`).

    python linkcheck.py BOOK
"""

import sys
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit


class _Page(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links: list[str] = []
        self.ids: set[str] = set()

    def handle_starttag(self, tag, attrs):
        for name, value in attrs:
            if name == "id" and value:
                self.ids.add(value)
            elif name in ("href", "src") and value:
                self.links.append(value)


def _parse(path: Path) -> _Page:
    page = _Page()
    page.feed(path.read_text(errors="replace"))
    return page


def broken(book: Path) -> list[str]:
    book = book.resolve()
    pages = {p: _parse(p) for p in book.rglob("*.html")}
    out = []
    for path, page in sorted(pages.items()):
        for link in page.links:
            url = urlsplit(link)
            if url.scheme or url.netloc or link.startswith("/"):
                continue
            target = (path.parent / unquote(url.path)).resolve() if url.path else path
            if target.is_dir():
                # a directory link is served as its index page
                target = target / "index.html"
            where = f"{path.relative_to(book)}: {link}"
            if not target.is_relative_to(book) or not target.is_file():
                out.append(f"{where}: no such file in the book")
            elif url.fragment and target.suffix == ".html":
                ids = pages[target].ids if target in pages else _parse(target).ids
                if unquote(url.fragment) not in ids:
                    out.append(f"{where}: no such anchor")
    return out


def main(argv=None) -> None:
    (book,) = argv if argv is not None else sys.argv[1:]
    errors = broken(Path(book))
    for e in errors:
        print(e, file=sys.stderr)
    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()
