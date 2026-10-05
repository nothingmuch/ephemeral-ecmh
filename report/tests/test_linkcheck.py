import linkcheck


def test_missing_files_and_anchors_are_reported(tmp_path):
    (tmp_path / "a.html").write_text(
        '<h2 id="x">X</h2><a href="b.html#y">b</a><a href="b.html#z">b</a>'
        '<a href="c.html">c</a><a href="#x">self</a><a href="https://e.org/#q">e</a>'
    )
    (tmp_path / "b.html").write_text('<p id="y">y</p>')
    assert linkcheck.broken(tmp_path) == [
        "a.html: b.html#z: no such anchor",
        "a.html: c.html: no such file in the book",
    ]


def test_a_file_outside_the_book_is_reported(tmp_path):
    book = tmp_path / "book"
    book.mkdir()
    (tmp_path / "outside.html").write_text("<p>x</p>")
    (book / "a.html").write_text('<a href="../outside.html">o</a>')
    assert linkcheck.broken(book) == [
        "a.html: ../outside.html: no such file in the book"
    ]


def test_a_directory_link_needs_an_index_page(tmp_path):
    (tmp_path / "runs" / "a").mkdir(parents=True)
    (tmp_path / "runs" / "b").mkdir()
    (tmp_path / "runs" / "a" / "index.html").write_text("<p>a</p>")
    (tmp_path / "x.html").write_text('<a href="runs/a/">a</a><a href="runs/b/">b</a>')
    assert linkcheck.broken(tmp_path) == ["x.html: runs/b/: no such file in the book"]
