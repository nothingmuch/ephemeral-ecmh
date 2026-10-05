{ inputs, lib, ... }:
{
  perSystem =
    { pkgs, pythonEnv, ... }:
    let
      repository = "https://github.com/nothingmuch/ephemeral-ecmh";
      # Links out of the book name the revision the site is built from; a
      # build of a modified tree names the commit it modifies.
      rev = inputs.self.rev or (lib.removeSuffix "-dirty" (inputs.self.dirtyRev or "main"));
      src = lib.fileset.toSource {
        root = ../.;
        fileset = lib.fileset.unions [
          ../README.md
          ../docs
          ../references.bib
          ../report/bibcheck.py
          ../report/bibliography.py
          ../report/book.py
          ../report/linkcheck.py
        ];
      };
      # nix build .#site: the README and docs/ as an mdbook, for GitHub Pages,
      # whose citations are checked against the bibliography and whose relative
      # links are checked against the rendered pages
      site =
        pkgs.runCommand "ecmh-site"
          {
            nativeBuildInputs = [
              pkgs.mdbook
              pkgs.mdbook-katex
              pythonEnv
            ];
          }
          ''
            python ${src}/report/bibcheck.py ${src}
            python ${src}/report/book.py ${src} book --repository ${repository} --rev ${rev}
            mdbook build book -d $out
            python ${src}/report/linkcheck.py $out
          '';
    in
    {
      packages.site = site;
      checks.site = site;
    };
}
