{
  perSystem =
    { pkgs, ... }:
    let
      # Lightweight reference models and test-vector generators that do not
      # need a full Sage install, and the benchmark report (report/).
      # Every dependency is a nixpkgs package; one that is not should come
      # through uv2nix, not pip or a venv, so that checks and CI resolve the
      # same environment. Scripts run under this output, not a system python.
      pythonEnv = pkgs.python3.withPackages (ps: [
        ps.galois
        ps.hypothesis
        ps.matplotlib
        ps.pandas
        ps.pytest
        ps.seaborn
      ]);
    in
    {
      _module.args = {
        inherit pythonEnv;
      };

      packages.python = pythonEnv;
    };
}
