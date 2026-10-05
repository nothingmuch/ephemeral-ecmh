{
  perSystem =
    { pkgs, ... }:
    let
      # seadata: SEA's modular polynomials (without them PARI computes its
      # own, far slower); galdata as sage ships it
      pariData = pkgs.symlinkJoin {
        name = "pari-data";
        paths = [
          pkgs.pari-galdata
          pkgs.pari-seadata-small
        ];
      };
    in
    {
      # Sage is used for curve parameter generation / verification (point
      # counting, GHS magic numbers, embedding degree) and as the executable
      # specification of the group law. PARI's gp is exposed separately for
      # quick one-off queries.
      packages.sage = pkgs.sage;
      packages.pari = pkgs.pari;

      # the crate's `pari` feature: build.rs links libpari from PARI_PREFIX,
      # and src/curvegen/pari.rs falls back to GP_DATA_DIR as it was at build time
      _module.args.pariEnv = {
        PARI_PREFIX = "${pkgs.pari}";
        GP_DATA_DIR = "${pariData}/share/pari";
      };
    };
}
