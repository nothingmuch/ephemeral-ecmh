{
  perSystem =
    {
      config,
      pkgs,
      pythonEnv,
      pariEnv,
      toolchains,
      ...
    }:
    let
      rustTools = with pkgs; [
        cargo-llvm-cov
        cargo-nextest
        cargo-show-asm
        cargo-sort
        hyperfine
        just
        rust-analyzer
      ];

      mkDevShell =
        craneLib: extra:
        craneLib.devShell (
          pariEnv
          // {
            packages = rustTools ++ [ config.treefmt.build.wrapper ] ++ extra;
          }
        );

      mathTools = [
        pythonEnv
        config.packages.sage
        config.packages.pari
      ];
    in
    {
      devShells = {
        # everything: rust + python reference models + sage/pari
        default = mkDevShell toolchains.stable mathTools;
        # rust only (no multi-GiB sage closure), e.g. for CI
        rust = mkDevShell toolchains.stable [ ];
        nightly = mkDevShell toolchains.nightly [ ];
        # rust with cross standard libraries, for reading other targets'
        # assembly (cargo asm --target ...; nix run .#asm-report)
        asm = mkDevShell toolchains.asm [ ];
        # math only: parameter generation and reference models
        math = pkgs.mkShell {
          packages = mathTools ++ [ config.treefmt.build.wrapper ];
        };
      };
    };
}
