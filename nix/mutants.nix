{
  perSystem =
    {
      pkgs,
      ccEnv,
      toolchains,
      ...
    }:
    let
      # nix run .#mutants -- [CARGO MUTANTS ARGS], e.g. -f src/field/gf2_109/mod.rs
      #
      # An app, not a check: a mutant costs a build and a test run, so the
      # whole crate takes hours. Tests build in the checks' profile.
      mutants = pkgs.writeShellApplication {
        name = "mutants";
        runtimeInputs = [
          toolchains.stable.cargo
          toolchains.stable.rustc
          pkgs.cargo-mutants
          pkgs.cargo-nextest
          pkgs.stdenv.cc
        ];
        runtimeEnv = ccEnv;
        text = ''exec cargo mutants --profile ci "$@"'';
      };
    in
    {
      packages.mutants = mutants;

      apps.mutants = {
        type = "app";
        program = "${mutants}/bin/mutants";
        meta.description = "Mutation-test the crate, or the files given, with cargo-mutants";
      };
    };
}
