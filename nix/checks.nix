{ inputs, lib, ... }:
{
  perSystem =
    {
      config,
      pkgs,
      hasCargo,
      commonArgs,
      toolchains,
      ...
    }:
    let
      rev = inputs.self.shortRev or "dirty";
      # Every Rust check builds with one toolchain, nightly, in one profile,
      # ci (Cargo.toml): dev's debug assertions and overflow checks, but
      # optimized, as the curvegen tests take 20x longer at opt-level 0.
      # Release code is what the benches build and time.
      craneLib = toolchains.nightly;
      ciArgs = commonArgs // {
        CARGO_PROFILE = "ci";
        # nightly's parallel front end: the lib's test binary, proptests
        # and all, is most of a check's compile
        RUSTFLAGS = "-Zthreads=8";
      };
      cargoArtifacts = craneLib.buildDepsOnly ciArgs;
      checkArgs = ciArgs // {
        inherit cargoArtifacts;
        version = rev;
        dontFixup = true;
        doInstallCargoArtifacts = false;
      };
      src = commonArgs.src;

      rustChecks = {
        # every test, the `pari` feature's included: one compile of the
        # crate rather than one per feature set (the feature only adds code)
        tests = craneLib.cargoNextest (checkArgs // { cargoNextestExtraArgs = "--no-tests=warn"; });

        clippy = craneLib.cargoClippy (
          checkArgs
          // {
            cargoClippyExtraArgs = "--all-targets --all-features -- -D warnings";
          }
        );

        doc = craneLib.cargoDoc (
          ciArgs
          // {
            inherit cargoArtifacts;
            cargoDocExtraArgs = "--no-deps --all-features";
            RUSTDOCFLAGS = "-D warnings";
          }
        );

        cargo-sort =
          pkgs.runCommand "cargo-sort-${rev}"
            {
              inherit src;
              nativeBuildInputs = [ pkgs.cargo-sort ];
            }
            ''
              cargo-sort --check --workspace "$src"
              mkdir -p $out
            '';
      };

      checks = lib.optionalAttrs hasCargo rustChecks;
    in
    {
      checks =
        checks
        // lib.optionalAttrs hasCargo {
          quick = pkgs.symlinkJoin {
            name = "quick-checks-${rev}";
            paths = with checks; [
              tests
              clippy
            ];
          };
          lint = pkgs.symlinkJoin {
            name = "lint-checks-${rev}";
            paths = with checks; [
              cargo-sort
              clippy
              doc
            ];
          };
          # every check: what each commit must pass (a link farm, as
          # bench-report's output is a file)
          all = pkgs.linkFarm "all-checks-${rev}" (
            {
              inherit (checks)
                tests
                clippy
                doc
                cargo-sort
                ;
            }
            // {
              inherit (config.checks) treefmt site;
            }
          );
        };
    };
}
