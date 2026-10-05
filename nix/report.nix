{ inputs, ... }:
{
  perSystem =
    {
      pkgs,
      pythonEnv,
      pariEnv,
      ccEnv,
      toolchains,
      commonArgs,
      ...
    }:
    let
      src = ../report;

      # Every criterion suite and examples/curvegen_times, in the bench
      # profile with the pari feature, for one named target CPU: crrl's
      # GF(2^127), gf2_109 and gf2_122 use carry-less multiplication only
      # when the target enables it (aes on aarch64, pclmulqdq on x86_64),
      # so the target is part of what a run times. `generic` is the
      # architecture's baseline, without it.
      # The first of each architecture's is the default.
      targets =
        {
          aarch64 = [
            [
              "apple-m4"
              "-C target-cpu=apple-m4"
            ]
            [
              "generic"
              "-C target-cpu=generic"
            ]
          ];
          x86_64 = [
            [
              "x86-64-v3"
              "-C target-cpu=x86-64-v3 -C target-feature=+pclmulqdq"
            ]
            [
              "generic"
              "-C target-cpu=x86-64"
            ]
          ];
        }
        .${pkgs.stdenv.hostPlatform.parsed.cpu.name};
      craneLib = toolchains.stable;
      benchBins =
        target: rustflags:
        let
          args =
            commonArgs
            // pariEnv
            // ccEnv
            // {
              pname = "ecmh-bench-bins-${target}";
              version = inputs.self.shortRev or inputs.self.dirtyShortRev or "unknown";
              RUSTFLAGS = rustflags;
              CARGO_PROFILE = "bench";
            };
        in
        craneLib.mkCargoDerivation (
          args
          // {
            cargoArtifacts = craneLib.buildDepsOnly (args // { doCheck = false; });
            nativeBuildInputs = [ pkgs.jq ];
            doInstallCargoArtifacts = false;
            buildPhaseCargoCommand = ''
              cargo bench $cargoExtraArgs --no-run --message-format=json-render-diagnostics >bench.json
            '';
            installPhaseCommand = ''
              mkdir -p $out/bin
              jq -r 'select(.reason == "compiler-artifact" and .executable != null
                  and (.target.kind == ["bench"] or .target.kind == ["example"]))
                | "\(.target.name) \(.executable)"' bench.json example.json |
                while read -r name exe; do install -m755 "$exe" "$out/bin/$name"; done
              # what bench-run records of the build: the features the flags
              # actually enable, e.g. pclmulqdq, avx2, aes
              echo ${target} >$out/target
              echo "$RUSTFLAGS" >$out/rustflags
              rustc --print cfg $RUSTFLAGS | sed -n 's/^target_feature="\(.*\)"$/\1/p' |
                paste -sd, - >$out/target-features
              rustc --version >$out/rustc
              echo ${inputs.self.rev or inputs.self.dirtyRev or "unknown"} >$out/source
            '';
          }
        );
      bench-bins-by-target = builtins.listToAttrs (
        map (t: {
          name = "bench-bins-${builtins.elemAt t 0}";
          value = benchBins (builtins.elemAt t 0) (builtins.elemAt t 1);
        }) targets
      );
      bench-bins = bench-bins-by-target."bench-bins-${builtins.elemAt (builtins.head targets) 0}";

      # nix run .#bench-report -- target/criterion out/
      bench-report = pkgs.writeShellApplication {
        name = "bench-report";
        runtimeInputs = [ pythonEnv ];
        text = ''exec python ${src}/bench_report.py "$@"'';
      };

      # nix run .#bench-run -- [--profile quick|full] [SUITE...] [-- CRITERION ARGS]
      bench-run = pkgs.writeShellApplication {
        name = "bench-run";
        runtimeInputs = [
          pkgs.jq
          pkgs.gawk
          # GNU date, for the milliseconds of the run's UUIDv7
          pkgs.coreutils
          bench-report
        ];
        runtimeEnv = pariEnv // {
          BENCH_BINS_BUILT = "${bench-bins}";
        };
        text = builtins.readFile ./bench-run.sh;
      };
    in
    {
      packages = {
        inherit bench-report bench-run bench-bins;
      }
      // bench-bins-by-target;

      apps.bench-run = {
        type = "app";
        program = "${bench-run}/bin/bench-run";
        meta.description = "Run the criterion suites serially on this machine, record it, and report";
      };

      apps.bench-report = {
        type = "app";
        program = "${bench-report}/bin/bench-report";
        meta.description = "Criterion results to facet-grid plots and a markdown/HTML report";
      };

      checks.bench-report =
        pkgs.runCommand "bench-report-tests"
          {
            nativeBuildInputs = [ pythonEnv ];
          }
          ''
            # matplotlib wants a writable config and font cache
            export HOME=$TMPDIR MPLCONFIGDIR=$TMPDIR/mpl
            cd ${src}
            python -m pytest -q -p no:cacheprovider tests
            touch $out
          '';
    };
}
