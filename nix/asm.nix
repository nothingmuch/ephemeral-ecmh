{ inputs, ... }:
{
  perSystem =
    {
      pkgs,
      pythonEnv,
      ccEnv,
      toolchains,
      commonArgs,
      ...
    }:
    let
      # nix run .#asm-report -- [--name NAME] [TARGET/CPU[+FEATURE...]...]
      asm-report = pkgs.writeShellApplication {
        name = "asm-report";
        runtimeInputs = [
          toolchains.asm.cargo
          toolchains.asm.rustc
          pkgs.stdenv.cc
          pkgs.coreutils
          pkgs.gnused
          pkgs.jq
          # llvm-cxxfilt, which demangles Rust's v0 symbols
          pkgs.llvmPackages.llvm
          pythonEnv
        ];
        # ccEnv for the cc wrapper, which links the build script
        runtimeEnv = {
          ASM_REPORT = ../report/asm_report.py;
        }
        // ccEnv;
        text = builtins.readFile ./source-id.sh + builtins.readFile ./asm-report.sh;
      };

      asm-listings = toolchains.asm.mkCargoDerivation (
        commonArgs
        // {
          pname = "ecmh-asm-listings";
          version = inputs.self.shortRev or inputs.self.dirtyShortRev or "source";
          cargoArtifacts = null;
          doInstallCargoArtifacts = false;
          nativeBuildInputs = [ asm-report ];
          ASM_SOURCE_REV = inputs.self.rev or inputs.self.dirtyRev or "unknown";
          ASM_SOURCE_DIRTY = if inputs.self ? rev then "false" else "true";
          buildPhaseCargoCommand = ''
            export ASM_RUNS="$TMPDIR/listings" MPLCONFIGDIR="$TMPDIR/mpl"
            asm-report --name matrix
          '';
          installPhaseCommand = ''
            mkdir -p "$out"
            cp -R "$TMPDIR/listings/matrix/"* "$out/"
          '';
          # Only assembly, JSON, TSV and Markdown; no executables to fix up.
          dontFixup = true;
        }
      );
    in
    {
      packages = { inherit asm-report asm-listings; };

      apps.asm-report = {
        type = "app";
        program = "${asm-report}/bin/asm-report";
        meta.description = "The asm probes' assembly for a matrix of targets and CPUs, summarised";
      };
    };
}
