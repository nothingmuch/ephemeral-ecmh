{ inputs, ... }:
{
  perSystem =
    { system, ... }:
    let
      pkgs = import inputs.nixpkgs {
        inherit system;
        overlays = [ inputs.rust-overlay.overlays.default ];
      };

      commonExtensions = [
        "clippy"
        "rust-analyzer"
        "rust-src"
      ];

      # stable builds the package, the bench executables and the dev shells:
      # everything the code needs (core::arch pclmulqdq / pmull intrinsics)
      # is stable. The checks build with nightly for its parallel front end
      # (nix/checks.nix); miri is for the nightly dev shell.
      rustToolchains = {
        stable = pkgs.rust-bin.stable.latest.default.override {
          extensions = commonExtensions ++ [ "llvm-tools" ];
        };
        # stable plus the standard libraries of the targets whose assembly
        # asm-report reads: compiling for them needs no linker, so no cross
        # toolchain. A separate toolchain, so the checks' stays as it is.
        asm = pkgs.rust-bin.stable.latest.default.override {
          extensions = commonExtensions;
          targets = [
            "aarch64-apple-darwin"
            "aarch64-linux-android"
            "x86_64-unknown-linux-gnu"
          ];
        };
        nightly = pkgs.rust-bin.selectLatestNightlyWith (
          t:
          t.default.override {
            extensions = commonExtensions ++ [
              "llvm-tools"
              "miri"
            ];
          }
        );
      };

      mkCraneLib = _: rust: (inputs.crane.mkLib pkgs).overrideToolchain rust;
      toolchains = builtins.mapAttrs mkCraneLib rustToolchains;

      # for the apps that run cargo with the cc wrapper outside a nix build,
      # where it doesn't know where libiconv is, which Rust's std links on
      # Darwin, and where rustc would target an older macOS (11.0 on
      # aarch64) than the stdenv's, which the wrapper and the dev shells use
      ccEnv = pkgs.lib.optionalAttrs pkgs.stdenv.hostPlatform.isDarwin {
        LIBRARY_PATH = "${pkgs.libiconv}/lib";
        MACOSX_DEPLOYMENT_TARGET = pkgs.stdenv.hostPlatform.darwinMinVersion;
      };
    in
    {
      _module.args = {
        inherit pkgs toolchains ccEnv;
      };
    };
}
