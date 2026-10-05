{ inputs, lib, ... }:
let
  # Expose the Rust build outputs when Cargo.toml is present.
  hasCargo = builtins.pathExists (inputs.self + "/Cargo.toml");
in
{
  perSystem =
    { toolchains, ... }:
    let
      craneLib = toolchains.stable;
      # The Cargo sources, and the recorded test fixtures that are not Rust.
      src = lib.fileset.toSource {
        root = ../.;
        fileset = lib.fileset.unions [
          (craneLib.fileset.commonCargoSources ../.)
          ../tests/fixtures
        ];
      };

      commonArgs = {
        inherit src;
        strictDeps = true;
      };

      cargoArtifactsRelease = craneLib.buildDepsOnly commonArgs;
    in
    {
      _module.args = {
        inherit
          hasCargo
          commonArgs
          cargoArtifactsRelease
          ;
      };

      packages = lib.optionalAttrs hasCargo {
        default = craneLib.buildPackage (commonArgs // { cargoArtifacts = cargoArtifactsRelease; });
      };
    };
}
