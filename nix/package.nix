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
      # The Cargo sources.
      src = lib.fileset.toSource {
        root = ../.;
        fileset = lib.fileset.unions [
          (craneLib.fileset.commonCargoSources ../.)
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
