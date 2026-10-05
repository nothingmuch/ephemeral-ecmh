{ inputs, ... }:
{
  imports = [ inputs.treefmt-nix.flakeModule ];

  perSystem = {
    treefmt = {
      projectRootFile = "flake.nix";

      programs.nixfmt.enable = true;
      programs.rustfmt.enable = true;
      programs.ruff-format.enable = true;
      programs.ruff-check.enable = true;
      programs.shellcheck.enable = true;
      programs.shfmt.enable = true;
      programs.just.enable = true;
      programs.taplo.enable = true;

      # .sage files are Python-with-preparser syntax; ruff cannot parse them.
      settings.global.excludes = [
        ".envrc"
        "*.sage"
      ];
    };
  };
}
