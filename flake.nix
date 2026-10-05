{
  description = "Ephemeral elliptic-curve multiset hash over a binary curve (~2^63 security)";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    flake-parts.url = "github:hercules-ci/flake-parts";
    rust-overlay = {
      url = "github:oxalica/rust-overlay";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    crane.url = "github:ipetkov/crane";
    treefmt-nix = {
      url = "github:numtide/treefmt-nix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };

  outputs =
    inputs:
    inputs.flake-parts.lib.mkFlake { inherit inputs; } {
      systems = [
        "aarch64-darwin"
        "x86_64-darwin"
        "aarch64-linux"
        "x86_64-linux"
      ];

      imports = [
        ./nix/toolchain.nix
        ./nix/python.nix
        ./nix/report.nix
        ./nix/mutants.nix
        ./nix/sage.nix
        ./nix/package.nix
        ./nix/site.nix
        ./nix/checks.nix
        ./nix/devshell.nix
        ./nix/treefmt.nix
        ./nix/validate-commits.nix
      ];
    };
}
