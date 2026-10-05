{
  perSystem =
    { pkgs, lib, ... }:
    let
      # The Go RIBLT benchmarks of yangl1996/riblt#3: upstream, which
      # checksums with the XOR of SipHash values, and David Buchanan's
      # fork, which replaces it by an ECMH over ristretto255 and then the
      # 64-bit mapping PRNG by ChaCha8. The fork's history contains
      # upstream's, so one source serves all three revisions.
      revisions = {
        upstream = {
          rev = "4afa6bc06cb2237d9ea273a51d97a7e05b3f573b";
          hash = "sha256-OSREJEtyI72Gd7voMchQ1Btm35bbALRdDmNpvfbn6HQ=";
          vendorHash = "sha256-78AgfdvZXI+vnAPl4jyrFwoq5qBR7AbU3ruv62ZA6w4=";
        };
        ristretto255 = {
          rev = "8c7dd499a8fdf2237e08240d17af5deddb1d428f";
          hash = "sha256-HwzwrnttJTQTRYJ/hU/RL3qvnv5gbL2ZNG+qOjN1dDs=";
          vendorHash = "sha256-Au/XSn37BHtrNWriqRgPT1w0+quKIolJbKE1tZCEV9I=";
        };
        ristretto255-chacha8 = {
          rev = "11521a4c73ebbbfe80050ec99a4f976eb9ce7147";
          hash = "sha256-dCezqydWNE5dH1nszNRDT2K7Y3nYLXQlo9qx5WOD/O0=";
          vendorHash = "sha256-Au/XSn37BHtrNWriqRgPT1w0+quKIolJbKE1tZCEV9I=";
        };
      };

      # The package's test binary, built as `go test -c` would.
      testBinary =
        name: r:
        pkgs.buildGoModule {
          pname = "riblt-go-${name}";
          version = builtins.substring 0 7 r.rev;
          src = pkgs.fetchFromGitHub {
            owner = "DavidBuchanan314";
            repo = "riblt-ecmh";
            inherit (r) rev hash;
          };
          inherit (r) vendorHash;
          doCheck = false;
          buildPhase = ''
            runHook preBuild
            go test -c -o riblt.test .
            runHook postBuild
          '';
          installPhase = ''
            runHook preInstall
            install -Dm755 riblt.test $out/bin/riblt-${name}.test
            runHook postInstall
          '';
        };

      binaries = lib.mapAttrs testBinary revisions;

      # nix run .#riblt-go-bench -- [go test flags], e.g. -test.bench EncodeAndDecode
      #
      # Runs every revision's benchmarks in turn, with the flags of the
      # issue's `go test -bench .` unless others are given.
      bench = pkgs.writeShellApplication {
        name = "riblt-go-bench";
        text = ''
          if [ "$#" -eq 0 ]; then set -- -test.run '^$' -test.bench .; fi
          ${lib.concatMapStrings (name: ''
            echo "# ${name} ${revisions.${name}.rev}"
            ${binaries.${name}}/bin/riblt-${name}.test "$@"
          '') (lib.attrNames revisions)}
        '';
      };
    in
    {
      packages = lib.mapAttrs' (n: v: lib.nameValuePair "riblt-go-${n}" v) binaries // {
        riblt-go-bench = bench;
      };

      apps.riblt-go-bench = {
        type = "app";
        program = "${bench}/bin/riblt-go-bench";
        meta.description = "Run the Go RIBLT benchmarks of yangl1996/riblt#3: upstream and the ristretto255 fork";
      };
    };
}
