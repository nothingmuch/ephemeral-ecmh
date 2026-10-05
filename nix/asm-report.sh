# shellcheck shell=bash
# The compiler's output for every configuration in a matrix of targets and
# CPUs, and report/asm_report.py's summary of the probes in it
# (src/asm_probes.rs):
#
#   nix run .#asm-report -- [--name NAME] [TARGET/CPU[+FEATURE...]...]
#
# Run from the repository root. A configuration is a target triple and an
# LLVM CPU, "default" for the target's baseline, with target features to
# add, e.g. aarch64-linux-android/cortex-x4+aes; they default to the
# matrix below. The run goes to asm-runs/NAME (NAME defaults to
# host-date), which must not exist yet:
#
#   meta.json             the toolchain, Cargo.lock's hash, the commit, and
#                         each configuration's flags and target features
#   asm/CONFIG.s          rustc's whole listing for each configuration
#   summary.md, probes.tsv, CONFIG/PROBE.s   asm_report.py's
#
# Builds are release builds of the library with the asm-probes feature,
# RUSTFLAGS for every crate, as bench-bins sets its target CPU, each
# configuration in its own target directory under target/asm, so a second
# run rebuilds only what changed. Nothing runs, so nothing here is timed.

# LLVM's CPU models: Apple's A14 to A17 cores (M1 and M2 share A14's and
# A15's, M3 A17's), M4 and M5; Android's little, big and prime Cortex
# cores and Snapdragon's Oryon, with +aes where the model leaves out the
# crypto extension that the phones' chips have (it is optional from Armv8
# on, and the GF(2^n) fields' PMULL needs it); and x86-64's distribution
# baselines (v3 has no PCLMULQDQ) and recent AMD and Intel cores. And the
# two bench-bins targets the models above miss: aarch64 generic (no crypto
# extension) and x86-64-v3+pclmulqdq.
default_configs=(
  aarch64-apple-darwin/default
  aarch64-apple-darwin/generic
  aarch64-apple-darwin/apple-a14
  aarch64-apple-darwin/apple-a15
  aarch64-apple-darwin/apple-a16
  aarch64-apple-darwin/apple-a17
  aarch64-apple-darwin/apple-m4
  aarch64-apple-darwin/apple-m5
  aarch64-linux-android/default
  aarch64-linux-android/cortex-a55+aes
  aarch64-linux-android/cortex-a78+aes
  aarch64-linux-android/cortex-a720+aes
  aarch64-linux-android/cortex-x4
  aarch64-linux-android/cortex-x4+aes
  aarch64-linux-android/cortex-x925+aes
  aarch64-linux-android/oryon-1
  x86_64-unknown-linux-gnu/default
  x86_64-unknown-linux-gnu/x86-64-v3
  x86_64-unknown-linux-gnu/x86-64-v3+pclmulqdq
  x86_64-unknown-linux-gnu/x86-64-v4
  x86_64-unknown-linux-gnu/znver4
  x86_64-unknown-linux-gnu/znver5
  x86_64-unknown-linux-gnu/sapphirerapids
  x86_64-unknown-linux-gnu/arrowlake
)

name=""
configs=()
while (($#)); do
  case $1 in
  --name)
    name=$2
    shift 2
    ;;
  -*)
    echo "asm-report: unknown option $1" >&2
    exit 2
    ;;
  *)
    configs+=("$1")
    shift
    ;;
  esac
done
((${#configs[@]})) || configs=("${default_configs[@]}")

[[ -f Cargo.toml && -f src/asm_probes.rs ]] || {
  echo "asm-report: run from the repository root" >&2
  exit 2
}

# Nix supplies the identity of its immutable source input. Runtime observations
# belong only to interactive runs; they would make a build output nondeterministic.
if [[ -n ${ASM_SOURCE_REV:-} ]]; then
  source_id() { printf '%s %s\n' "$ASM_SOURCE_REV" "${ASM_SOURCE_DIRTY:-false}"; }
  host=""
  name=${name:-matrix}
else
  host=$(hostname -s 2>/dev/null || hostname)
  name=${name:-$host-$(date +%Y%m%d-%H%M)}
fi
run=${ASM_RUNS:-asm-runs}/$name
if [[ -e $run ]]; then
  echo "asm-report: $run exists; pick another --name" >&2
  exit 2
fi
mkdir -p "$run/asm"
run=$(cd "$run" && pwd)

# a configuration's rustc flags, one per line
flags() {
  local target=${1%%/*} cpu=${1#*/} f
  local features=${cpu#*+}
  [[ $features == "$cpu" ]] && features=""
  cpu=${cpu%%+*}
  [[ $cpu == default ]] || echo "-Ctarget-cpu=$cpu"
  for f in ${features//+/ }; do echo "-Ctarget-feature=+$f"; done
  [[ $target == x86_64-* ]] && echo "-Cllvm-args=--x86-asm-syntax=intel"
  return 0
}

read -r commit dirty < <(source_id)
listings=()
entries=()
for config in "${configs[@]}"; do
  target=${config%%/*}
  mapfile -t rustflags < <(flags "$config")
  echo "asm-report: $config (${rustflags[*]})" >&2
  dir=target/asm/${config//\//-}
  RUSTFLAGS="${rustflags[*]}" CARGO_TARGET_DIR=$dir \
    cargo rustc --locked --quiet --release --lib --features asm-probes \
    --target "$target" -- --emit asm
  s=("$dir/$target/release/deps/"ephemeral_ecmh-*.s)
  # a fresh build leaves the last one's listing, which is current
  if ((${#s[@]} != 1)) || [[ ! -f ${s[0]} ]]; then
    echo "asm-report: expected one listing for $config, found ${s[*]}" >&2
    exit 1
  fi
  out=$run/asm/$config.s
  mkdir -p "$(dirname "$out")"
  cp "${s[0]}" "$out"
  listings+=("$config=$out")
  features=$(rustc --print cfg --target "$target" "${rustflags[@]}" |
    sed -n 's/^target_feature="\(.*\)"$/\1/p' | paste -sd, -)
  entries+=("$(jq -n --arg config "$config" --arg flags "${rustflags[*]}" \
    --arg features "$features" \
    '{config: $config, rustflags: $flags, target_features: ($features | split(","))}')")
done
read -r built _ < <(source_id)
if [[ $built != "$commit" ]]; then
  echo "asm-report: the sources changed during the run ($commit, then $built)" >&2
  exit 1
fi

jq -n \
  --arg name "$name" \
  --arg date "$(if [[ -z ${ASM_SOURCE_REV:-} ]]; then date -u +%Y-%m-%dT%H:%M:%SZ; fi)" \
  --arg host "$host" \
  --arg rustc "$(rustc -vV)" \
  --arg cargo_lock "$(sha256sum Cargo.lock | cut -d' ' -f1)" \
  --arg commit "$commit" \
  --argjson dirty "$dirty" \
  --argjson configs "$(printf '%s\n' "${entries[@]}" | jq -s .)" \
  '{name: $name, date: $date, host: $host, rustc: $rustc,
    cargo_lock_sha256: $cargo_lock, commit: $commit, dirty: $dirty,
    profile: "release", features: ["asm-probes"], configs: $configs}
    | if .date == "" then del(.date, .host) else . end' \
  >"$run/meta.json"

python "$ASM_REPORT" "$run" "${listings[@]}"
echo "asm-report: $run/summary.md" >&2
