# shellcheck shell=bash
# Criterion suites, one after another, into a run directory that records the
# machine, then the report:
#
#   nix run .#bench-run -- [--profile quick|full] [SUITE...] [-- CRITERION ARGS]
#
# SUITE defaults to every suite. The run goes to bench-runs/ID under the
# current directory, or under BENCH_RUNS. ID is a fresh UUIDv7, so runs
# from any machine can be published side by side without choosing names,
# and sort by when they started; meta.json describes the run.
#
# The quick profile (the default) times with the benches' shared
# configuration, benches/common's config(); full passes criterion's own
# defaults as flags, 3 s warm-up, 5 s measurement, 100 samples and 100k
# resamples, for about five times as long. CRITERION ARGS override either.
# The run directory holds:
#
#   meta.json    machine, toolchain, source, the group suite's certified
#                fixtures (r, cofactor, automorphism group order), and the
#                RIBLT plan
#   riblt-plan.tsv  with riblt: each family's hash and scope, chosen from the
#                group operations the suites before it timed
#   load.tsv     the load averages every BENCH_LOAD_INTERVAL seconds
#                (default 5, the interval at which both kernels update
#                them), with what ran then: idle for BENCH_IDLE seconds
#                (default 30) before the first suite and after the last,
#                each suite, and the curvegen searches; on Linux also the
#                threads runnable at the sample
#   bins         a GC root for the executables it timed
#   criterion/   criterion's output (CRITERION_HOME)
#   log/SUITE    each suite's output
#   curvegen.csv with curvegen: single-shot searches per seed and family
#                (examples/curvegen_times.rs), CURVEGEN_SEEDS of them
#                (default 4 quick, 16 full); all families took about
#                four minutes a seed on an Apple M4
#   report/      bench-report on the run
#
# The executables are packages.bench-bins, which nix builds before this
# starts, from the flake's own source, for the architecture's default
# target CPU (apple-m4, x86-64-v3 with pclmulqdq); or BENCH_BINS, another
# such output: another target's (bench-bins-generic, the baseline without
# carry-less multiplication) or another revision's:
#
#   BENCH_BINS=$(nix build --no-link --print-out-paths .?rev=REV#bench-bins-generic)
#
# Nothing builds while anything is timed, and the suites run serially:
# don't build or bench anything else on the machine meanwhile.

default_suites=(compare group)
profile=quick
suites=()
while (($#)); do
  case $1 in
  --profile)
    profile=$2
    shift 2
    ;;
  --)
    shift
    break
    ;;
  -*)
    echo "bench-run: unknown option $1" >&2
    exit 2
    ;;
  *)
    suites+=("$1")
    shift
    ;;
  esac
done
((${#suites[@]})) || suites=("${default_suites[@]}")
case $profile in
quick)
  profile_args=()
  ;;
full)
  profile_args=(--warm-up-time 3 --measurement-time 5 --sample-size 100 --nresamples 100000)
  ;;
*)
  echo "bench-run: unknown profile $profile; quick or full" >&2
  exit 2
  ;;
esac
# criterion rejects a repeated option, so the profile's give way to the
# command line's
args=()
for ((i = 0; i < ${#profile_args[@]}; i += 2)); do
  [[ " $* " == *" ${profile_args[i]}"[\ =]* ]] || args+=("${profile_args[@]:i:2}")
done
set -- "${args[@]}" "$@"
criterion_args="$*"

os=$(uname -s)
arch=$(uname -m)
host=$(hostname -s 2>/dev/null || hostname)
# RFC 9562 version 7: 48 bits of Unix milliseconds, the version, the
# variant, and 74 random bits
uuid7() {
  local ms r
  ms=$(printf '%012x' "$(date +%s%3N)")
  r=$(od -An -N10 -tx1 /dev/urandom | tr -d ' \n')
  printf '%s-%s-7%s-%x%s-%s\n' "${ms:0:8}" "${ms:8:4}" "${r:0:3}" \
    $((0x${r:3:1} & 3 | 8)) "${r:4:3}" "${r:7:12}"
}
id=$(uuid7)
run=${BENCH_RUNS:-bench-runs}/$id
# a run directory holds one run: criterion keeps what it finds there, so a
# second run into it would mix in the first's results under new metadata
if [[ -e $run ]]; then
  echo "bench-run: $run exists" >&2
  exit 2
fi

bins=${BENCH_BINS:-$BENCH_BINS_BUILT}
for s in "${suites[@]}"; do
  [[ -x $bins/bin/$s ]] || {
    echo "bench-run: $bins has no suite $s" >&2
    exit 2
  }
done
target=$(<"$bins/target")
rustflags=$(<"$bins/rustflags")
target_features=$(<"$bins/target-features")
# Every feature the target CPU turns on, by the host's name for it: on a
# CPU without one the executables may stop at an illegal instruction.
# Features with no name here are not checked: neon and the x86-64 v1
# baseline are architectural, and lor, pan, pmuv3, ras, vh and ssbs are
# system-level, never emitted by codegen (the M4 reports FEAT_SSBS 0
# although apple-m4 enables ssbs).
host_name() {
  if [[ $os == Darwin && $arch == arm64 ]]; then
    case $1 in
    aes) echo FEAT_AES FEAT_PMULL ;; bf16) echo FEAT_BF16 ;;
    bti) echo FEAT_BTI ;; crc) echo FEAT_CRC32 ;; dit) echo FEAT_DIT ;;
    dotprod) echo FEAT_DotProd ;; dpb) echo FEAT_DPB ;; dpb2) echo FEAT_DPB2 ;;
    fcma) echo FEAT_FCMA ;; fhm) echo FEAT_FHM ;; flagm) echo FEAT_FlagM ;;
    fp16) echo FEAT_FP16 ;; frintts) echo FEAT_FRINTTS ;; i8mm) echo FEAT_I8MM ;;
    jsconv) echo FEAT_JSCVT ;; lse) echo FEAT_LSE ;; lse2) echo FEAT_LSE2 ;;
    paca | pacg) echo FEAT_PAuth ;; rcpc) echo FEAT_LRCPC ;;
    rcpc2) echo FEAT_LRCPC2 ;; rdm) echo FEAT_RDM ;; sb) echo FEAT_SB ;;
    sha2) echo FEAT_SHA256 ;; sha3) echo FEAT_SHA3 FEAT_SHA512 ;;
    sme) echo FEAT_SME ;; sme2) echo FEAT_SME2 ;;
    # no sysctl on CPUs without them, so they read as absent
    sve) echo FEAT_SVE ;; sve2) echo FEAT_SVE2 ;;
    esac
  elif [[ $arch == arm64 || $arch == aarch64 ]]; then
    case $1 in
    aes) echo aes pmull ;; bf16) echo bf16 ;; bti) echo bti ;; crc) echo crc32 ;;
    dit) echo dit ;; dotprod) echo asimddp ;; dpb) echo dcpop ;;
    dpb2) echo dcpodp ;; fcma) echo fcma ;; fhm) echo asimdfhm ;;
    flagm) echo flagm ;; fp16) echo fphp asimdhp ;; frintts) echo frint ;;
    i8mm) echo i8mm ;; jsconv) echo jscvt ;; lse) echo atomics ;;
    lse2) echo uscat ;; paca) echo paca ;; pacg) echo pacg ;;
    rcpc) echo lrcpc ;; rcpc2) echo ilrcpc ;; rdm) echo asimdrdm ;;
    sb) echo sb ;; sha2) echo sha1 sha2 ;; sha3) echo sha3 sha512 ;;
    sve) echo sve ;; sve2) echo sve2 ;; sme) echo sme ;;
    esac
  else
    case $1 in
    sse3) echo pni ;; sse4.1) echo sse4_1 ;; sse4.2) echo sse4_2 ;;
    lzcnt) echo abm ;; bmi1) echo bmi1 ;; cmpxchg16b) echo cx16 ;;
    pclmulqdq | aes | avx | avx2 | bmi2 | f16c | fma | movbe | popcnt | xsave | \
      ssse3 | adx | sha | vpclmulqdq | vaes | gfni | avx512*) echo "${1//./_}" ;;
    esac
  fi
}
host_has() {
  if [[ $os == Darwin && $arch == arm64 ]]; then
    [[ $(sysctl -n "hw.optional.arm.$1" 2>/dev/null) == 1 ]]
  elif [[ $os == Darwin ]]; then
    # Intel macs name them in upper case, sse4.1 as SSE4.1, abm as LZCNT
    local f=$1
    case $f in pni) f=sse3 ;; sse4_1) f=sse4.1 ;; sse4_2) f=sse4.2 ;; abm) f=lzcnt ;; bmi1) f=bmi1 ;; esac
    sysctl -n machdep.cpu.features machdep.cpu.leaf7_features machdep.cpu.extfeatures 2>/dev/null |
      tr '[:upper:]' '[:lower:]' | grep -qw -- "$f"
  else
    grep -Eq "^(flags|Features)[[:space:]]*:.* $1( |$)" /proc/cpuinfo
  fi
}
IFS=, read -ra recorded <<<"$target_features"
for f in "${recorded[@]}"; do
  for h in $(host_name "$f"); do
    if ! host_has "$h"; then
      echo "bench-run: $bins is built for $target, which needs $f; this CPU lacks it ($h)" >&2
      exit 2
    fi
  done
done

# Every id the filter selects must be one the report classifies: refuse
# now, not after the timing. Outside a pipeline, errexit stops at the
# first suite that cannot list.
listed=$(mktemp)
trap 'rm -f "$listed"' EXIT
for s in "${suites[@]}"; do
  "$bins/bin/$s" --bench --list "$@" >>"$listed"
done
sed -n 's/: benchmark$//p' "$listed" | bench-report --check-ids || exit 2

mkdir -p "$run/log"
run=$(cd "$run" && pwd)
# the run keeps the executables it timed, and the libpari they load
nix-store --add-root "$run/bins" --realise "$bins" >/dev/null
export CRITERION_HOME=$run/criterion

if [[ $os == Darwin ]]; then
  cpu=$(sysctl -n machdep.cpu.brand_string)
  cores=$(sysctl -n hw.ncpu)
  memory=$(sysctl -n hw.memsize)
  os_version="macOS $(sw_vers -productVersion)"
  governor=""
else
  # awk, not head: under pipefail a closed pipe fails the script
  cpu=$(awk -F': ' '/^(model name|Model)[[:space:]]*:/ { print $2; exit }' /proc/cpuinfo)
  # aarch64 Linux names the core by part number alone
  [[ -n $cpu ]] || cpu=$(awk -F': ' '/^CPU part/ { print "CPU part " $2; exit }' /proc/cpuinfo)
  cores=$(nproc)
  memory=$(($(sed -n 's/^MemTotal: *\([0-9]*\) kB/\1/p' /proc/meminfo) * 1024))
  os_version=$(sed -n 's/^PRETTY_NAME="\{0,1\}\([^"]*\)"\{0,1\}$/\1/p' /etc/os-release 2>/dev/null || true)
  governor=$(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || true)
  if [[ -n $governor && $governor != performance ]]; then
    echo "bench-run: the CPU frequency governor is $governor, not performance;" \
      "times will be noisier" >&2
  fi
fi
# the flake's revision, which ends in -dirty for uncommitted edits
commit=$(<"$bins/source")
dirty=false
[[ $commit == *-dirty ]] && dirty=true
group_fixtures='{}'

# the 1, 5 and 15 minute load averages and, where the kernel reports it,
# the threads runnable now other than the reader
loadavg() {
  if [[ $os == Darwin ]]; then
    sysctl -n vm.loadavg | awk '{ gsub(/[{}]/, ""); print $1 "\t" $2 "\t" $3 "\t" }'
  else
    awk '{ split($4, r, "/"); print $1 "\t" $2 "\t" $3 "\t" r[1] - 1 }' /proc/loadavg
  fi
}
phase=$run/log/phase
during() {
  printf '%s' "$1" >"$phase.new"
  mv "$phase.new" "$phase"
}
sample() {
  printf 'time\tphase\tload1\tload5\tload15\trunnable\n'
  while :; do
    printf '%s\t%s\t%s\n' "$(date +%s)" "$(<"$phase")" "$(loadavg)"
    sleep "${BENCH_LOAD_INTERVAL:-5}"
  done
}
idle() {
  during "idle $1"
  sleep "${BENCH_IDLE:-30}"
}

meta() {
  jq -n \
    --arg id "$id" \
    --arg started "$started" \
    --arg finished "${1-}" \
    --arg host "$host" \
    --arg os "$os_version" \
    --arg kernel "$(uname -sr)" \
    --arg arch "$arch" \
    --arg cpu "$cpu" \
    --argjson cores "$cores" \
    --argjson memory "$memory" \
    --arg governor "$governor" \
    --arg rustc "$(<"$bins/rustc")" \
    --arg target "$target" \
    --arg rustflags "$rustflags" \
    --arg bins "$bins" \
    --arg target_features "$target_features" \
    --arg commit "$commit" \
    --argjson dirty "$dirty" \
    --arg suites "${suites[*]}" \
    --arg profile "$profile" \
    --arg criterion_args "$criterion_args" \
    --argjson group_fixtures "$group_fixtures" \
    '{id: $id, started: $started, finished: $finished, host: $host,
      os: $os, kernel: $kernel, arch: $arch, cpu: $cpu, cores: $cores,
      memory: $memory, governor: $governor, rustc: $rustc,
      target: $target, rustflags: $rustflags, bins: $bins, target_features: ($target_features | split(",")),
      commit: $commit, dirty: $dirty, suites: ($suites | split(" ")),
      profile: $profile, criterion_args: $criterion_args, group_fixtures: $group_fixtures}
     | with_entries(select(.value != ""))' >"$run/meta.json"
}

started=$(date -u +%Y-%m-%dT%H:%M:%SZ)
meta
echo "bench-run: $run ($cpu; $bins, for $target)" >&2

during "idle before"
sample >"$run/load.tsv" &
sampler=$!
trap 'rm -f "$listed"; kill "$sampler" 2>/dev/null || true' EXIT
idle before

for s in "${suites[@]}"; do
  echo "bench-run: $s" >&2
  during "$s"
  # --bench, as cargo passes it: criterion times only when it's given
  "$bins/bin/$s" --bench "$@" 2>&1 | tee "$run/log/$s"
  if [[ $s == group ]]; then
    # Read declarations emitted by the executable whose timings were recorded.
    # Missing declarations remain unknown; source revision alone implies none.
    group_fixtures=$(jq -Rn '
      [inputs | select(startswith("group-fixture\t")) | split("\t")
        # name, r as the certificate verified it, a probable prime (a string:
        # it exceeds a double),
        # the cofactor and the automorphism group order rho can use
        | if length == 5 and .[1] != "" and (.[2] | test("^r=[0-9]+$"))
             and (.[3] | test("^cofactor=[0-9]+$"))
             and (.[4] | test("^automorphisms=[0-9]+$"))
          then {key: .[1], value: {r: .[2][2:],
                                   cofactor: (.[3][9:] | tonumber),
                                   automorphisms: (.[4][14:] | tonumber)}}
          else error("invalid group-fixture declaration") end]
      | group_by(.key)
      | map(if (map(.value) | unique | length) == 1 then .[0]
            else error("conflicting group-fixture declarations") end)
      | from_entries' "$run/log/$s")
  fi
  meta
done

idle after
kill "$sampler"
wait "$sampler" 2>/dev/null || true
rm -f "$phase"
meta "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
bench-report "$run" "$run/report"
