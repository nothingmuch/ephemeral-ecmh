---
type: Code Generation Report
title: Field and curve code generation
description: Static instruction counts of field and point operations across 24 compiler targets and CPU models, from the asm-report run at revision 2d6c7e3c.
tags: [codegen, assembly, carry-less-multiplication, aarch64, x86-64]
sources:
  - id: asm-run
    resource: nix run .#asm-report at revision 2d6c7e3c, rustc 1.98.1, LLVM 22.1.8, release profile, --features asm-probes (24 configurations, 239 probes)
    title: asm-report matrix
generated: { by: claude-code/claude-opus-5-5 }
contributors: [codex/gpt-6, claude-code/claude-fable-5-1]
status: draft
---

# Field and curve code generation

This report examines optimized field arithmetic and curve operations for
implementation effects that could confound the ECMH and RIBLT comparison. It
records observations from the source revision and target configurations below.
Instruction counts identify code for inspection; they do not measure execution
time or establish a ranking of fields. Branch counts are static branch sites,
not executed branches. Native timing evidence is reported separately in the
benchmark results.

**Run:**
- `nix run .#asm-report` at 2d6c7e3c (rustc 1.98.1, LLVM 22.1.8, release
  profile, `--features asm-probes`).
- 24 configurations × 239 probes. Each probe is an `#[inline(never)]` wrapper in
  `src/asm_probes.rs` around one field operation, or one group add, sub, equals,
  prepare or decode, for each of 29 group representations.
- "Closure" is the probe plus every crate-local function it reaches. Calls out
  of the crate are listed but not counted: crrl's `set_invert`, Plonky3's
  `try_inverse`, the allocator. Merged functions are reported as aliases.
- The listings, `meta.json`, `summary.md` and `probes.tsv` are archived outside
  the repository, in `ephemeral-ecmh-archive/asm-runs/matrix-20261003`.
- `nix build .#asm-listings` builds the same matrix from the flake source and
  records the source revision in place of the host and date. A listing at any
  revision can therefore be rebuilt without the archive.
- Probe and row names use the benchmark report's family identifiers: the curve
  model (`binary`, `weier`, `edwards`, `twisted`) followed by the field, `127`
  for $\mathrm{GF}(2^{127})$ or for $\mathbb{F}_p$ with $p = 2^{127} - 1$,
  `61x2` for $\mathrm{GF}(p^2)$ with $p = 2^{61} - 1$, `goldilocks2` for
  $\mathrm{GF}(p^2)$ over the Goldilocks prime, and so on. [Candidate
  constructions](constructions.md) describes each family. All counts are of
  static instructions.

## Interpretation

The inspected instruction sequences are consistent with the stated arithmetic
counts. The following implementation and target effects qualify the
corresponding measurements. This inspection does not establish that each
implementation is optimal, or that the measured ordering generalizes to other
processors.

| finding | rows affected | status |
|---|---|---|
| Binary fields get carry-less multiply only from build flags | every $\mathrm{GF}(2^n)$ row off the carry-less targets | the binary-field timings hold only for builds with PMULL/PCLMULQDQ; `bench-bins-generic` measures the other backend |
| The single-item Edwards `prepare` performs an inversion | `group/*/prepare` for edwards127/107/61x2 | intended; the comparison uses the batched rows |

Other observations explained by the inspected code:
- The unscaled binary form has the fewest carry-less products per add and the
  most per equality test (§4).
- $\mathrm{GF}(2^{122})$'s multiply needs more carry-less products than
  $\mathrm{GF}(2^{127})$'s, from its representation (§3).
- fp127's short sqrt is a loop (§6).
- The Goldilocks panic stubs are unreachable bounds checks, and the Goldilocks²
  add's branches are upstream's rare-carry paths (§5, §6).
- The merged $\mathrm{GF}(2^{122})$ $\lambda$ adds follow from the formulas
  (§4).

## 1. Binary-field backend selection

| configuration | PMULL/PCLMULQDQ in `gf2_127::mul` | `binary127::add` closure |
|---|---|---|
| Apple, every model, `default` (apple-m1) included | 4 | 239 |
| Apple `generic` | 0 | 4707 |
| Android `default`, `cortex-x4` | 0 | 4707, 3324 |
| Android `cortex-{a55,a78,a720,x4,x925}+aes`, `oryon-1` | 4 | 243–255 |
| x86-64 `default`, `v3`, `v4` | 0 | 5951, 5885, 5176 |
| x86-64 `v3+pclmulqdq` | 4 | 225 |
| `znver4`, `znver5`, `sapphirerapids`, `arrowlake` | 4 | 193–225 |

- **Android:** In the inspected LLVM configurations, the Cortex models leave out
  the optional crypto extension. Without it, the binary-field implementations
  compile to their portable backends, carry-less products from integer
  multiplications of masked operands (crrl's `gfb254_m64` for
  $\mathrm{GF}(2^{127})$, `soft.rs` for $\mathrm{GF}(2^{109})$ and
  $\mathrm{GF}(2^{122}){}$): 13 to 19 times the instruction count. Among the
  inspected Android models, only `oryon-1` enables it by default. A deployment
  using the carry-less backend must establish hardware support and enable
  `+aes`, either as a baseline requirement or through runtime detection and
  dispatch.
- **x86-64:** none of the distribution baselines (v1–v4) includes PCLMULQDQ.
  Adding `+pclmulqdq` alone to v3 changes only the binary-field probes: the
  multiplies use PCLMULQDQ with the counts of the named cores, and
  `binary127::add` falls from 5885 to 225 instructions. The feature, not the CPU
  model, selects the backend.
- **Benchmarks:** `bench-bins` compiles for `apple-m4` and for
  `x86-64-v3+pclmulqdq`, `bench-bins-generic` for `generic` and `x86-64`; all
  four are in this matrix. A portable build uses the path of the last two.

For comparison, `edwards127::add` is 361–370 on every Android model, 362 on M4
and 524 on Zen 4. These are code sizes on two different paths, not speeds.

## 2. Apple: the configurations compiled are byte-identical through M4

Compiled: `aarch64-apple-darwin` at `default` (apple-m1), apple-a14, a15, a16,
a17, m4 and m5, and `generic`. No iOS target was built: the A-series models were
compiled for macOS. These are compiler-tuning probes, not iOS ABI coverage or
phone measurements.

- **default through apple-m4:** the whole listings are byte-identical.
- **apple-m5:** 14 of the 239 probe closures differ: the binary127 equals in
  every representation, and the decode, equals or prepare of the Edwards and
  Weierstrass families. In each, LLVM uses `umin`/`umax` where the others use
  `cmp` + `cset`/`csel`, with the neighboring instructions rescheduled.
- **generic:** no crypto extension, so the binary fields take the portable path
  (`binary127::add` 4707, as on Android `default`). The CPU model also changes
  the prime-curve code: 149 probes differ from apple-m4 in their closure
  metrics, among them the Edwards, $a = -1$ and Jacobian adds, by a few percent
  (edwards127::add 373 against 362). The generic benchmark subset therefore
  carries prime-field control rows.

Identical code removes the need to review each configuration's code separately,
not the need to measure each machine.

## 3. Multiplies per binary field (M4; PCLMULQDQ counts on x86 are the same)

| field | PMULL per mul / square | where the products go | `mul` closure |
|---|---|---|---|
| $\mathrm{GF}(2^{127})$, crrl | 4 / 2 | schoolbook $2 \times 2$; the trinomial folds with shifts | 21 |
| $\mathrm{GF}(2^{109})$, local | 6 / 4 | 4 products; two PMULL folds (pentanomial, elements kept unreduced in 128 bits) | 23 |
| $\mathrm{GF}(2^{122})$ tower, local | 7 / 6 | 3 Karatsuba products; 2 folds per $\mathrm{GF}(2^{61})$ half, paired in one vector | 30 |
| $\mathrm{GF}(2^{61})$, local | 3 / 3 | 1 product, 2 folds | 15 |

- **The code matches `src/field/gf2_122/pmull.rs`.** The standalone probe uses 4
  instructions to materialize the reduction constant; the inspected inlined
  callers hoist them.
- **Critical path:** $\mathrm{GF}(2^{122})$'s is three dependent PMULLs
  (product, then two folds); $\mathrm{GF}(2^{127})$'s is one PMULL plus shifts.
- **Per extended point add:** 44 PMULLs over $\mathrm{GF}(2^{127})$, 68 over
  $\mathrm{GF}(2^{109})$, 82 over $\mathrm{GF}(2^{122})$. Pentanomial reduction
  requires more PMULL instructions than Karatsuba saves in the multiplication
  step.

These counts arise from this representation and reduction schedule; they do not
establish a lower bound for arithmetic in the field, and they are not timings.

## 4. Binary fields and representations: $\mathrm{GF}(2^{109})$ against $\mathrm{GF}(2^{127})$, dense against GLS over $\mathrm{GF}(2^{122})$

$\mathrm{GF}(2^{109})$ and $\mathrm{GF}(2^{127})$ both hold an element in two
64-bit words and differ in the modulus. Degree 109 has no irreducible trinomial,
so the pentanomial $z^{109} + z^5 + z^4 + z^2 + 1$ folds with two PMULLs where
the trinomial $z^{127} + z^{63} + 1$ folds with shifts
([§3](#3-multiplies-per-binary-field-m4-pclmulqdq-counts-on-x86-are-the-same)),
and its $\sqrt{z}$ is dense, so a square root costs a multiplication where
$\sqrt{z} = z^{64} + z^{32}$ costs two shifts. The PMULL counts do not order
the timings: in the M4 native run at 7d4c6aa1 the $\mathrm{GF}(2^{109})$ add,
at 68 PMULLs against 44, took 12.8 against 14.6 ns at throughput, while its
square root took 2.89 against 1.26 ns.

PMULLs per operation's closure, M4 (x86 carry-less targets equal):

| operation | coordinates | $\mathrm{GF}(2^{127})$ | $\mathrm{GF}(2^{109})$ | $\mathrm{GF}(2^{122})$ | $\mathrm{GF}(2^{122})$ GLS |
|---|---|---|---|---|---|
| add | extended | 44 | 68 | 82 | 80 |
| add | $\lambda$ (and $w$ encoding) | 46 | 74 | 93 | 93 |
| add | unscaled | 32 | 50 | 61 | 61 |
| equals | extended | 8 | 12 | 14 | 14 |
| equals | $\lambda$ (and $w$ encoding) | 4 | 6 | 7 | 7 |
| equals | unscaled | 12 | 18 | 21 | 21 |
| prepare | extended | 10 | 16 | 20 | 18 |
| prepare | $w$ encoding, unscaled | 0 | 0 | 0 | 0 |

- **Formulas:** the add counts equal the formulas' costs at the per-field counts
  of §3: extended (Pornin's complete addition) $8M + 2S + 2 m_\beta$, unscaled
  (the same formula with $\beta$ cancelled) $7M + 2S$. Equality is $1M$ in
  $\lambda$ coordinates, $2M$ extended and $3M$ unscaled, which compares $S u$
  with $v X Z$ ($w^2 = S/(XZ) = v/u$).
- **Unscaled:** the fewest products per add and the most per equality test.
  RIBLT adds once per item and cell (the cost model's $k(m)$ adds per item);
  equality tests are paid per difference during peeling, which the decision
  model leaves out. The counts therefore point toward the unscaled form for this
  workload; the timing rows decide.
- **Prepare:** the $w$ encoding and unscaled addends are the decoded affine
  point itself, so preparation is a copy (3–6 instructions). The $\lambda$
  form's prepare computes $\lambda = x + y/x$: one inversion (crrl's
  `set_invert`, uncounted, over $\mathrm{GF}(2^{127})$; 170–220 PMULLs in the
  crate over $\mathrm{GF}(2^{109})$ and $\mathrm{GF}(2^{122}){}$).
- **Merged:** rustc merges `binary122_lambda::add`, `sub` and `equals`
  (`Curve<Gf>`) with the `binary122_gls_w` ones (`Curve<Gf61>`). `Curve<K>`
  selects the representation of $B$ and its derived constants $b$ and $\beta$;
  $\lambda$ mixed addition uses none of them, so identical code is expected.
- **GLS:** only the extended add and prepare differ, by the multiplication by
  $\beta$: 82 versus 80 PMULLs, 20 versus 18. These counts do not imply a
  corresponding speedup; the GLS endomorphism concerns scalar multiplication,
  outside the primary RIBLT workload.

## 5. Point adds: branches and stack

Static counts over the add's closure; ranges are over the configurations of each
architecture.

- **Branches:**
  - Extended and unscaled binary adds, and the Edwards and $a = -1$ adds over
    fp127, fp107, $\mathrm{GF}(p^2)$ for $p = 2^{61} - 1$ and $p = 2^{64} - 59$,
    are branch-free.
  - $\lambda$ adds have 4–5 branches: $x_P = x_Q$ and the identity, exceptional
    cases of the formulas.
  - Weierstrass projective adds (Renes–Costello–Batina) have one branch, the
    test for an $O$ addend, which has no affine representation (two over fp107
    and on x86 over $\mathrm{GF}(p^2){}$). Jacobian adds have 4–7 and include
    the inlined doubling (746–880 instructions on M4, against 528–670 for the
    complete projective form), so only timings can compare them.
  - twisted128's add has 4 branches on M4, from fp128's carry handling, which
    its source documents as variable-time; its inputs are public.
  - twisted-goldilocks2's add has 61 branch sites on aarch64 and 20 on x86. They
    are Plonky3's Goldilocks reductions, which branch on carries and borrows
    that upstream marks rare (`branch_hint`). They are not exceptional cases of
    the point formulas, and a static count does not say how often they are
    taken.
- **Stack:**
  - The binary addition probes on carry-less targets have at most two stack
    accesses (on aarch64, none outside the $\mathrm{GF}(2^{122})$ adds; on x86,
    a push and pop of one callee-saved register); their field elements use
    vector registers. The portable path spills: 111–767 accesses on x86 without
    PCLMULQDQ.
  - The prime-curve additions spill on x86: 34–110 stack accesses, against 8–31
    on aarch64 apart from weier61x2 (22–55). Register pressure from the
    multi-limb representations is a platform-dependent implementation effect in
    this comparison.

## 6. Prime fields and extensions

- **Multiplication:** the M4 prime multiply rows have no calls out of the crate.
  Their reduction sequences differ; counts alone do not determine which field is
  faster.
- **fp127 sqrt:** the 71-instruction body includes a loop initialized to 125.
  The 26× difference in static closure size relative to fp128 (1878) does not
  give the ratio of executed instructions or timings.
- **Goldilocks:**
  - Inversion is an external call to Plonky3's `try_inverse`, so closure counts
    omit its body.
  - The base field's and Goldilocks²'s `panic_bounds_check` paths come from the
    upstream Tonelli–Shanks table indexing, which Goldilocks²'s square root
    reaches through the base field's. The source establishes `i < m <= 32`, so
    they are unreachable, not a runtime cost.
- **Extensions:** $\mathrm{GF}(p^2)$ for $p = 2^{61} - 1$ and $p = 2^{64} - 59$,
  and Goldilocks², each carry an Edwards or $a = -1$ group and, for
  $p = 2^{61} - 1$, a Weierstrass group; their decode closures (1.1–1.6k on M4)
  include the square root in $\mathrm{GF}(p^2)$.

## 7. Hash-to-curve side

- **`prepare`:**
  - Preparation for the variants whose affine representation is already the
    addend consists of copying it: 3–6 instructions (binary $w$ encoding and
    unscaled, Weierstrass).
  - For extended binary forms it costs 73–100 on M4 (computing $T$).
  - For Edwards, the single-item probe is 1.7–3.6k instructions and allocates.
    It delegates to `to_cached_batch` on one item. The inversion is intrinsic to
    the coordinate conversion (the affine is the Montgomery point); the Vec is
    an implementation choice. RIBLT's encode and batched peel use
    `prepare_batch`, which amortizes the inversion.
- **`decode`:**
  - $\mathrm{GF}(2^{127})$'s is 193–274 with PMULL; its inversion is crrl's
    `set_invert`, uncounted here.
  - $\mathrm{GF}(2^{122})$'s and $\mathrm{GF}(2^{109})$'s are 640–950 on
    carry-less targets, all in the crate (half-trace and inversion). On x86
    without PCLMULQDQ they are 1900–2210 and 1340–1490.

## 8. Coverage limitations

- The inspected matrix has no probes for sqrt_ratio or the whole hash-to-addend
  batch path.
- No iOS targets, and no measured phones: core-model names do not identify a
  Snapdragon or Dimensity device.
- Out-of-crate callees are not counted (Run, above).
