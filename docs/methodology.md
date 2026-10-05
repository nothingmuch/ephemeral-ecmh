---
type: Chapter
title: Measurement methodology
description: How operation latency, throughput, batching, and benchmark provenance are measured and reported.
tags: [benchmarks, methodology, nix]
sources:
  - id: bibliography
    resource: ../references.bib
    title: Bibliography, annotated in literature.md
generated: { by: claude-code/claude-opus-5-5 }
contributors: [codex/gpt-6, claude-code/claude-fable-5-1]
status: draft
---

# Measurement methodology

The benchmarks report operation latency, throughput with independent accumulators, and batch costs. Batch size and representation conversion are part of the workload. Each result records its source revision, compiler, target features, machine, and measurement method.
With Nix installed, the checks run from the repository root:

```sh
nix flake check
```

The book builds with:

```sh
nix build .#site
```

A reference measurement run:

```sh
nix run .#bench-run -- --name experiment compare
```
Each run is identified by a fresh UUIDv7, which sorts by start time and needs no
coordination between machines. The executables come from a Nix derivation built
for a named target CPU: `bench-bins` is the architecture's default (apple-m4 on
aarch64, x86-64-v3 with PCLMULQDQ on x86_64), and `bench-bins-generic` is the
baseline without carry-less multiplication. `BENCH_BINS` selects another such
build. The runner refuses executables that need a feature the host lacks, times
the suites serially, and records metadata and raw Criterion output under
`bench-runs/<id>/`, together with a report. Other builds and measurements should
not run during timing. To show whether they did, the runner samples the load
averages every 5 seconds, the interval at which the kernel updates them, into
`load.tsv`, beginning and ending with 30 seconds idle that measure the machine's
own load; the report gives each phase's mean and recovers the threads runnable
from successive 1-minute averages where the kernel does not report them.

A report can be regenerated from an existing run whose benchmark identifiers the
reporter recognizes; runs recorded under an earlier schema need the reporter of
their revision:

```sh
nix run .#bench-report -- bench-runs/ID report-output
```

Published runs are kept in the repository in a compact form, one directory
`results/<id>/` per run. `benchmarks.csv` has one row per benchmark:
Criterion's group, function and parameter, the estimate the report reads (the
slope where Criterion fits one, otherwise the mean), its 95% confidence
interval, and the declared throughput. `meta.json` is the run's metadata without
the host name, `curvegen.csv`, present when the run timed curve selection, the
searches per seed, and `load.tsv` the load samples. Criterion's output, the
figures and the report are not kept; the classification into layers, families
and operations is not stored either, and is applied by the reporter of the
revision that renders the run. Only a run of a committed source is exported,
since no published revision reproduces one with uncommitted edits. Each run is
exported into its own directory under `results/`, which must not exist yet:

```sh
nix run .#bench-report -- --export bench-runs/ID results/ID
```

The reporter renders the directory as it renders the run's Criterion output,
to the same tables and figures, except that only the default estimate is
available. The site build renders each run under `results/` as a chapter of
the book's evidence, whose first paragraph states the commit, machine, build
and start time recorded in `meta.json`.

## Reading a run report

Values are Criterion's estimate and 95% confidence interval, per element where
the benchmark declares an element count and otherwise per iteration. A table
colours each value by its ratio to its column's minimum on matplotlib's inferno
scale, from pale yellow at the minimum to purple at 16 times or more. The
minimum is bold when no other interval overlaps its own; otherwise the values
that overlap it are tied, in italics. Overlap leaves an order unresolved. A tag
after a row's name marks work unlike the rest of its column: addition formulas
that are not complete (they branch on the identity and on doubling), or a term
another suite timed. The XOR baseline is excluded from column minima. Tooltips
name the benchmark each value reads.

A modeled cost sums separately measured operations, and its interval combines
their half-widths in quadrature; it has no independently measured coverage, and
neither has a range over curve-selection seeds, which spans their point
estimates. The insertion model uses the mean number $k(m)$ of coded symbols an
item maps to and
compatible representations: a recipe either hashes, prepares the hash's output,
and adds the prepared addend, or hashes straight to the addend and skips the
preparation. Pornin's map returns the extended binary addend, and its variants
to $(x, \lambda)$ and $(u, v)$ the λ and $u$ addends; Elligator 2 and SSWU
return affine points, which the measured preparation takes. Preparation and
addition are timed on try-and-increment inputs, so a recipe with another map is
a model, not an end-to-end measurement. Every recipe also adds the salted
map digest and the walk of the mapping's indices below $m$, timed alone by
riblt.mapping and the same for every family, which is why the model's ratios
between families are smaller than those of their hashes and additions. The
model omits cell allocation, the key XOR and the counts; the measured RIBLT
encoding includes them.

Preparation occurs once per source item, addition once per cell update, and
subtraction negates an addend before adding it. Encoding and decoding occur per
transmitted and received cell. Throughput runs eight independent accumulators,
latency one dependent chain. A batched operation reports time per element at
its batch size, 1024 unless a sweep states another. Batch inversion
(Montgomery's trick) replaces $n$ inversions with one and $3(n - 1)$
multiplications; binary hashing shares them across the pending candidates of
each round, and preparation, normalization and decoding where the
representation needs an inversion. The prime-field hashes are bound by a square
root per item, which batching does not share. An encoder can batch all its
source items, a decoder only the cells of the current peeling wave. The binary
w-codec families hash to reusable addends, so their preparation is a copy.

## Contributing a run

A run on another machine is added to `results/` by pull request. The runner
records the commit of the checkout it runs in, and only a run of a committed
source is exported, so the run happens in a clean clone of a published commit:

```sh
git clone https://github.com/nothingmuch/ephemeral-ecmh
cd ephemeral-ecmh
nix build --no-link .#bench-bins
nix run .#bench-run -- --profile full
```

The run's directory, `bench-runs/ID`, is named by its UUIDv7; the report
describes the machine and build from `meta.json`. Building beforehand keeps
compilation out of the timed window. The full profile (3 s warm-up, 5 s
measurement and 100 samples per benchmark) takes several hours over all suites;
runs under the default quick profile are preliminary. On x86_64 the native build
needs AVX2 and PCLMULQDQ; `BENCH_BINS` and suite arguments select another build
or a subset, as above. The run is then exported:

```sh
nix run .#bench-report -- --export bench-runs/ID results/ID
```

The export is committed on a branch, whose pull request adds only
`results/ID/`:

```sh
git switch -c runs/ID
git add results/ID
git commit -m "runs: ID"
```
