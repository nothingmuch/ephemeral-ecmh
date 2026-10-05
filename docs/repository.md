---
type: Chapter
title: Repository organization
description: Where the field arithmetic, curves, group capabilities, benchmarks, reports, published results and reference programs live.
tags: [repository, source-code]
generated: { by: claude-code/claude-opus-5-5 }
contributors: [codex/gpt-6, claude-code/claude-fable-5-1]
status: draft
---

# Repository organization

- [Field arithmetic](../src/field/mod.rs): representations, arithmetic backends,
  and batch inversion.
- [Curves](../src/curve/mod.rs): group laws, coordinate systems, hash maps, and
  encodings.
- [Group capabilities](../src/group.rs): the operations consumed by
  [ECMH](../src/ecmh.rs), the RIBLT workload, and the benchmarks.
- [Benchmarks](../benches/) and [reporting tools](../report/).
- `results/`: the published benchmark runs, each in the compact form described
  in [Measurement methodology](methodology.md), from which the site renders one
  chapter per run. `bench-runs/`, where runs are recorded, is not under version
  control.
- [Tests](../tests/): known-answer tests against the reference programs'
  fixtures, point counts against PARI, and the certificate corpus. Every
  recorded fixture opens with how it was produced: the generated
  `tests/common/*.rs` with the Sage command that rewrites them, the corpus
  `tests/fixtures/cert_corpus.txt` with the test that records and checks it,
  and the report's `report/tests/fixtures/*.tsv` with the run and derivation
  they come from, as far as these were recorded.
- [Nix](../nix/): the pinned toolchains and the derivations for the package,
  checks, benchmark executables, reports, assembly listings, and this site;
  `nix run .#validate-commits -- RANGE` runs the checks at every commit of a
  range, as the [continuous integration](../.github/workflows/) does for each
  commit a push or pull request adds.
- [Documents](../docs/): the chapters; `docs/index.md` gives their order.
- [Reference programs](../sage/): the known-answer-test generators and order
  proofs the Rust fixtures come from (`kat*.sage`, `prove_orders.sage`,
  `agm_vectors.sage`, `sieve.sage`, `gf2_*.sage`), and the surveys behind the
  choice of fields and curves: candidate field sizes (`field_sizes.sage`), GLS
  fields (`gls_survey.sage`), Weil descent (`ghs.sage`, `ghs122.sage`), and
  complete Edwards curves over pseudo-Mersenne primes (`edwards_search.py`).
