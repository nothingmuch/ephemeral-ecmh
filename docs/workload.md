---
type: Chapter
title: Workload
description: The multiset digest, the operations a RIBLT performs on it, and their counts in repeated reconciliation.
tags: [riblt, workload, cost-model]
sources:
  - id: bibliography
    resource: ../references.bib
    title: Bibliography, annotated in literature.md
generated: { by: claude-code/claude-opus-5-5 }
contributors: [codex/gpt-6, claude-code/claude-fable-5-1]
status: draft
---

# Workload

## Comparison with the reference implementations

The Go benchmarks of
[yangl1996/riblt#3](https://github.com/yangl1996/riblt/issues/3) run at pinned
revisions of [riblt-ecmh]: the upstream XOR checksum, the ristretto255 ECMH with
the reference 64-bit mapping, and the same with a ChaCha8 mapping. All three are
run by `nix run .#riblt-go-bench`. `riblt.stream` follows their
BenchmarkEncodeAndDecode. The decoder holds $d/2$ items of its own and the $d$
common ones, and the encoder holds $d/2$ of its own and the $d$ common ones;
symbols are 64 bytes, and insertion is outside the timed routine. `riblt.encode`
per item corresponds to BenchmarkSketchAddSymbol at smaller $m$. ristretto255
runs as a family of the RIBLT suites, so that a small curve and the fixed group
are measured in one harness.

The implementations differ in ways that affect the comparison:

- Items are hashed with tagged, salted SHA-256 (two digests for ristretto255),
  where riblt-ecmh hashes the bare symbol with SHA-512. The upstream XOR
  checksum is SipHash-2-4 under a fixed key; `xor-siphash.64` keys it with the
  salt.
- The rateless decoder recovers an item with the addend its purity check
  prepared; the ristretto255 version of riblt-ecmh hashes the item to the curve
  again when it recovers it.
- The ristretto255 addend is the point itself, since curve25519-dalek does not
  expose its cached form.
- The Rust default mapping uses xoshiro256++ seeded by the whole 32-byte map
  digest. Rows marked `map=chacha8` use ChaCha8, and `map=mcg64` rows use the
  reference 64-bit generator. The Go ChaCha8 version of riblt-ecmh replaced
  that generator because two items with equal 64-bit seeds have equal
  mappings, and the decoder cannot separate them
  ([yangl1996/riblt#3](https://github.com/yangl1996/riblt/issues/3);
  [Known weaknesses](ecc_security.md#known-weaknesses)). Even Rust's ChaCha8
  comparison uses different output from Go's ChaCha8Rand, so their cells
  differ. An unmarked mapping row uses the default of its recorded source
  revision; historical runs are not relabelled when the default changes.

The checksum work is hashing items to addends, adding addends to cells, checking
purity while peeling, and encoding or decoding points. Rows that differ only in
the group, in the same Rust harness with the same mapping, workload, batching
and build, measure checksum substitution, including the choice of hash and
point representation. A time ratio against Go riblt-ecmh also includes the
mapping generator, the reuse policies above and the port from Go to Rust;
these implementation gains are not attributed to this project. Replacing
ChaCha8 with xoshiro256++ is of the same kind: it is chosen to make the
checksum's cost clearer in the end-to-end figures. Matched `map=chacha8` rows
provide a comparison of the generator's effect once the new clean runs are
measured; historical timings do not quantify this change.
