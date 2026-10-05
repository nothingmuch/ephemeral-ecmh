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

For a namespace $N$, let $H_N$ map an item to a group element. The digest of a
multiset $S$ is

```text
D_N(S) = sum of H_N(x), with multiplicity, over x in S.
```

Insertion and removal add or subtract a hashed point. A prepared representation
of $H_N(x)$ is reused for every coded symbol the item maps to. The measured
costs are hashing to a reusable addend, signed accumulation, and encoding or
decoding accumulated points.

The encoding must represent the additive closure of the hash outputs, including
the identity, negatives, and sums, because coded symbols hold arbitrary sums.
Prime order is not required. A subgroup or quotient group is sufficient if
hashing, equality, accumulation, and canonical encoding all refer to the same
group. Not every bit string needs to be a valid encoding.

The [RIBLT workload](../src/riblt.rs) implements encoding and peeling with an
XOR item field, an integer count, and an ECMH checksum. It is a cost model for
these operations, not a specification of a wire protocol.

## Cost in repeated reconciliation

The RIBLT mapping includes an item in coded symbol $i$ with probability
$1/(1 + i/2)$, so an item enters $k(m) = 2(H_{m+1} - 1)$, about $2 \ln m$, of
the first $m$ symbols, and symbol $i$ holds about $n/(1 + i/2)$ of $n$ items:
the first symbols are dense. Reconciling a difference of $d$ items takes $m$ of
about $1.35d$ symbols [yang-gilad-alizadeh-2024]. Under three assumptions,

1. weights are expected weights,
2. every round reconciles the same number $d$ of new items, and
3. each round rebuilds its symbols from the full set,

one round over a set of $n$ items costs $n k(1.35d)$ additions, about
$2n \ln(1.35d)$. Hashing is per item and its result can be kept, so a set that
grows to $N$ items by rounds of $d$ new items costs $N$ hashes in all, but

```text
sum over rounds of n_j k(1.35d), about (N^2 / d) ln(1.35d) additions.
```

Under this model the number of additions grows quadratically in $N$ and the
number of hashes linearly: a round performs $d$ hashes against $2n \ln(1.35d)$
additions. Which dominates the elapsed time also depends on their measured unit
costs. Decoding adds about $d k(m)$ subtractions per round, and one equality
test per examined symbol, lower order than rebuilding when $n$ is much larger
than $d$. If instead each party keeps its symbol sums between rounds and adds
only new items (assumption 3 dropped), encoding falls to about $2N \ln(1.35d)$
additions, linear in $N$, with one hash per $2 \ln(1.35d)$ additions; decoding
is then of the same order. Kept hashes and sums are valid only within one
namespace: a new beacon value changes the salt, the curve and the mapping, and
everything is recomputed.

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
