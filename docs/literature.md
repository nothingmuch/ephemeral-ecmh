---
type: Annotated Bibliography
title: Annotated literature
description: For each key of references.bib, what the work contributes to the evaluation of ECMH-based RIBLT checksums over per-namespace curves; in the book, headed by the formatted entry.
tags: [bibliography, ecmh, riblt, elliptic-curves, point-counting, weil-descent]
sources:
  - id: bibliography
    resource: ../references.bib
generated: { by: claude-code/claude-opus-5-5 }
contributors: [codex/gpt-6, claude-code/claude-fable-5-1]
status: draft
---

# Annotated literature

One entry per key in [`references.bib`](../references.bib), stating why the
work bears on the question this repository asks: how efficient an ECMH checksum
for RIBLT cells can be made when the curve is derived from a short-lived public
salt (or a shared secret), with 48–64-bit collision resistance and 16-byte
identifiers. Within a section, entries are ordered by role: foundational
results, attacks and bounds, methods used here, alternatives considered. The
chapters cite works by key; in the book, each entry is headed by its citation
as formatted from the bibliography, and each key cited links to its entry.

The entries keep four distinctions apart: an encoding, a hash indifferentiable
from a random oracle, and a weak encoding with a collision-resistance
reduction; a field primitive, a curve map, and a prepared accumulator operand;
fixed-time and variable-time implementation; and setup reusable per field
against setup per epoch. The module boundaries and the benchmark rows follow
the same distinctions.

## Curve generation, point counting, Weil descent and discrete-logarithm security

### Generic discrete-logarithm bounds and attacks

- `bernstein-lange-2012`, `bernstein-lange-2013`. With a precomputed table of
  size about $N^{2/3}$, each further discrete logarithm in a group of size $N$
  costs about $N^{1/3}$ operations. A fixed curve would let an attacker amortize
  this across epochs. Deriving the curve from the salt prevents direct reuse of
  curve-specific tables; field-level arithmetic and work on the candidate space
  remain reusable. These two papers motivate a curve per namespace.

### Point counting and curve generation

- `ansi-x9-62-1998`, `rfc5639`, `baigneres-et-al-2015`,
  `lenstra-wesolowski-2017`. Curves derived from public randomness. X9.62's
  verifiably random curves and Brainpool's verifiably pseudo-random ones hash a
  published seed, which `bernstein-et-al-2015` shows still leaves the seed's
  author a choice. The Million Dollar Curve takes its seed from several national
  lotteries and turns it into parameters by a deterministic filter, for one
  curve meant to replace P-256 and Curve25519. Lenstra and Wesolowski's trx
  derives a stream of curves from the output of their beacon, unicorn: the first
  acceptable candidate in a fixed enumeration, with data to show that the
  earlier candidates fail, which is the shape of this repository's certificate.
  Both provide curves for general use at the 128-bit level, where trx takes
  hours per curve; here a new beacon value selects a smaller curve for each
  namespace.
- `bonneau-clark-goldfeder-2015`. Bitcoin block hashes as a beacon: a lower
  bound on each block's min-entropy, and the cost to a miner of biasing the
  output by withholding blocks, which is the manipulation this study leaves to
  the protocol.

## Hash-to-curve, encodings, ECMH and inversion

### Security contract and foundational constructions

- `bellare-micciancio-1997`. Establishes the incremental multiset-hashing
  framework underlying ECMH. It supplies the algebraic security question for
  cell checksums; a fast map and group operation must be assessed within that
  question rather than by mapping speed alone.
- `wegman-carter-1981`. Defines almost-universal hash families and applies them
  to authentication and to testing sets for equality under insertion and
  deletion. The salted projection of an item's unsalted identifier is such a
  family: its collision bound holds for identifiers fixed before the salt, and
  after the salt the identifiers' being hash outputs leaves only a generic
  search ([Adversary](problem.md#adversary)).
- `wagner-2002`. The generalized birthday algorithm combines many independently
  chosen terms into an additive relation; for an $n$-bit modulus and about
  $2^{\sqrt{n}}$ lists it takes about $2^{2\sqrt{n}}$ operations, subexponential
  in $n$. Its list-merging step needs a group representation in which partial
  sums can be matched on bits, as in $\mathbb{Z}/2^n$ or $\mathbb{Z}/r$.
  Elliptic-curve groups have no such representation, so without discrete
  logarithms the $k$-sum bound is not the attack cost for these checksums; with
  them, the checksum is AdHash modulo $r$ and the algorithm applies.
- `phan-wagner-2006`. Collisions in randomize-then-combine hashes of
  pair-chained message blocks: repeated blocks, coinciding padding and cyclic
  chaining give distinct messages with equal hashes, in several cases with
  equal multisets of pair terms, so that modular addition fails as XOR does.
  They follow from the encoding of a message into terms, not from the
  combining operation; here each item is hashed whole to one term and the
  digest of a multiset is meant to be independent of order, so the attacks do
  not carry over. The paper also records a PCIHF variant proposed with its
  modulus raised from $2^{160}$ to $2^{1600}$ in response to Wagner's algorithm.
- `maitin-shepard-et-al-2016`. Gives the binary-curve ECMH construction, its
  weak-encoding collision-security reduction and separate mapping,
  group-operation and batching measurements. It distinguishes a sufficient ECMH
  encoding from a hash indifferentiable from a random oracle; applying the
  reduction to this repository's adapted map still requires establishing that
  map's hypotheses and reduction loss.
- `boneh-lynn-shacham-2001`. MapToGroup, an early reference for counter-based
  rejection sampling into a curve group. Here the shared driver `curve::h2c`
  samples digest-derived candidate encodings until decoding succeeds; the
  weak-encoding condition of `maitin-shepard-et-al-2016`, a sufficient route
  to collision resistance, may be established from each codec's accepted
  preimages, including identity and torsion handling. The variable-time path
  is evaluated under the project's public-input assumption.
- `rfc9380`. Specifies hash-to-field, domain separation, mappings and cofactor
  handling for hash-to-curve suites. It provides a reference decomposition for
  the implementation and benchmarks, while its standardized suites do not
  directly certify arbitrary ephemeral curves or the binary maps used here.
- `bernstein-et-al-2013`. Introduces Elligator and the distinction between
  curve-point encodings and uniform-looking bit strings. Elligator 2 is an
  applicable mapping reference for suitable odd-characteristic curves with
  rational two-torsion, but its encoding alone must not be labeled a full
  random-oracle hash.

### Field primitives used by mapping and decoding

- `bernstein-2005`. Poly1305 evaluates a message as a polynomial in a clamped
  key over $\mathbb{F}_{2^{130} - 5}$, a field chosen for fast integer
  arithmetic, and bounds the differential probability by the polynomial's
  degree. The projection of item identifiers evaluates a two-block message in
  the same way, unclamped and without the pad; its field is measured beside
  $\mathrm{GF}(2^{127})$ and $\mathbb{F}_{2^{127} - 1}$, which are cheaper
  where their multiplications are ([`hash`](../src/hash.rs)).

## Field and group arithmetic, set reconciliation and baselines

### Set reconciliation and the cell checksum

- `yang-gilad-alizadeh-2024`. Rateless IBLTs: a coded-symbol stream whose
  communication overhead approaches 1.35 symbols per set difference, and the
  item-to-cell mapping `src/riblt.rs` implements. Each cell carries a
  checksum used to decide whether the cell holds exactly one item; the
  question of what that checksum must resist, raised in yangl1996/riblt#3,
  is the question this repository answers.
- `minsky-trachtenberg-zippel-2003`, `dodis-et-al-2008`, `naumenko-et-al-2019`.
  Algebraic set reconciliation: characteristic-polynomial interpolation, the
  PinSketch BCH syndrome, and its use in Erlay (Minisketch), proposed for
  Bitcoin's transaction relay and not deployed. Decoding is exact and needs no
  checksum, so the purity question does not arise, but decoding costs time
  quadratic in the difference and the capacity is fixed in advance. This is an
  alternative to rateless IBLTs; the case for RIBLTs rests on rateless,
  linear-time decoding.
- `riblt-ecmh`. An independent Go implementation of RIBLT with an
  elliptic-curve multiset-hash checksum: the nearest prior implementation
  of the combination evaluated here.
- `meyer-scherer-2024`. Range-based set reconciliation over conventional
  hashes. It names censorship resistance, that a malicious party cannot
  generate sets which two honest nodes fail to reconcile, and universality,
  that a node's main computation does not depend on its peer's set; the
  adversary of [Adversary](problem.md#adversary) is the one their first
  property resists. It observes that RIBLTs meet adversarial sets by
  randomizing the hash for each session, which repeats the hashing in every
  session; a salt shared by an epoch's sessions amortizes that work.
- `hoyte-2023`. Range-based reconciliation as deployed in negentropy. Its
  fingerprint threat model is the third party of
  [Adversary](problem.md#adversary), who inserts crafted elements so that others
  fail to synchronize. It reports collisions of 256-bit XOR fingerprints in two
  seconds by Gaussian elimination and of sums modulo $2^{256}$ in about 28 hours
  with a $k$-dimensional birthday solver, the kind of problem to which known
  discrete logarithms reduce the checksum, and judges ECMH collisions
  infeasible but ECMH slower, with a single reference implementation.
- `corallo-2016`. Compact block relay keys short transaction identifiers with
  SipHash under a hash of the block header and a nonce. The header commits to
  the transactions it identifies, so they cannot be chosen against the key.
  Items reconciled within an epoch can be, which is why a short keyed hash
  does not serve here although a block hash still bounds the epoch.
- `krohn-freedman-mazieres-2004`, `juels-et-al-2015`,
  `boneh-freeman-katz-waters-2009`. Integrity of coded data against corrupt
  symbols from participants: homomorphic hashes of a publisher's blocks verify
  each block of a rateless erasure code on receipt, Falcon codes authenticate
  LT-coded symbols with a MAC under a key the sender and receiver share, and
  signatures on a linear subspace let the nodes of a random linear network code
  derive and verify signatures on combinations. The first uses a
  discrete-logarithm homomorphic hash, as the cell checksum does. Each
  authenticates data against a source or a key from which every legitimate
  symbol derives; the items of a reconciled set are authored by any party, the
  adversary included, and carry no such credential, so corrupt symbols from a
  participating peer remain outside the adversary of
  [Adversary](problem.md#adversary).
- `bolt-7`, `chen-et-al-2026`. Lightning gossip, a public set of channel and
  node announcements that any node can extend, and an evaluation of IBLT-based
  reconciliation for it, which leaves adversarial conditions to future work.

### Multiset hashes in deployment, and fixed-group baselines

- `fips-180-4`, `aumasson-bernstein-2012`. SHA-256, the digest behind
  `src/hash.rs` and every map, and SipHash. XOR of either is the
  non-adversarial baseline every curve is compared against.

### Binary curve models and addition formulas

- `pornin-2022`. Complete formulas for ordinary binary curves with $|E| = 2r$ in
  $(x, s)$ coordinates, the extended $(X : S : Z : T)$ accumulators of
  `curve::binary::extended`, the $w$ encoding of `wcodec`, and the mixed
  addition that `unscaled` reduces to $7M + 2S$ by moving the curve constant
  into the addend, against $8M + 2S$ for the implemented incomplete $\lambda$
  formula.
- `pornin-2023`. The same formulas over
  $mathrm{GF}(2^{254}) = mathrm{GF}(2^{127})[u]$ (GLS254) and the map to the
  curve.
