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

- `pollard-1978`. The rho method: $\sqrt{\pi r/2}$ expected group operations and
  constant memory. The reference attack behind every "rho" figure in the report.
- `wiener-zuccherato-1998`, `duursma-gaudry-morain-1999`. Rho runs $\sqrt{2}$
  faster with the negation map and $\sqrt{2m}$ faster on curves with an
  automorphism group of order $2m$ (Koblitz curves, curves with efficiently
  computable endomorphisms). The report divides the reference work factor by
  these factors where they apply: by $\sqrt{2}$ for every family, and by a
  further $\sqrt{2}$ for the GLS family over $\mathrm{GF}(2^{122})$, whose
  endomorphism group has order 4.
- `bernstein-lange-2012`, `bernstein-lange-2013`. With a precomputed table of
  size about $N^{2/3}$, each further discrete logarithm in a group of size $N$
  costs about $N^{1/3}$ operations. A fixed curve would let an attacker amortize
  this across epochs. Deriving the curve from the salt prevents direct reuse of
  curve-specific tables; field-level arithmetic and work on the candidate space
  remain reusable. These two papers motivate a curve per namespace.

### Weil descent and index calculus over extension fields

- `gaudry-hess-smart-2002`. The GHS construction: the descent gives a
  hyperelliptic curve over $\mathrm{GF}(q)$ of genus about $2^{m-1}$, $m$ the
  "magic number" determined by the minimal polynomial of $\sqrt{b}$ under
  Frobenius. `sage/ghs.sage` computes $m$ for the binary families and the genus
  that results.
- `menezes-qu-2001`. For $n$ prime in $[160, 600]$ the GHS genus is too large
  for any attack. The degree 127 lies outside that range and is the exceptional
  small prime degree: 2 has order 7 modulo 127, so $x^{127} - 1$ splits into
  degree-7 factors over $\mathbb{F}_2$ and the minimal genus is $2^6$, so the
  descent must be examined for $\mathrm{GF}(2^{127})$.
- `galbraith-hess-smart-2002`. Isogenous curves have equal order and the
  attack transfers along isogenies, so the vulnerable set is a union of isogeny
  classes, not of curves. `sage/ghs122.sage` checks the class, and the $m = 127$
  statement counts classes.
- `hess-2003`. Generalizes GHS: with $\sqrt{b} = \gamma_1 \cdot \gamma_2$ for
  $\gamma_1, \gamma_2$ in the span of a degree-7 factor's kernel and
  $\mathbb{F}_2$, the descent still has genus 127 or 128. Over
  $\mathrm{GF}(2^{127})$ this reaches 4537 Frobenius classes, 65 of them of
  order $2 \cdot \text{prime}$ and so admissible to the certificate, a heuristic
  $2^{-52}$ of the admissible isogeny classes (Hess's own estimate
  $s \cdot q^{2d}/(q^{n/2} \cdot n)$ gives $2^{-52.3}$). Accepted without a
  check: across $T$ epochs the probability is at most $T$ times this heuristic
  figure, and a descent compromises only its own epoch.
- `hess-2005`. Survey chapter covering the above; cited where the problem
  statement states the Weil-descent conclusion without the derivation.
- `hankerson-karabina-menezes-2009`. Security analysis of GLS curves over
  $\mathrm{GF}(2^{2m})$: the GHS descent to $\mathrm{GF}(2^m)$ and the $n = 2$
  index calculus leave the discrete logarithm at about $2^m$, the rho cost. It
  is the reference for the $\mathrm{GF}(2^{122})$ family's security statement.
- `gaudry-2009`. Index calculus for $E$ over $\mathrm{GF}(q^n)$ in
  $\widetilde{O}(q^{2-2/n})$ for fixed $n \ge 2$: $\widetilde{O}(q)$ at $n = 2$,
  which ties rho for every quadratic-extension family here (fp61x2, fp64x2,
  goldilocks2, and $\mathrm{GF}(2^{122})$ as $\mathrm{GF}(2^{61})^2$). The
  reason degree 2 is admissible and degree $\ge 3$ is not.
- `diem-2011`. The discrete logarithm over $\mathrm{GF}(q^n)$ is subexponential
  when $n$ grows with $\log q$ within fixed bounds. It does not apply at $n = 2$
  or at prime $n = 127, 109, 107$; cited to say so.
- `joux-vitse-2012`. Cover and decomposition attacks made practical for
  $\mathrm{GF}(q^6)$-type fields; the boundary case showing where
  composite-degree extensions become dangerous, and that $n = 2$ is not among
  them.

### Point counting and curve generation

- `schoof-1995`. Schoof's algorithm and the Elkies–Atkin improvements (SEA):
  polynomial-time point counting over any finite field. PARI's `ellsea` runs it
  for the odd-field families, with the Elkies-prime early abort the `pari`
  module uses.
- `satoh-2000`, `mestre-2000`, `gaudry-2002`. Canonical lifts and Mestre's AGM
  for characteristic 2: counting a curve over $\mathrm{GF}(2^{127})$ in
  milliseconds. The Rust `agm` counter is a canonical lift of this kind and
  PARI's `F2xq_ellcard` is Harley's variant of the AGM; `gaudry-2002` is the
  comparison that established the method.
- `satoh-skjernaa-taguchi-2003`. The SST variant of the canonical lift, one
  of the methods `gaudry-2002` compares; the method implemented is the AGM.
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
- `pari-gp`, `sagemath`. The point-counting and primality oracle for
  fixtures and the `pari` feature, and the system in which the reference
  programs under `sage/` are written; Sage counts points through PARI.

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
  them, the checksum is AdHash modulo $r$ and the algorithm applies ([Known
  weaknesses](ecc_security.md#known-weaknesses)).
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
  discrete logarithms reduce the checksum ([Known
  weaknesses](ecc_security.md#known-weaknesses)), and judges ECMH collisions
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

- `certicom-2010`, `libsecp256k1`. secp256k1's parameters and the library
  timed in `benches/compare.rs` (compressed-point parsing, ElligatorSwift
  decoding, Jacobian accumulation): the cost of ECMH on a fixed 256-bit
  curve.
- `bernstein-2006`, `hamburg-2015`, `rfc9496`, `curve25519-dalek`. Curve25519,
  the Decaf quotient that removes its cofactor, the ristretto255 specification,
  and the implementation timed in `benches/compare.rs`: the second fixed 256-bit
  baseline. Decaf's quotient, which removes a cofactor of 4, is the model for
  the `twisted` family's encoding of $E/\langle T \rangle$, which removes only
  the 2-torsion point $T$.
- `fips-180-4`, `aumasson-bernstein-2012`. SHA-256, the digest behind
  `src/hash.rs` and every map, and SipHash. XOR of either is the
  non-adversarial baseline every curve is compared against.

### Binary fields and their dependence on carry-less multiplication

- `crrl`. Pornin's library: the $\mathrm{GF}(2^{127})$ backends
  (`field::gf2_127`) and the GLS254 code from which `curve::binary127::map` and
  the $\mathrm{GF}(2^{122})$ tower are adapted.

### Binary curve models and addition formulas

- `pornin-2022`. Complete formulas for ordinary binary curves with $|E| = 2r$ in
  $(x, s)$ coordinates, the extended $(X : S : Z : T)$ accumulators of
  `curve::binary::extended`, the $w$ encoding of `wcodec`, and the mixed
  addition that `unscaled` reduces to $7M + 2S$ by moving the curve constant
  into the addend, against $8M + 2S$ for the implemented incomplete $\lambda$
  formula.
- `pornin-2023`. The same formulas over
  $\mathrm{GF}(2^{254}) = \mathrm{GF}(2^{127})[u]$ (GLS254) and the map to the
  curve. The $\mathrm{GF}(2^{122}) = \mathrm{GF}(2^{61})[u]$ tower and the GLS
  family (`binary.122-gls`) test whether this shape transfers at half the size.
- `pornin-2024`. The general method: a prime-order group abstraction with
  complete formulas from any curve of even order, most efficient for
  $|E| \equiv 2 \pmod{4}$. Both characteristics are covered, so it states once
  the structure the binary families ($E[r] = 2E$) and the `twisted` quotient
  use.
- `oliveira-et-al-2014`, `oliveira-et-al-2013`. $\lambda$-projective
  coordinates, with mixed addition $8M + 2S$ and no curve constant: the
  `curve::binary::lambda` accumulators. On curves with a dense $B$, where
  Pornin's extended formulas need two full multiplications by $\beta = B^{1/4}$,
  this costs less than the extended accumulator, at the price of incompleteness;
  the unscaled formulas, at $7M + 2S$, move $\beta$ into the addend and cost
  less than $\lambda$'s $8M + 2S$ for the addition itself.
- `galbraith-lin-scott-2009`. GLS curves: over $\mathrm{GF}(q^2)$, quadratic
  twists of curves defined over $\mathrm{GF}(q)$, with an endomorphism $\psi$,
  $\psi^2 = -1$. For binary curves, $a = u$ and $B$ in $\mathrm{GF}(2^{61})$, as
  in `binary.122-gls`. The endomorphism speeds scalar multiplication, which this
  workload does not perform; what the family gains here is a curve constant in
  $\mathrm{GF}(2^{61})$, so each multiplication by it costs two word products
  instead of three. Its rho penalty is in the security section.
- `bernstein-lange-farashahi-2008`. Binary Edwards curves: the complete model in
  characteristic 2. For $n \ge 3$ every ordinary curve over $\mathrm{GF}(2^n)$
  is birationally equivalent to a complete binary Edwards curve (Theorem 4.3),
  so the model is available for every certified binary curve and the exclusion
  is one of cost. Mixed addition costs $13M + 3S + 3D$ (projective
  $21M + 1S + 4D$; $16M + 1S + 4D$ for $d_1 = d_2$, with no mixed variant of
  lower cost) against $7M + 2S$ for the unscaled formulas of `pornin-2022`. No
  direct map to the model with a demonstrated advantage was found, and generic
  projective equality costs about $4M$ against $1M$ to $3M$ on the $(x, s)$
  representations; neither is a lower bound for a cached or quotient
  representation. Excluded on these inspected operation counts under equal field
  arithmetic; not implemented or measured, so the exclusion is a prioritization,
  not a lower bound. The $d_1 = d_2$ forms, of lower cost, have
  $\mathrm{Tr}(d_1^2 + d_2) = 0$ in the birational Weierstrass model, and so
  $4 \mid |E|$; the certified binary curves, with $|E| = 2r$ and $r$ odd,
  require $d_1 \ne d_2$.
- `kim-lee-negre-2014`. Revisits binary Edwards arithmetic and keeps
  completeness on complete curves ($\mathrm{Tr}(d_2) = 1$). For general
  $d_1, d_2$, projective addition costs $15M + 2S + 4D$ (the paper's baseline
  for the 2008 formulas, from the EFD, is $18M + 3S + 6D$), mixed addition
  $13M + 2S + 2D$ and doubling $2M + 5S + 3D$. The $d_1 = d_2$ forms (projective
  addition $14M + 1S + 2D$, mixed $13M + 1S + 2D$) and the complete differential
  addition-and-doubling at $5M + 4S + D$ exist exactly when the curve has a
  rational point of order 4 (Lemma 1), so not on the certified curves; the
  differential formulas also assume a fixed difference $\omega_0 \ne 0$ (and
  $\omega_0 \ne 1$ when $d_1 = d_2$). Mixed addition stays near twice the
  $7M + 2S$ of `pornin-2022` with dense constants, and the differential result
  does not implement an arbitrary cell update.
- `bernstein-2009`, `koziel-azarderakhsh-mozaffari-2015`, `koziel-bec-small`,
  `farias-albertini-barreto-2018`, `loiseau-fournier-2018`,
  `hajra-karati-sen-2026`. The fast binary Edwards formulas are differential, in
  the coordinate $w = x + y$, which identifies $P$ with $-P$: $w(Q + P)$ is
  computed from $w(P)$, $w(Q)$ and $w(Q - P)$, as a Montgomery ladder supplies
  it. Batch Binary Edwards gives the $d_1 = d_2$ mixed form at $5M + 2D + 4S$
  with $d$ chosen sparse so that multiplication by $d$ reduces to a few shifts
  and additions (`bernstein-2009`, section "Differential addition and
  doubling"); the later papers correct and re-cost these $w$-formulas
  ($5M + 1D + 4S$ mixed, `koziel-azarderakhsh-mozaffari-2015`;
  `koziel-bec-small` is its Sage check), apply them with $d_1 = d_2$
  (`farias-albertini-barreto-2018`, Section 4), add the co-Z technique
  (`loiseau-fournier-2018`, Algorithm 2), and vectorize the ladder step with
  VPCLMULQDQ on a $d_1 = d_2$ curve with sparse $d$ and cofactor 4
  (`hajra-karati-sen-2026`, Table 2). A cell update adds an independent hashed
  point with no known difference, so none of these formulas applies to it; for
  generic dense parameters $D$ counts as a full $M$, and the sparse-parameter
  specialization they rely on is not established for seed-derived curves.

### Odd-characteristic curve models and addition formulas

- `edwards-2007`, `bernstein-lange-2007`. The Edwards normal form and its
  complete addition law when $d$ is a non-square: the `edwards` family, with
  $d$ drawn per salt.
- `bernstein-et-al-2008`, `hisil-et-al-2008`. Twisted Edwards curves and
  extended coordinates. With $a = -1$ a square, the addition is complete and
  a cached addend costs $7M$; this is the `twisted` family's accumulator, and
  the reason it requires $q \equiv 1 \pmod{4}$.
- `renes-costello-batina-2016`. Complete projective formulas for
  prime-order short Weierstrass curves, $11M + 2 m_b$ for mixed addition:
  the `weier` family's accumulator, and the cost of completeness when the
  order has no cofactor to quotient away.
- `efd`. The Explicit-Formulas Database: the named formulas
  (madd-2007-bl, dbl-2001-b) in `curve::weier::jacobian` and the operation
  counts quoted in the module documentation.
- `montgomery-1987`. Montgomery curves and simultaneous inversion. The inversion
  method is `field::batch`, which every batch sum and batch normalization uses;
  the Montgomery model is the form in which `edwards` points travel and are
  summed in batches by the affine chord law.

### Prime fields

- `polygon-zero-2022`, `plonky3`. The Goldilocks prime $2^{64} - 2^{32} + 1$ and
  the implementation `field::goldilocks2` uses for its base field and
  binomial extension. Its 2-adicity of 32, chosen for FFTs, makes square
  roots slower than over $2^{64} - 59$; the `goldilocks2` rows measure that
  cost.
- `bruestle-gafni-2023`, `haboeck-levit-papini-2024`. BabyBear and the Mersenne
  prime $2^{31} - 1$, the 31-bit fields of current STARK systems. Reaching
  $2^{122}$ points over a 31-bit prime needs an extension of degree 4, which
  admits the heuristic index-calculus asymptotic $\widetilde{O}(q^{3/2})$ of
  `gaudry-2009`, below the generic $q^2$. Group size alone therefore does not
  justify the target security, and these fields are excluded from the
  comparison.
