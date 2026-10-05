---
type: Chapter
title: Problem statement
description: The adversary, the beacon-salted setting and its reference work factor, per-namespace curves against precomputation, curve selection and certificates, design criteria and scope.
tags: [threat-model, precomputation, curve-generation, point-counting, certificates, weil-descent]
sources:
  - id: bibliography
    resource: ../references.bib
    title: Bibliography, annotated in literature.md
generated: { by: claude-code/claude-opus-5-5 }
contributors: [codex/gpt-6, claude-code/claude-fable-5-1]
status: draft
---

# Problem statement

## Adversary

The sets considered accept items from any party, and the peers that hold them
reconcile with one another. The mempools of Bitcoin nodes, for whose transaction
relay Erlay proposes set reconciliation [naumenko-et-al-2019], and the channel
and node announcements of Lightning gossip [bolt-7], for which IBLT-based
reconciliation has been evaluated [chen-et-al-2026], are such sets: their items
are adversarial inputs. Neither proposal is deployed; both networks relay by
flooding. Meyer and Scherer call the property required of reconciliation
censorship resistance: a malicious party cannot generate sets that two honest
nodes fail to reconcile [meyer-scherer-2024]. Hoyte states the same adversary
for the fingerprints of range-based reconciliation: a third party inserts
crafted elements so that innocent ones fail to synchronize [hoyte-2023].

The adversary is a third party. It can evaluate the item hash once the salt is
published, and it can place items of its choosing in the sets of honest peers,
for instance by broadcasting them. It authors items whose hashes satisfy a
relation, aiming to make two honest peers fail to reconcile, or obtain a wrong
difference, when the items lie in the symmetric difference of their sets; items
common to both sets cancel. The RIBLT's mapping of items to coded symbols admits
a similar attack that does not involve the checksum ([Known
weaknesses](ecc_security.md#known-weaknesses)).

The items, the salt and the curve are public, so the timing of the computation
reveals nothing that the adversary could not compute itself, and implementations
may run in variable time. Hashing may therefore use try-and-increment
[boneh-lynn-shacham-2001], whose number of attempts depends on the item. RFC
9380 specifies deterministic maps instead and recommends against rejection
sampling because its timing has leaked secret inputs [rfc9380, Section 1]; it
requires constant-time hashing when the inputs to an encoding are secret
[rfc9380, Section 10.3], which they are not here.

The item hash assumes secure hash functions of three kinds, which differ in what
they bind and in how long their outputs live. An item is identified by a 256-bit
digest that does not depend on the salt: tagged SHA-256 of its encoding, or an
identifier the application already holds, such as a transaction's txid. The
identifier is computed once and serves every namespace. The adversary chooses
it, through the preimage, and only its collision resistance, about $2^{128}$,
binds an item to it.

Hashing to the curve and the RIBLT's mapping depend on the salt, although the
curve already changes with the beacon value. Each map starts from a field
element and a sign, taken from at most 128 bits of a digest, before the curve
enters, and two items whose elements and signs agree hash to the same point on
every curve that accepts that element. Without the salt, such a pair would cost
about $2^{64}$ digests once and serve in every namespace; with it, the search is
repeated for each namespace. An unsalted mapping would likewise let colliding
schedules be searched for before the beacon ([Known
weaknesses](ecc_security.md#known-weaknesses)). The two keys serve different
ends. The hash to the curve keeps the public salt, under which an item's point
identifies it, and resists collisions, in every session of the namespace. The
schedule is better keyed by a secret that each pair of peers shares, where one
is available: no third party can then search for colliding schedules, and a peer
whose neighbours have converged receives independent streams from them rather
than one stream repeated. The benchmarks derive both keys from one salt, which
changes no cost.

Neither use needs a salted 256-bit digest of each item. Once per namespace,
salted SHA-256, modelled as a random oracle, derives keys from the beacon value,
or the schedule's from a secret, under a tag for each use; once per item and
namespace, a universal hash under these keys projects the identifier to 128 bits
[wegman-carter-1981]. With the identifier read as two elements $m_1, m_2$ of a
field of about $2^{127}$ elements and a nonzero key $k$, the projection
$U_k(m_1, m_2) = m_1 k^2 + m_2 k$ is Poly1305's polynomial on a two-block
message, without its clamping or one-time pad [bernstein-2005]. Two identifiers
fixed before the salt collide for at most two keys, with probability about
$2^{-126}$. Once the salt is published, the adversary can evaluate the
projection, but its inputs are SHA-256 outputs, which it does not choose: a
random oracle followed by a fixed function whose fibres are of nearly equal size
is a random oracle into the field, and a colliding pair costs a generic search
over about $2^{64}$ identifiers, as for a salted 128-bit digest. The
projection's linearity gives no shortcut. Items are mapped to the curve
independently and added only as points, so a relation among projections is a
relation among points only through a symmetry of the map, such as
$\mathrm{Map}(-u) = -\mathrm{Map}(u)$ for simplified SWU and Elligator 2
[rfc9380, bernstein-et-al-2013]; to exploit one before the salt, two identifiers
must agree up to sign as field elements in both halves, about 254 bits. The
field is therefore chosen by cost: $\mathrm{GF}(2^{127})$ where the target has
carry-less multiplication, $\mathbb{F}_{2^{127} - 1}$ where it does not, with
Poly1305's $\mathbb{F}_{2^{130} - 5}$ measured as a comparison
([`hash`](../src/hash.rs)). The attempts of a try-and-increment map after the
first use further keys over the same identifier, so their candidates are linear
in the first attempt's rather than independent; the deterministic maps, Pornin's
[pornin-2023], Elligator 2 and simplified SWU, use only the first. For the
mapping, 128 bits suffice as well: two items stall a rateless decoder
indefinitely only if their schedules agree entirely, which after the salt
requires equal seeds, about $2^{64}$ identifiers. Schedules that agree on a
prefix only delay decoding, at a cost that does not depend on the seed's width;
that is a property of the sketch, not of its checksum. The implementation
digests items with salted SHA-256 for both uses by default.
[`Riblt::projected`](../src/riblt.rs) keys cells by the items' identifiers and
replaces both digests by projections, each under keys of its own and in a field
chosen for its cost; both are benched ([Costs by
lifetime](methodology.md#costs-by-lifetime)).

What need not be uniform on the group, or indifferentiable from a random oracle,
is the map from digests to points: relations must be hard to find within the
subset of the group that it reaches. Maitin-Shepard et al. give a sufficient
condition. Let the group have a direct factor of prime order $r$ whose
complement can be sampled efficiently, and let the output of a random oracle be
mapped to the group by an $(\alpha, \beta)$-weak encoding: one whose preimages
can be counted efficiently and sampled uniformly, of which no point has more
than $\alpha$, and of which a uniformly random point has on average at least
$\alpha/\beta$. Then finding a multiset collision with multiplicities below $r$
is as hard as computing discrete logarithms, at a cost of simulating the oracle
of about $\beta$ scalar multiplications per query and half the success
probability [maitin-shepard-et-al-2016, Section 4.1, Theorem 1, and Appendix A].
An encoding that misses a fixed fraction of the points can qualify: their
characteristic-2 Shallue–van de Woestijne encoding misses about $9/32$ of them
and has $\beta$ close to 3 [maitin-shepard-et-al-2016, Section 4.2]. Whether the
maps used here meet these hypotheses is open ([Known
weaknesses](ecc_security.md#known-weaknesses)).

A participating peer can always prevent reconciliation by misrepresenting its
set; no checksum prevents this, and it is not considered. Corrupt coded symbols
are the pollution attack of network coding, in which forwarded packets are
linear combinations of a source's packets, and there the source authenticates
the space that the combinations must lie in. Homomorphic hashes of a publisher's
blocks let a receiver verify each encoded block of a rateless erasure code as it
arrives [krohn-freedman-mazieres-2004], authenticated LT codes such as Falcon
codes verify each symbol with a message authentication code under a key shared
by sender and receiver [juels-et-al-2015], and signatures on the linear subspace
spanned by a source's packets let any node derive and verify signatures on their
combinations in random linear network coding [boneh-freeman-katz-waters-2009]. A
coded symbol of a RIBLT is likewise a sum, and its checksum a homomorphic hash
of the items it holds, but none of these schemes applies: each authenticates
data against a source or a key from which every legitimate symbol derives,
whereas the items of a reconciled set are authored by any party, the adversary
included, and carry no signature or key that distinguishes an honest item from
an adversarial one. Whoever supplies a curve certificate stands in a separate
relation to its verifier ([Parameter selection and
verification](#parameter-selection-and-verification)).

## Setting and reference work factor

In the universal setting the item hash is public, and the checksum must resist
relation finding without a secret. Keying the hash with a secret the peers
coordinate, as the Rateless IBLT paper proposes [yang-gilad-alizadeh-2024,
Section 4.3], defeats the precomputation of relations among item hashes but
repeats the hashing of every item for each key; randomized for each
reconciliation session, it does so in every session [meyer-scherer-2024, Section
II]. This repository considers namespaces whose salt is published by a
randomness beacon, such as a block hash [bonneau-clark-goldfeder-2015], and
whose checksums are shared by all sessions within a bounded epoch. Bitcoin's
compact block relay keys the short identifiers of a block's transactions with a
hash of the block header and a nonce [corallo-2016]. The header commits to those
transactions, so they cannot have been chosen in relation to the key: as in the
keyed model, the set is fixed before the key is revealed, and a short keyed hash
serves. Items reconciled within an epoch carry no such guarantee: the adversary
may author them after the salt is published, so the checksum must resist
relation finding for the duration of the epoch, which the regular publication of
new blocks bounds. A forgery affects reconciliation within one namespace only.
This study assumes that the protocol validates the items it obtains
independently of reconciliation, as Bitcoin nodes validate transactions, so that
a failed or wrong decoding costs a repeated reconciliation and admits no invalid
item. A failure of soundness for one curve, even a discrete logarithm computed
in its group, is then a failure of liveness, and a transient one: the next
beacon value selects a new curve. The short horizon and the limited scope of
failure motivate a reference work factor between $2^{48}$ and $2^{64}$ group
operations, with the upper part of the range preferred; the families compared
here lie between $2^{52}$ and $2^{63}$ ([Security
considerations](ecc_security.md)). This range is an application assumption under
evaluation; salting does not by itself imply it.

Maitin-Shepard et al. [maitin-shepard-et-al-2016, Section 4.1] reduce collision
finding with bounded multiplicities to the discrete-logarithm problem in a large
prime-order component of the checksum group, under the random-oracle and
encoding assumptions stated there. For a prime component of order $r$,
$\sqrt{r}$ group operations serve here as the generic reference cost. An order
near $2^{125}$ gives a reference cost near $2^{62.5}$. Several curve models over
fields of about 127 bits admit 16-byte encodings of such groups. The security
assessment must also cover the hash maps and parameter selection actually used.

## Precomputation and per-namespace curves

A public salt over a fixed checksum group serves here as a point of comparison,
not as a design under consideration: the salt changes the item-to-point map, but
the group, and any precomputation against it, serves every epoch. Bernstein and
Lange [bernstein-lange-2012] give a discrete-logarithm algorithm that, after
about $r^{2/3}$ operations of precomputation, uses a table of about $r^{1/3}$
entries and about $r^{1/3}$ group operations per logarithm;
[bernstein-lange-2013] discusses the cost model for such precomputation. For $r$
near $2^{126}$ these scales are $2^{84}$, $2^{42}$, and $2^{42}$ respectively.
They are leading-order estimates for a time–memory tradeoff, not measured costs
and not the cost of a RIBLT forgery. Their significance depends on offline
resources, storage, and the number of targets.

Each beacon value therefore determines both the salt and the curve. Subject to
the entropy of the beacon and to the parameter-selection procedure, this
prevents direct reuse of tables computed for earlier curves. It does not exclude
all useful precomputation, nor does it guarantee the generic work factor in
every namespace.

## Parameter selection and verification

Candidate curves are derived from the beacon value by hashing, and a namespace
uses the first candidate whose group order is an admissible cofactor times a
prime $r$. Deriving curves from public randomness has precedent in verifiably
random seeds [ansi-x9-62-1998, rfc5639], lottery draws [baigneres-et-al-2015],
and a beacon whose outputs seed a stream of curves [lenstra-wesolowski-2017];
those curves serve general use at the 128-bit level, whereas here each beacon
value selects a smaller one for its namespace. A participant may determine that
curve by counting points itself, or may check a certificate: for the accepted
candidate, $r$ and a hashed point $P$ whose cofactor multiple $Q$ is not $O$ and
satisfies $rQ = O$, and for each earlier candidate a witness that its order is
inadmissible, as trx proposes for its own enumeration. The verifier accepts any
witness that satisfies the certificate rules. This repository's prover derives
its witnesses from hashed points, so its certificate is determined by the seed
and the rejection policy (how labels are chosen, and whether by factoring).
Checking a certificate requires scalar multiplications, a probable-prime test,
and an embedding-degree bound, but no point counting.

The probable-prime test is Miller–Rabin to 24 fixed bases, so the primality of
$r$ is probable, not proven. A candidate with the required prime-order structure
has no valid rejection witness, so a certificate cannot skip it. The prover
supplies $r$ for the accepted candidate; the verifier checks its admissible
interval, the condition $rQ = O$ with $Q$ not $O$, and the 24 Miller–Rabin
bases. The second condition establishes that the order of $Q$ divides $r$. If
$r$ is prime, $Q$ has order $r$, and the family's cofactor and the Hasse bound
then establish $|E| = hr$. The fixed-base test does not prove that $r$ is prime.

The direct checker and certificate verifier share the order-parameter
predicates: the probable-prime test, the Hasse interval at the required
cofactor, exclusion of $r = q$, and the embedding-degree bound. Direct checking
takes the full order from a trusted point counter. Certificate verification
uses the nonidentity point and Hasse argument above instead; its order
conclusion remains conditional on $r$ being prime.

These formats do not certify every possible outcome of direct checking. An
order-admissible candidate can fail the fixed certificate-point check, and the
formats have no rejection witness for an anomalous order or a failed embedding
bound. A family without a full-order rejection witness also needs a factor for
a rejected composite. Direct selection can proceed without those witnesses;
the prover reports failure rather than skipping an uncertifiable candidate.
Whenever certification succeeds, its selected index agrees with direct
selection under the same candidate criteria and filters that reject only
candidates failing those criteria.

The families differ in the cost of point counting. Two implementations count
points here. PARI/GP [pari-gp] is a computer algebra system for number theory: a
C library, PARI, with an interpreter, GP. The optional `pari` feature links the
library and calls it for point counting and integer factorization
([`curvegen::pari`](../src/curvegen/pari.rs)). SageMath (Sage) [sagemath], in
which the reference programs are written, also counts points through PARI. The
other implementation is this repository's own, in Rust, for binary fields only;
it has no dependency on PARI.

Over binary fields, Satoh's canonical lift [satoh-2000] and Mestre's
arithmetic–geometric mean, in the form compared and combined by Gaudry
[gaudry-2002], count points in polynomial time using 2-adic arithmetic alone.
The repository implements this for $\mathrm{GF}(2^{127})$,
$\mathrm{GF}(2^{109})$ and $\mathrm{GF}(2^{122})$
([`curvegen::agm`](../src/curvegen/agm/mod.rs)), counting the dense 122-bit
curves in a pentanomial basis and the Galbraith–Lin–Scott (GLS) curves
[galbraith-lin-scott-2009] through their subfield curve over
$\mathrm{GF}(2^{61})$. A GLS constant is $B = \beta^4$ for $\beta$ in
$\mathrm{GF}(2^{61})$, so GLS candidates are drawn from about $2^{61}$
constants, against about $2^{122}$ for the dense family; the accepted curves are
a subset of each. For all three fields the Rust prover produces certificates
without PARI: a candidate whose odd part is composite with no prime factor below
the trial-division bound ($2^{16}$ for $\mathrm{GF}(2^{122})$, $2^{10}$
otherwise) is rejected by its full order rather than by a factor.

The known-answer certificates were produced with Sage, which derives its
witnesses as the Rust prover does, from hashed points
([`sage/kat_common.sage`](../sage/kat_common.sage)). For the 109- and 122-bit
families it labels a rejection by the smallest prime factor of its odd part
(below $2^{16}$ in the 122-bit scripts, without a bound in the 109-bit one).
Without PARI, the Rust prover reproduces the indices and orders of these
certificates and the 122-bit certificates byte for byte; with PARI supplying the
factors, it reproduces the 109-bit certificates byte for byte as well.

Over prime fields, the corresponding method is the Schoof–Elkies–Atkin algorithm
[schoof-1995], which requires modular polynomials. The repository does not
implement it: prime-field and $\mathrm{GF}(p^2)$ curves are counted by PARI,
which handles curves with complex multiplication by a small discriminant
separately before falling back to it. A participant in a prime-field namespace
who does not run PARI or an equivalent system must therefore rely on another
party's certificate to establish which curve is in use, whereas a binary curve
can be derived and checked independently by a small, self-contained program. In
a setting where every participant derives the curve from a public beacon, this
is an advantage of binary curves, separate from their arithmetic cost.

The cost of selection is paid once per namespace. Each run report's Curve
selection section gives, per family and seed of `curvegen.csv`, the wall-clock
time to find a curve from a seed, to find and certify it, and to verify the
certificate, and how the prover's time divides among the phases of the search.
The number of candidates before an admissible order is geometric, so the time to
find varies severalfold between seeds. On an Apple M4 the Rust implementation
finds a $\mathrm{GF}(2^{127})$ curve in about 20 ms on average, an order of
magnitude faster than PARI on the same seeds; PARI finds the prime-field and
$\mathrm{GF}(p^2)$ curves in 0.25 to 2 s on average; and certificates verify in
under a millisecond on average.

Binary fields add a structural check that prime fields do not need: Weil descent
through subfields. For the prime extension degrees 127 and 109 the only proper
subfield is $\mathbb{F}_2$ [menezes-qu-2001]. For $m = 109$, 2 has order 36
modulo 109 and every descent has genus at least $2^{35}$. For $m = 127$, 2 has
order 7, and the generalized Gaudry–Hess–Smart (GHS) descent
[gaudry-hess-smart-2002, hess-2003], which extends to isogenous curves
[galbraith-hess-smart-2002], reaches curves of admissible order with heuristic
probability about $2^{-52}$ per epoch: 65 isogeny classes of order $2r$, against
none among the 36 classes the basic descent reaches
([`sage/ghs.sage`](../sage/ghs.sage)). No check is made; a descent compromises
only the epoch whose curve it reaches. The composite degree 122 is assessed in
[`curvegen::select122`](../src/curvegen/select122.rs).

## Design criteria

These assumptions differ from those of fixed-curve design:

- The operations that dominate the cost are one hash-to-curve per item, followed
  by one point addition or subtraction per coded symbol the item maps to, plus
  encoding and equality tests during peeling. Scalar multiplication is not among
  them; only curve selection and certification use it.
- Curve parameters vary by namespace. An optimization tied to particular
  parameters must include its setup cost and its amortization within one
  namespace. Constants fixed by the field or the curve model remain available.
- The arithmetic processes public data, so variable-time implementations are
  acceptable ([Adversary](#adversary)).
- Point counting and parameter certification occur once per namespace and are
  amortized, but must remain tractable for each new beacon value ([Parameter
  selection and verification](#parameter-selection-and-verification)).
- Field and curve structure require separate analysis for attacks below the
  reference work factor, including applicable Weil-descent methods
  [gaudry-hess-smart-2002, hess-2005, gaudry-2009, diem-2011, joux-vitse-2012].
  The existence of a descent does not by itself establish such an attack.

[Candidate constructions](constructions.md) lists the constructions compared, at
128-bit and shorter encodings. Their arithmetic costs are measured here; their
security conditions are assessed separately.

## Questions outside this study

Which beacon value salts which data, and the interval in which a participant can
withhold or delay a value it has observed, are questions for the protocol that
uses the checksum; they are not addressed here. Whether the checksum covers sets
or multisets changes the bounds on relation finding of [Setting and reference
work factor](#setting-and-reference-work-factor) but not the constructions
compared. Targets without carry-less multiplication, WebAssembly among them, are
represented only by the portable build described in [Measurement
methodology](methodology.md). How the RIBLT performs when a peer reconciles with
several others, and how far a third party can delay decoding with items whose
schedules agree, are questions for the sketch rather than for its checksum.
