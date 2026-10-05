# Ephemeral ECMH

## Authorship and disclaimer

The author formulated the research question, defined the threat model
([Adversary](docs/problem.md#adversary)) and security arguments not related to
curve selection, and directed the language models that carried out the work.
The implementation and evaluation were carried out almost entirely by Claude
Opus 5.5, with review and contributions from Codex Astra 6 and Claude Fable
5.1.

The author is not an expert in elliptic curve cryptography and has only a
vague familiarity with binary elliptic curves. This report was written for the
author as much as under his direction: the models assembled most of the
security arguments for curve selection and the literature they rest on, with
most of the cited works being new to the author. It is shared before the author
is able to verify it fully, and has not been reviewed by a specialist.

## Synopsis

This repository evaluates elliptic-curve multiset hashes (ECMH)
[maitin-shepard-et-al-2016] over small curves as the cell checksum of rateless
invertible Bloom lookup tables (RIBLTs) [yang-gilad-alizadeh-2024].

When the sets being reconciled contain potentially adversarial inputs, RIBLT
can be susceptible to collisions. This can be addressed using a secret shared
by the peers engaging in reconciliation, or by utilizing ECMH.

The aim is to relax the requirement for a shared secret and mitigate
collisions caused by adversarial inputs, using only a public beacon value as a
salt, but without paying the full cost of ECMH on standard curves.

Weaker curves can be randomly selected periodically based on the public
randomness beacon, and this project investigates the suitability and
performance of such curves (and their underlying fields) in the context of
RIBLT over sets with potentially adversarial items. The question is whether such
a checksum attains a work factor between $2^{48}$ and $2^{64}$ group operations
for each curve, and which field, curve model, and representation minimize its
cost per item.

## Motivation

Set reconciliation lets two peers that hold similar sets exchange only their
difference. In peer-to-peer networks the sets often accept items from any party:
Erlay proposes set reconciliation for the transaction relay of Bitcoin nodes
[naumenko-et-al-2019], and IBLT-based reconciliation has been evaluated for the
channel and node announcements of Lightning gossip [bolt-7, chen-et-al-2026];
neither is deployed, and both networks relay by flooding. Their items may be
adversarial inputs. The adversary considered here is a third party that authors
items and places them in the sets of honest peers, for instance by broadcasting
them, aiming to make two honest peers fail to reconcile or obtain a wrong
difference; resistance to it is what Meyer and Scherer call censorship
resistance [meyer-scherer-2024]. A participating peer can always disrupt
reconciliation by sending corrupt data or misrepresenting its set; that is
unavoidable and not considered ([Adversary](docs/problem.md#adversary)).

A RIBLT [yang-gilad-alizadeh-2024] recognizes a coded symbol that holds a single
item by a checksum. In the reference implementation the checksum is the XOR of
item hashes. An adversary who can evaluate the item hash constructs two distinct
sets with equal checksums by linear algebra over $\mathbb{F}_2$, at negligible
cost and with any hash function
([yangl1996/riblt#3](https://github.com/yangl1996/riblt/issues/3)). Addition
modulo $2^{64}$ in place of XOR has the same weakness, by lattice reduction and
Wagner's algorithm [wagner-2002].

A secret key defeats the attack when the adversary cannot evaluate the keyed
hash. This is the model of short keyed hashes such as SipHash
[aumasson-bernstein-2012], designed against hash flooding of hash tables.
Bitcoin's compact block relay keys the 48-bit short identifiers of a block's
transactions with a hash of the block header and a nonce [corallo-2016]: the
header commits to the transactions, so the set is fixed before the key is
revealed and no transaction can be chosen against it. Open sets have no such
ordering. Items keep arriving after a key shared by all peers is published, and
the attack on the XOR checksum applies whatever the hash. A key kept per
connection, as Erlay's salts are [naumenko-et-al-2019], or a secret key
coordinated by the two peers, as the Rateless IBLT paper proposes against
injected items, repeats the hashing of every item for each key, and the
checksums of coded symbols no longer serve every peer, though their sums do
[yang-gilad-alizadeh-2024, Section 4.3]. Meyer and Scherer count that cost
against randomizing per session [meyer-scherer-2024].

The hashing can instead be shared within an epoch. In a namespace salted by a
beacon value, such as a key derived from a block hash, an adversary who can
neither predict nor influence that value cannot search for collisions before the
namespace begins. Unlike in compact block relay, it can still author items
afterwards, so the checksum must resist relation finding without a secret until
the next beacon value; with a block hash as the beacon, new blocks determine how
long that is.

An ECMH sums hashed points in an elliptic-curve group. For maps to the curve
that satisfy the hypotheses of the known reductions, its collision resistance
reduces to the discrete-logarithm problem in that group
[maitin-shepard-et-al-2016]; whether the maps used here satisfy them is open
([Known weaknesses](docs/ecc_security.md#known-weaknesses)). The curves in use
for this purpose have about $2^{256}$ points and about 128-bit security, and
hashing an item to them costs about 110 to 195 times the SHA-256 hash of the XOR
checksum (see [Findings](#findings)). A shorter horizon and a failure confined
to one namespace suggest that a work factor of $2^{48}$ to $2^{64}$ suffices,
provided the application detects a sketch collision by other means and recovers,
for example by re-keying or by another reconciliation method. This is an
application assumption.

A small fixed curve is exposed to precomputation [bernstein-lange-2012]. With
the discrete logarithms of the item hashes, the checksum becomes AdHash modulo
the group order, and the attack on the XOR checksum carries over through
Wagner's algorithm or lattice reduction; without them, the group does not expose
the residues modulo its order that these algorithms need. Hoyte reports such
collisions modulo $2^{256}$ in about 28 hours on eight cores [hoyte-2023]. The
group orders here are below $2^{128}$, where the leading-order estimate for
Wagner's algorithm is smaller by a factor of $2^{9.5}$ to $2^{11.5}$, so a
discrete-logarithm break makes the attack practical. The curve is therefore
derived from the beacon value as well, deterministically, by a public
candidate-selection procedure. A participant must be able to derive and check
the curve for each new beacon value.

## Design

The checksum uses hashing to the curve, point addition, and negation; scalar
multiplication occurs only in curve selection. The items are public, so the
arithmetic may run in variable time and items may be hashed by
try-and-increment; the map from an item's SHA-256 digest to the group need not
be uniform, only such that relations are hard to find within the subset it
reaches ([Adversary](docs/problem.md#adversary)). Over binary fields, the
complete formulas of Pornin [pornin-2022, pornin-2023] have an operation count
for addition independent of the extension degree $m$, negation is one field
addition, and a dense curve constant enters Pornin's extended coordinates only
through two multiplications by $\beta = B^{1/4}$ per addition. Scaling the
addends by the inverse of $b = B^{1/2}$ when they are decoded removes these,
leaving seven field multiplications and two squarings ($7M + 2S$) per mixed
addition.

The curve of a namespace is the first candidate derived from the beacon value
whose group order is an admissible cofactor times a prime $r$. A participant can
check that choice without counting points: a certificate gives $r$ and a
nonidentity point $Q$ with $rQ = O$ for the accepted candidate, and a witness of
an inadmissible order for each earlier candidate. Its conclusion is conditional
on the primality of $r$, which is tested to 24 fixed Miller–Rabin bases, not
proven ([Parameter selection and
verification](docs/problem.md#parameter-selection-and-verification)).

Random short Weierstrass and Edwards curves over prime fields are the
conservative alternatives: their field multiplication uses the integer
multiplier of every target, and a prime field has no subfield over which a Weil
descent applies.

With $m = 61$, $\mathrm{GF}(2^{122}) = \mathrm{GF}(2^{61})^2$ admits
Galbraith–Lin–Scott (GLS) curves, whose constants range over about $2^{61}$
values and whose automorphism lowers the rho cost by a factor $\sqrt{2}$ beyond
that of the negation map. A $\mathrm{GF}(2^{61})$ coordinate fits one 64-bit
word, and points have 16-byte encodings; the 107-bit and 109-bit families have
14-byte encodings. Weil descent constrains the choice of binary field.
[Candidate constructions](docs/constructions.md) lists the constructions
compared and their operation counts.

In a RIBLT, an item's point is added to about $2 \ln m$ of the first $m$ coded
symbols, and the item is hashed once per namespace, so additions outnumber
hashes ([Workload](docs/workload.md)). The simulator reproduces the forgery of
yangl1996/riblt#3 against the XOR checksum, and every curve group rejects the
forged cells.

## Findings

The figures below are from runs of source revision ef9e82d3 in the full profile:
on an Apple M4, run `01a1196a-88e4-7cdd-97ef-f46747dd8ab0`, unless they are
given for an Intel Xeon W-1390P, run `01a11cd2-6ca3-73a7-8501-245c98af23fc`,
built for x86-64-v3 with PCLMULQDQ. Those of the portable build are from run
`01a11a5b-9046-7637-b5af-cfa347a1a47d` on the M4, built for
`-C target-cpu=generic`, and run `01a12137-de75-762e-8d58-604336dcc781` on the
Xeon, built for `-C target-cpu=x86-64`. Batched hashing processes 1024 items per
batch, addition costs are throughput figures, and RIBLT encoding is per item,
for $n = 3500$ items into a sketch of $m = 1350$ cells.

Gains from the Go-to-Rust port or the mapping-generator change (xoshiro256++ in
place of ChaCha8) are not attributed to this project; matched comparisons in the
Rust harness measure checksum substitution ([comparison
conditions](docs/workload.md#comparison-with-the-reference-implementations)).

- The binary families over $\mathrm{GF}(2^{109})$ and $\mathrm{GF}(2^{127})$
  lead. Their fastest representation, unscaled accumulators with Pornin's
  complete formulas, encodes into a RIBLT sketch at 282 and 288 ns per item and
  adds in 9.8 and 10.7 ns, against 165 ns per item for the XOR baseline and 5.84
  µs for ristretto255. The incomplete λ-projective variants add in 12.0 and 14.9
  ns and encode at 303 and 326 ns; the extended coordinates, also complete, at
  324 and 329 ns.
- $\mathrm{GF}(2^{109})$ encodes a point in 14 bytes, against 16 for
  $\mathrm{GF}(2^{127})$, so each coded symbol's checksum is two bytes shorter,
  at a nominal rho cost of $2^{53.8}$ against $2^{62.8}$. Both fields occupy two
  64-bit limbs, and the smaller one is only slightly faster on the M4 and
  slightly slower on the Xeon; [Code generation](docs/codegen.md) (section 4)
  compares their reductions, which differ by modulus.
- Pornin's map is the faster hash, batched or not. To the unscaled addend, in
  batches, it takes 60.7 and 57.2 ns per item over $\mathrm{GF}(2^{109})$ and
  $\mathrm{GF}(2^{127})$, against 81.2 and 85.8 ns for try-and-increment; one
  item at a time, 533 and 461 ns against 874 and 869 ns, mainly for the field
  inversion. The RIBLT benchmarks therefore hash every binary family by Pornin's
  map, as the plan derived from the group operations selects ([The RIBLT
  plan](docs/methodology.md#the-riblt-plan)). Summing two maps costs 916 ns over
  $\mathrm{GF}(2^{127})$.
- Over $\mathrm{GF}(2^{122})$ the GLS-shaped and dense families encode at 306
  and 308 ns; with extended coordinates the GLS family adds about 6% faster
  (16.8 against 17.8 ns), at a nominal rho cost 0.5 bits lower ($2^{59.8}$
  against $2^{60.3}$).
- The prime-field families hash an item in 0.62 to 1.92 µs, mainly for the
  square root, which batching does not share, and add in 12.5 to 45.4 ns across
  the implemented formulas. The fastest of them in RIBLT encoding, twisted
  Edwards over $\mathrm{GF}(p^2)$ with $p = 2^{61} - 1$, takes 900 ns per item.
- XOR of SHA-256 digests, the insecure baseline, costs 24.2 ns per item and 0.38
  ns per cell. ristretto255 hashes an item in 4.72 µs and secp256k1 in 4.49 µs,
  185 to 195 times the baseline.
- The binary-field figures require carry-less multiplication (PMULL or
  PCLMULQDQ) to be enabled in the build; the portable backend executes 13 to 19
  times the instructions ([Code generation](docs/codegen.md)). Built without it,
  $\mathrm{GF}(2^{127})$ stays in the class of the prime-field curves: its
  families encode into the sketch at 1.90 to 2.43 µs per item, against 1.42 to
  1.68 µs for those over $\mathbb{F}_{2^{127} - 1}$. The families over
  $\mathrm{GF}(2^{122})$ take 2.90 to 4.12 µs and those over
  $\mathrm{GF}(2^{109})$ 4.53 to 6.06 µs, up to about four times the prime-field
  curves. On the Xeon $\mathrm{GF}(2^{127})$ takes 3.60 to 4.73 µs, about twice
  the 1.88 to 2.24 µs over $\mathbb{F}_{2^{127} - 1}$, $\mathrm{GF}(2^{122})$
  4.45 to 6.11 µs and $\mathrm{GF}(2^{109})$ 7.69 to 10.15 µs.
- Hashing an item anew in each namespace need not repeat SHA-256
  ([Adversary](docs/problem.md#adversary)): projecting its 32-byte identifier,
  the item's SHA-256 at 33.1 ns, to 128 bits under the namespace's key costs 5.1
  ns over $\mathrm{GF}(2^{127})$, 9.6 ns over $\mathbb{F}_{2^{127} - 1}$ and
  16.5 ns over $\mathbb{F}_{2^{130} - 5}$, against 33.0 ns for the salted
  SHA-256 of the identifier. Keyed by projections of identifiers computed
  beforehand, the unscaled accumulators over $\mathrm{GF}(2^{109})$ and
  $\mathrm{GF}(2^{127})$ encode at 222 and 227 ns per item, against 282 and 288
  ns by salted SHA-256 of the items ([Costs by
  lifetime](docs/methodology.md#costs-by-lifetime)). Without carry-less
  multiplication the projection over $\mathrm{GF}(2^{127})$ costs 60.7 ns on the
  M4 and 114 ns on the Xeon, more than the salted hash, and the portable build's
  RIBLT plan projects over $\mathbb{F}_{2^{127} - 1}$ instead.
- Over four seeds, finding a $\mathrm{GF}(2^{127})$ curve took 2.4 to 23.6 ms
  and certifying it, search included, 2.4 to 23.8 ms. Verifying a certificate
  counts no points and takes 55 to 68 µs across the families, 45 to 55% of it in
  the embedding-degree bound, or up to 0.59 ms when it also checks every
  rejected candidate. For the prime-field families the implemented prover counts
  points with PARI.
- On the Xeon the unscaled accumulators over $\mathrm{GF}(2^{109})$ and
  $\mathrm{GF}(2^{127})$ encode at 422 and 414 ns per item and add in 13.5 and
  13.2 ns, against 218 ns per item for the XOR baseline and 7.69 µs for
  ristretto255. Pornin's map remains the faster hash: in batches it takes 110
  and 98.0 ns per item against 125 ns for try-and-increment, and one item at a
  time 657 and 249 ns against 1154 and 401 ns, over $\mathrm{GF}(2^{127})$ about
  half the M4's cost. The prime-field families hash an item in 0.70 to 2.21 µs,
  and the fastest of them encodes at 1.18 µs. SHA-256 takes 53.4 ns, so
  ristretto255 and secp256k1, at 6.03 and 6.17 µs, hash at 113 to 116 times the
  baseline. Projecting an identifier costs 9.7 to 26.0 ns by field, against 51.4
  ns for its salted SHA-256, and keyed by projections the unscaled accumulators
  encode at 331 and 320 ns. Verifying a certificate takes 70 to 91 µs.
- The nominal rho costs of the families range from $2^{52.3}$ to $2^{63.3}$. The
  documented attacks obtain discrete logarithms and then solve additive
  relations; their costs and limitations are given in [Known
  weaknesses](docs/ecc_security.md#known-weaknesses).

## License

[MIT](./LICENSE)
