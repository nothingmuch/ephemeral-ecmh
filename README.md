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
hashing an item to them costs 150 to 200 times the SHA-256 hash of the XOR
checksum (see [Preliminary findings](#preliminary-findings)). A shorter horizon
and a failure confined to one namespace suggest that a work factor of $2^{48}$
to $2^{64}$ suffices, provided the application detects a sketch collision by
other means and recovers, for example by re-keying or by another reconciliation
method. This is an application assumption.

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

## Preliminary findings

The figures below are from run `01a10a32-10fd-76e5-bc0c-18195b7bccb0` on an
Apple M4, of source revision 9c592abc, in the quick profile, which takes fewer
samples than the full one; clean full-profile runs will replace them, and the
findings will be revised against those runs. Batched hashing processes 1024
items per batch, addition costs are throughput figures, and RIBLT encoding is
per item, for $n = 3500$ items into a sketch of $m = 1350$ cells.

Gains from the Go-to-Rust port or the mapping-generator change (xoshiro256++ in
place of ChaCha8) are not attributed to this project; matched comparisons in the
Rust harness measure checksum substitution ([comparison
conditions](docs/workload.md#comparison-with-the-reference-implementations)).

- The binary families over $\mathrm{GF}(2^{109})$ and $\mathrm{GF}(2^{127})$
  lead. Their fastest representation, unscaled accumulators with Pornin's
  complete formulas, encodes into a RIBLT sketch at 274 and 279 ns per item and
  adds in 9.8 and 10.6 ns, against 166 ns per item for the XOR baseline and 5.92
  µs for ristretto255. The incomplete λ-projective variants add in 11.9 and 15.0
  ns and encode at 292 and 319 ns; the extended coordinates, also complete, at
  315 and 326 ns.
- $\mathrm{GF}(2^{109})$ encodes a point in 14 bytes, against 16 for
  $\mathrm{GF}(2^{127})$, so each coded symbol's checksum is two bytes shorter,
  at a nominal rho cost of $2^{53.8}$ against $2^{62.8}$. Both fields occupy two
  64-bit limbs, and the smaller one is only slightly faster; [Code
  generation](docs/codegen.md) (section 4) compares their reductions, which
  differ by modulus.
- Pornin's map is the faster hash, batched or not. To the unscaled addend, in
  batches, it takes 59.1 and 57.9 ns per item over $\mathrm{GF}(2^{109})$ and
  $\mathrm{GF}(2^{127})$, against 75.5 and 79.8 ns for try-and-increment; one
  item at a time, 525 and 479 ns against 886 and 881 ns, mainly for the field
  inversion. The RIBLT benchmarks therefore hash every binary family by
  Pornin's map, as the plan derived from the group operations selects ([The
  RIBLT plan](docs/methodology.md#the-riblt-plan)). Summing two maps costs 957
  ns over $\mathrm{GF}(2^{127})$.
- Over $\mathrm{GF}(2^{122})$ the GLS-shaped and dense families encode at the
  same 302 ns; with extended coordinates the GLS family adds about 6% faster
  (16.8 against 17.8 ns), at a nominal rho cost 0.5 bits lower ($2^{59.8}$
  against $2^{60.3}$).
- The prime-field families hash an item in 0.66 to 2.01 µs, mainly for the
  square root, which batching does not share, and add in 12.4 to 45.4 ns across
  the implemented formulas. The fastest of them in RIBLT encoding, twisted
  Edwards over $\mathrm{GF}(p^2)$ with $p = 2^{61} - 1$, takes 917 ns per item.
- XOR of SHA-256 digests, the insecure baseline, costs 23.5 ns per item and 0.39
  ns per cell. ristretto255 hashes an item in 4.99 µs and secp256k1 in 4.70 µs,
  199 to 212 times the baseline.
- The binary-field figures require carry-less multiplication (PMULL or
  PCLMULQDQ) to be enabled in the build; the portable backend executes 13 to 19
  times the instructions ([Code generation](docs/codegen.md)).
- Over four seeds, finding a $\mathrm{GF}(2^{127})$ curve took 2.4 to 24.1 ms
  and certifying it, search included, 2.5 to 24.4 ms. Verifying a certificate
  counts no points and takes 56 to 68 µs across the families, or up to 0.62 ms
  when it also checks every rejected candidate. For the prime-field families the
  implemented prover counts points with PARI.
- The nominal rho costs of the families range from $2^{52.3}$ to $2^{63.3}$. The
  documented attacks obtain discrete logarithms and then solve additive
  relations; their costs and limitations are given in [Known
  weaknesses](docs/ecc_security.md#known-weaknesses).

## License

[MIT](./LICENSE)
