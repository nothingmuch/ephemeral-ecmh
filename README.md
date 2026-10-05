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
hashing an item to them costs far more than the SHA-256 hash of the XOR
checksum. A shorter horizon
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

## License

[MIT](./LICENSE)
