---
type: Chapter
title: Background
description: The RIBLT checksum, its linearity over F_2 (yangl1996/riblt#3), why AdHash does not repair it, and keyed and group-based alternatives.
tags: [riblt, checksum, xhash, adhash, ecmh]
sources:
  - id: riblt-issue-3
    resource: https://github.com/yangl1996/riblt/issues/3
    title: "yangl1996/riblt#3: Adversarial threat model"
  - id: bibliography
    resource: ../references.bib
    title: Bibliography, annotated in literature.md
generated: { by: claude-code/claude-opus-5-5 }
contributors: [codex/gpt-6, claude-code/claude-fable-5-1]
status: draft
---

# Background

## The RIBLT checksum

A RIBLT reconciles two sets: the decoder recovers their symmetric difference
rather than reading back stored values. A coded symbol holds three sums over the
items mapped to it: the XOR of the items, an item count, and a checksum. The
checksum is the XOR of fixed-width item hashes. Decoding peels a *pure* symbol,
which holds exactly one item, and uses the checksum to recognize purity.

In the construction studied here the checksum is an ECMH: a coded symbol holds
the signed sum of the items' salted hashes to an elliptic curve, in place of the
XOR of fixed-width hashes. The item field remains the XOR of the item bytes. A
variant in which each item is itself a curve point,
a 128-bit identifier recovered from the point sum, is not implemented.

## Linearity of the checksum

[yangl1996/riblt#3](https://github.com/yangl1996/riblt/issues/3) observes that
this checksum is linear over $\mathbb{F}_2$. An adversary who can evaluate the
item hash can find a set of distinct, non-colliding items whose hashes XOR to a
target value. Gaussian elimination over slightly more candidates than the hash
width suffices, at negligible cost. The weakness lies in the XOR combination
(XHASH [bellare-micciancio-1997]) rather than in the hash function, so it
persists with SHA-256. The issue demonstrates two constructed sets that the
decoder wrongly reports as equal. Whether such relations can also suppress one
chosen item across every coded symbol it maps to is left open there. In either
case the checksum provides no resistance to relation finding.

Replacing XOR with addition modulo $2^{64}$ (AdHash [bellare-micciancio-1997])
does not remove the weakness. Relations with small integer coefficients among
hash values modulo a 64-bit integer are found efficiently by lattice reduction
[maitin-shepard-et-al-2016, Appendix B]. Given an unrestricted number of terms,
Wagner's generalized birthday algorithm [wagner-2002] finds them in time
subexponential in the size of the modulus. The security of AdHash rests on the
hardness of a weighted knapsack problem, which requires a modulus of thousands
of bits [bellare-micciancio-1997, wagner-2002]; a variant of one incremental
hash of this kind was proposed with its modulus raised from $2^{160}$ to
$2^{1600}$ in response to Wagner's algorithm [phan-wagner-2006, Section 2]
([Known weaknesses](ecc_security.md#known-weaknesses)).

## Keyed and group-based checksums

Keying the item hash with a secret prevents the attack by parties that cannot
evaluate the keyed hash. Short keyed hashes such as SipHash
[aumasson-bernstein-2012] serve this purpose against hash flooding of hash
tables, where the key does not leave one process: the adversary neither holds
the key nor obtains keyed hashes of candidates of its choosing. Under such a key
a RIBLT gives up universality: coded symbols computed under one key cannot serve
peers that do not share it. This repository assumes less. The salt and the set
are public, the set accepts items from anyone, the adversary included, and the
salt is a value of a randomness beacon outside the adversary's control, renewed
often enough that discrete logarithms in the checksum group cannot be computed
while it is in use ([Setting and reference work
factor](problem.md#setting-and-reference-work-factor)).

The issue discussion proposes computing the checksum in a finite abelian group
in which relations with small coefficients among hash outputs are hard to find.
ECMH is such a construction. An implementation over ristretto255 [riblt-ecmh]
reports 256-bit checksums in place of 64-bit ones and a higher end-to-end cost.
It also replaces the pseudorandom generator that maps items to coded symbols, so
its timings do not isolate the cost of the group.
