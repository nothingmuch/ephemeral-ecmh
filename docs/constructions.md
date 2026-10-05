---
type: Chapter
title: Candidate constructions
description: The fields, curve models, accumulators, encodings and hash maps compared, and the models excluded.
tags: [elliptic-curves, binary-fields, edwards, weierstrass, hash-to-curve]
sources:
  - id: bibliography
    resource: ../references.bib
    title: Bibliography, annotated in literature.md
generated: { by: claude-code/claude-opus-5-5 }
contributors: [codex/gpt-6, claude-code/claude-fable-5-1]
status: draft
---

# Candidate constructions

| Construction | Field | Representation choices | Sources |
|---|---|---|---|
| Ordinary binary curves | $\mathrm{GF}(2^{127})$, $\mathrm{GF}(2^{109})$ | Extended, unscaled or lambda-projective accumulation; subgroup encodings, by $x$ or by Pornin's $w$ | [pornin-2022, oliveira-et-al-2013, crrl, pornin-2024] |
| Binary curves over a quadratic tower | $\mathrm{GF}((2^{61})^2)$ | Dense or subfield (GLS) curve constants; the same accumulators and encodings | [galbraith-lin-scott-2009, pornin-2023, oliveira-et-al-2013] |
| Complete Edwards curves, $a = 1$ | $\mathbb{F}_p$, $p = 2^{127} - 1$ or $2^{107} - 1$; $\mathrm{GF}(p^2)$, $p = 2^{61} - 1$ | Montgomery affine points, Edwards accumulators, cached addends | [bernstein-lange-2007, hisil-et-al-2008, montgomery-1987] |
| Short Weierstrass curves, $a = -3$ | The same three fields | Complete projective or Jacobian accumulation | [renes-costello-batina-2016, efd] |
| Complete twisted Edwards curves, $a = -1$ | $\mathbb{F}_p$, $p = 2^{128} - 275$; $\mathrm{GF}(p^2)$, $p = 2^{61} - 1$, $2^{64} - 59$ or $2^{64} - 2^{32} + 1$ | Quotient by the rational 2-torsion subgroup, a construction related to Decaf; canonical 16-byte encoding | [bernstein-et-al-2008, hisil-et-al-2008, hamburg-2015, pornin-2024] |
| ristretto255 (reference) | $\mathbb{F}_p$, $p = 2^{255} - 19$ | Upstream implementation | [rfc9496, bernstein-2006, curve25519-dalek] |
| secp256k1 (reference) | $\mathbb{F}_p$, $p = 2^{256} - 2^{32} - 977$ | Upstream implementation | [certicom-2010, libsecp256k1] |

The 107- and 109-bit families give 14-byte encodings. Their nominal rho costs
range from $2^{52.3}$ to $2^{53.8}$ group operations ([Security
considerations](ecc_security.md)). They show the cost of trading security margin
for width. The two 256-bit groups are fixed-curve references at conventional
security levels.

The variable-curve families hash by try-and-increment [boneh-lynn-shacham-2001].
Alternatives are Elligator 2 [bernstein-et-al-2013] on every Edwards and twisted
Edwards family, simplified SWU [rfc9380] on every Weierstrass family, and
Pornin's map [pornin-2023, crrl] on every binary family, each wherever its
preconditions hold; [`curve::h2c`](../src/curve/h2c/mod.rs) tabulates the cells
and the reason for each empty one. These are map-level alternatives. A single
map is not indifferentiable from a random oracle, which the checksum does not
require; the [Adversary](problem.md#adversary) section states what it does
require, and a sufficient condition that is not established for these maps or
for try-and-increment. The fixed-curve references use the maps named by their
benchmarks.
