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

Binary Edwards curves [bernstein-lange-farashahi-2008] are excluded on inspected
operation counts under equal field arithmetic; they are not implemented or
measured. Every ordinary binary curve is birationally equivalent to a complete
binary Edwards curve [bernstein-lange-farashahi-2008, Theorem 4.3], so the
exclusion is one of cost, not availability. Their complete mixed addition costs
$13M + 3S + 3D$, and $13M + 2S + 2D$ with the formulas of Kim, Lee and Negre
[kim-lee-negre-2014]; the complete unscaled formulas used here [pornin-2022]
cost $7M + 2S$. $M$, $S$ and $D$ denote a field multiplication, a squaring and a
multiplication by a curve constant, counted as a full multiplication for generic
dense parameters. Mixed addition is the most frequent operation of the workload
([Workload](workload.md)). The binary Edwards formulas with lower operation
counts [bernstein-2009] are differential: they compute $P + Q$ from $P$, $Q$ and
$P - Q$, as a Montgomery ladder supplies them, and a coded-symbol update adds an
independent hashed point with no known difference.

The quadratic extensions $\mathrm{GF}(p^2)$ fit 16-byte encodings with one-word
base arithmetic. Goldilocks [polygon-zero-2022] uses Plonky3's implementation
[plonky3]; its 2-adicity of 32 makes square roots, and so decoding and hashing,
slower than over $2^{64} - 59$; the benchmark report's Goldilocks rows measure
the difference. The 31-bit primes of current proof systems, M31
[haboeck-levit-papini-2024] and BabyBear [bruestle-gafni-2023], would need
extensions of degree 4 to reach the same group size. Degree-4 extensions admit
the heuristic index-calculus asymptotic $\widetilde{O}(q^{3/2})$ [gaudry-2009],
below the generic $q^2$, so group size alone does not justify the target
security there; these fields are excluded from the comparison.

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
benchmarks. The RIBLT benchmarks hash each family by the construction its run
measured cheapest ([The RIBLT plan](methodology.md#the-riblt-plan)).
