---
type: Chapter
title: Security considerations
description: Nominal rho costs per family, extension-field attacks, the limits of the security evidence, and the known weaknesses.
tags: [security, discrete-logarithm, rho, index-calculus, adhash]
sources:
  - id: bibliography
    resource: ../references.bib
    title: Bibliography, annotated in literature.md
generated: { by: claude-code/claude-opus-5-5 }
contributors: [codex/gpt-6, claude-code/claude-fable-5-1]
status: draft
---

# Security considerations

Encoded width, group order, and attack cost are distinct quantities. A 128-bit
encoding does not imply 128-bit collision resistance. Generic estimates must be
read together with extension-field structure, endomorphisms, the parameter
selection procedure, and the distribution of the hash map. A discrete-logarithm
estimate is not by itself a proof of hardness for every multiset relation
problem.

The table gives the nominal cost of rho with the negation map,
$\sqrt{\pi r / 4}$ group operations [pollard-1978, wiener-zuccherato-1998],
divided by a further $\sqrt{2}$ for the GLS automorphism group of order 4
[duursma-gaudry-morain-1999], for the largest prime $r$ dividing each family's
group order. The "Families compared" table of each run report computes the same
quantity from the $r$ recorded by that run's group suite:

| Family | Group order | $\log_2$ of rho cost |
|---|---|---|
| Binary, $\mathrm{GF}(2^{127})$ | $2r$, $r$ near $2^{126}$ | 62.8 |
| Binary, $\mathrm{GF}(2^{122})$, dense constant | $2r$, $r$ near $2^{121}$ | 60.3 |
| Binary, $\mathrm{GF}(2^{122})$, GLS | $2r$, $r$ near $2^{121}$ | 59.8 (automorphism group of order 4) |
| Binary, $\mathrm{GF}(2^{109})$ | $2r$, $r$ near $2^{108}$ | 53.8 |
| Weierstrass, $\mathbb{F}_p$, $p = 2^{127} - 1$ | $r$ | 63.3 |
| Edwards, $\mathbb{F}_p$, $p = 2^{127} - 1$ | $4r$ | 62.3 |
| Twisted Edwards, $\mathbb{F}_p$, $p = 2^{128} - 275$; $\mathrm{GF}(p^2)$, $p = 2^{64} - 59$ or Goldilocks | $4r$ | 62.8 |
| Weierstrass, $\mathrm{GF}(p^2)$, $p = 2^{61} - 1$ | $r$ | 60.8 |
| Edwards and twisted Edwards, $\mathrm{GF}(p^2)$, $p = 2^{61} - 1$ | $4r$ | 59.8 |
| Weierstrass, $\mathbb{F}_p$, $p = 2^{107} - 1$ | $r$ | 53.3 |
| Edwards, $\mathbb{F}_p$, $p = 2^{107} - 1$ | $4r$ | 52.3 |

For curves over $\mathrm{GF}(p^2)$, index calculus over the extension
[gaudry-2009] costs about $p$ group operations up to logarithmic factors, the
same as rho; the GHS descent for $\mathrm{GF}(2^{122}) = \mathrm{GF}(2^{61})^2$
[hankerson-karabina-menezes-2009] likewise gives no gain over rho. Candidates
over $\mathrm{GF}(p^2)$ whose $j$-invariant lies in $\mathbb{F}_p$ are excluded
as a conservative policy.

The setting is an ephemeral public namespace, with salt and curve derived from a
beacon value. It is distinct from a secret-keyed namespace, in which keyed XOR
checksums suffice under the adversary assumptions of [Keyed and group-based
checksums](background.md#keyed-and-group-based-checksums). Conclusions for one
setting do not transfer to the other. The implementation provides [tagged,
salted hashing](../src/hash.rs).

The [parameter-selection code](../src/curvegen/mod.rs) and the [reference
programs](../sage/) state the acceptance conditions for each family. Tests
combine field reference computations, independent affine group laws, Sage/PARI
fixtures, canonical-encoding checks, signed accumulation, and complete workload
checks. Sampled tests are evidence, not proof, and probable-prime tests do not
certify primality.

## Known weaknesses

Each weakness below is accepted at the reference work factor, stated as an
assumption on the application, or left open, as indicated.

**Generic work factor.** The 107- and 109-bit families have rho costs of
$2^{52.3}$ to $2^{53.8}$, in the lower part of the reference range.
Curve-specific precomputation is limited, not excluded, by deriving a curve per
beacon value ([Precomputation and per-namespace
curves](problem.md#precomputation-and-per-namespace-curves)).

**Weil descent.** A $\mathrm{GF}(2^{127})$ curve is reachable by the generalized
GHS descent with heuristic probability about $2^{-52}$ per beacon value, and no
check excludes it ([Parameter selection and
verification](problem.md#parameter-selection-and-verification)).

**Probable primality.** The verifier tests $r$ with Miller–Rabin to 24 fixed
bases. Composites passing all of them can be constructed, so no error bound
holds against a prover who chooses $r$. The conclusion that the group order is
the cofactor times $r$ is conditional on the primality of $r$.

**Pairing reductions.** The embedding degree of $r$ is checked to exceed
$2^{20}$, which excludes the reductions of Menezes–Okamoto–Vanstone and
Frey–Rück to that bound; the degree is not otherwise determined. The check, a
baby-step giant-step search for $q^k = 1 \bmod r$, takes about $2^{11}$ products
modulo $r$, which the curve-generation benchmarks time alone
(`curvegen/embedding`). Anomalous curves, of order $q$, are rejected.
