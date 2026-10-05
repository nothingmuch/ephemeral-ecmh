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

**Discrete logarithms reduce the checksum to an insecure AdHash.** Write $G$ for
a generator of the prime-order component and $h_x$ for the discrete logarithm of
the hash of item $x$, so that its prime-order part is $h_x G$. An adversary who
knows $h_x$ for a collection of items converts a relation $\sum a_x H(x) = O$
into the congruence $\sum a_x h_x \equiv 0 \pmod{r}$, with a further condition
modulo the cofactor, which is at most 4. The checksum is then AdHash
[bellare-micciancio-1997] with modulus $r$ [maitin-shepard-et-al-2016, Section
8], the group order and not a modulus chosen for AdHash: for the families
compared, $r$ lies between about $2^{105}$ and $2^{127}$, where AdHash's
security requires thousands of bits (below). The following are leading-order
estimates of the resulting attacks:

- A nonzero relation with coefficients in $\lbrace -1, 0, 1 \rbrace$ gives two
  distinct sets with equal checksums. Wagner's algorithm [wagner-2002] finds one
  in about $2^{2\sqrt{\log_2 r}}$ operations on about as many hashed items,
  $2^{20.5}$ to $2^{22.5}$ here, for a relation of about $2^{\sqrt{\log_2 r}}$
  terms; list sizes and sorting contribute further factors. A relation with
  coefficients in $\lbrace 0, 1 \rbrace$ and right side the logarithm of a
  target checksum $T$ gives one set with checksum $T$; with coefficients in
  $\lbrace -1, 0, 1 \rbrace$ it gives two sets whose checksums differ by $T$.
- Relations with small integer coefficients give the corresponding multiset
  relations. The estimate of [maitin-shepard-et-al-2016, Appendix B] gives a
  coefficient vector of Euclidean norm about $2^{2\sqrt{\log_2 r \log_2 c}}$
  from $\sqrt{\log_2 r / \log_2 c}$ items, for a lattice reduction algorithm of
  Hermite factor constant $c$. With LLL ($c$ about 1.021), this is a
  heuristically expected norm of 12 to 15 from 59 to 65 items, and so no
  coefficient larger in magnitude. Such relations are forgeries only where the
  application admits multiplicities; they are not in general set relations.

For comparison, additive hashes with 128-bit security against Wagner's algorithm
need 2688 bits of output when a relation may have at most $2^{25}$ terms, and
4160 bits without that bound [mihajloska-gligoroski-samardjiska-2015, Section 4;
meyer-2022, Section V-B]. Against lattice reduction on multisets, at the Hermite
factor constant $c$ of about 1.007 that 128-bit security requires, AdHash needs
a modulus of over 400,000 bits [maitin-shepard-et-al-2016, Appendix B]. For an
additive fingerprint modulo $2^{256}$, where the leading-order estimate above is
about $2^{32}$ operations, Hoyte reports a $k$-dimensional birthday solver that
finds collisions in about 28 hours on eight cores, with 60 GB of memory and 1.5
TB of storage [hoyte-2023]. The checksum does not expose the residues modulo $r$
that these algorithms need until the logarithms are known. Once they are, the
modulus is $r$, at most $2^{127}$, and the estimate falls by a factor of
$2^{9.5}$ to $2^{11.5}$, low enough that a relation among the residues can be
found in practice.

In these attacks the cost is dominated by computing the logarithms; this bounds
the cost of a forgery from above, not from below. For $L$ logarithms in one
group, with $L$ well below $r^{1/4}$, rho with retained distinguished points
costs on the order of $\sqrt{L}$ times one logarithm [bernstein-lange-2012,
Section 1, after Kuhn and Struik]: about $2^3$ times the tabulated cost for the
multiset relation above, and $2^{10}$ to $2^{11}$ times for Wagner's algorithm,
as order-of-magnitude estimates. Whether lattice methods find set relations from
fewer logarithms is not assessed. The attack of
[yangl1996/riblt#3](https://github.com/yangl1996/riblt/issues/3) on the XOR
checksum is the case of this structure in which the relation is found by linear
algebra over $\mathbb{F}_2$ at negligible cost. The step from a relation to a
decoding failure, which depends on the cells the items map to, is not assessed
here ([Linearity of the checksum](background.md#linearity-of-the-checksum)).

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

**Hash to curve.** Items are hashed with tagged, salted SHA-256 and a counter
([`hash`](../src/hash.rs)). A deterministic map takes one 128-bit half of the
first digest (`hash_to_curve_map1`) or sums the images of both halves
(`hash_to_curve_map2`); try-and-increment takes the two halves of each digest in
order and increments the counter when both fail. The reduction of a candidate to
the field and the sign convention are defined per map. No RFC 9380 hash-to-curve
suite is implemented. Uniformity on the group is not required
([Adversary](problem.md#adversary)), but the weak-encoding hypotheses of the
reduction of [maitin-shepard-et-al-2016] are not established for these maps.

**Mapping-seed collisions in the RIBLT.** This weakness belongs to the RIBLT's
mapping, not to the checksum, and lies outside the research question; it is
stated because the RIBLT comparison rows use both mappings. The cells an item
maps to are determined by its digest under the mapping tag
([`riblt`](../src/riblt.rs)), through the gap law of the reference
implementation [yang-gilad-alizadeh-2024]. Two items whose generators start in
the same state map to the same cells in the same order for every prefix of the
coded symbols. If both lie in the symmetric difference, every cell they map to
holds both: after the other items are peeled, its count is even and its checksum
is the sum or difference of their two points, so no cell containing them becomes
pure and neither item is recovered, unless a checksum collision or false
singleton occurs independently. A finite peel returns failure; a rateless
decoder that requests coded symbols until decoding succeeds does not terminate.
No checksum collision is required, and the group's work factor does not enter.
The reference implementation seeds a 64-bit multiplicative generator with the
item's 64-bit hash. If the mapping input is public, two items with equal seeds
are found by a birthday search over about $2^{32}$ digests, below the rho cost
of every curve family compared, and an item sharing a chosen item's seed costs
about $2^{64}$. Between honest items, a coincidence has probability about
$n^2 / 2^{65}$ for $n$ items in the difference. Buchanan's implementation
[riblt-ecmh] replaces the generator with ChaCha8, citing this case
([yangl1996/riblt#3](https://github.com/yangl1996/riblt/issues/3), comment of
2026-02-24). The implemented default seeds xoshiro256++ with the whole 256-bit
map digest. Its state is the digest itself, so equal generator states require
equal digests, at a cost of about $2^{128}$ for a pair, with one exception in
$2^{256}$: the all-zero digest is reseeded from SplitMix64 and shares a state
with the one digest that encodes that state. Since the generator runs through
its $2^{256} - 1$ nonzero states in one cycle, distinct states yield distinct
sample streams; whether two streams place their items in the same cells of a
given table is the gap recurrence's matter, as for any generator. The generator
need not be cryptographic: the digest is already a salted hash, and the mapping
is public. The gap law is unchanged, and ChaCha8 (`ChaCha8`) and the 64-bit
generator (`Mcg64`) remain as comparisons. Nor does the seed need 256 bits: a
128-bit seed puts equal seeds at about $2^{64}$ identifiers after the salt.
`Riblt::projected` seeds the generator with the whole projection of the item's
identifier under keys of the mapping's own ([Adversary](problem.md#adversary)),
in a field chosen for its cost apart from the checksum's. Its halves are a
linear map of the identifier's halves under two keys, invertible unless a key
is zero or the two are equal, so equal seeds require equal identifiers; under
$\mathbb{F}_{2^{130} - 5}$, whose halves are truncated to 128 bits, they
require a coincidence of 256 bits. Either costs about $2^{128}$. The
projection's cost is benched, but it is not the default. The mapping need not
depend only on the public salt
of the checksum: a schedule keyed by a secret shared between the peers prevents
a third party from searching offline for equal seeds, even with a 64-bit
generator, while the checksum's namespace remains public. The coincidence
between honest items is unchanged by the key, since it does not depend on who
can compute seeds. Schedules that agree only on a prefix delay a rateless
decoder without stalling it; that delay, and any other a third party induces
through the schedule, concerns the sketch rather than its checksum and is not
assessed ([Questions outside this
study](problem.md#questions-outside-this-study)).

**Side channels.** Hashing, field inversion, and batched inversion run in
variable time. This is acceptable only because the hashed items are public in
the setting considered ([Adversary](problem.md#adversary)); the code is
unsuitable for secret inputs. Power and electromagnetic leakage are not
considered.

**Assumptions on the application.** The beacon value is assumed unpredictable
until the epoch begins; withholding, delay, and grinding by its producer, whose
cost Bonneau, Clark and Goldfeder bound for block hashes
[bonneau-clark-goldfeder-2015], are outside this study ([Questions outside this
study](problem.md#questions-outside-this-study)). A checksum collision
causes a decoding error within one namespace. The protocol is assumed to
validate the items it obtains independently, so that such an error, or a curve
whose discrete logarithms are computed, is a failure of liveness until the next
beacon value selects a new curve ([Setting and reference work
factor](problem.md#setting-and-reference-work-factor)); recovery is left to the
protocol.
