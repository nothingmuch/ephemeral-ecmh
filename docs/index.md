---
okf_version: "0.2"
---

# Study

* [Background](background.md) - The RIBLT checksum, its linearity over $\mathbb{F}_2$ (yangl1996/riblt#3), why AdHash does not repair it, and keyed and group-based alternatives.
* [Problem statement](problem.md) - The adversary, the beacon-salted setting and its reference work factor, per-namespace curves against precomputation, curve selection and certificates, design criteria and scope.
* [Workload](workload.md) - The multiset digest, the operations a RIBLT performs on it, and their counts in repeated reconciliation.
* [Candidate constructions](constructions.md) - The fields, curve models, accumulators, encodings and hash maps compared, and the models excluded.
* [Security considerations](ecc_security.md) - Nominal rho costs per family, extension-field attacks, the limits of the security evidence, and the known weaknesses.
* [Measurement methodology](methodology.md) - What the benchmarks measure, the decision table and dominance rule, how to run checks, benchmarks and reports, and how results are published.

# Reference

* [Annotated literature](literature.md) - For each key of references.bib, what the work contributes to the evaluation of ECMH-based RIBLT checksums over per-namespace curves; in the book, headed by the formatted entry.
