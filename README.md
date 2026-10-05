# Ephemeral ECMH

## Authorship and disclaimer

The author formulated the research question, defined the threat model
and security arguments not related to
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

## License

[MIT](./LICENSE)
