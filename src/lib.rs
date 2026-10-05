//! Elliptic-curve multiset hashing and RIBLT operation costs.
//!
//! The study targets reference work factors between 2^48 and 2^64 group
//! operations and compact curve-point identifiers. A codec must cover the
//! additive closure of the hash's output range, including the identity and
//! signed sums. A namespace is the set of checksums computed under one
//! beacon value, which determines both the hash salt and the curve; an
//! epoch is an application's instance of a namespace. Curve selection is
//! amortized across a namespace; the principal costs are variable-time
//! hashing to reusable addends, addition, subtraction, and encoding on
//! public data.
//!
