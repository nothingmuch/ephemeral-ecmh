use super::criteria::OrderError;
/// How to witness a counted order's rejection, or certify its acceptance.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Verdict {
    /// #E = cofactor * r with r passing the probable-prime policy.
    Accept(u128),
    /// The order fails criteria for which this format has no rejection witness.
    Inadmissible(OrderError),
    /// Rejected by an odd-order or family-specific cofactor witness.
    Reject(u128),
    /// Rejected by a point of order the smallest prime factor of this odd
    /// composite (which has none below 2^10).
    Composite(u128),
}
