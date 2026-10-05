//! Tagged, salted SHA-256 with a per-salt midstate.
//!
//! `SHA256(SHA256(tag) || SHA256(tag) || salt || 0^32 || msg || ctr_le32)`.
//! The first 128 bytes are exactly two compression blocks, so each salt
//! costs two compressions once, and a message of up to 51 bytes costs one.

use sha2::{Digest, Sha256};

#[derive(Clone)]
pub struct Salted(Sha256);

impl Salted {
    pub fn new(tag: &[u8], salt: &[u8; 32]) -> Self {
        let t = Sha256::digest(tag);
        let mut h = Sha256::new();
        h.update(t);
        h.update(t);
        h.update(salt);
        h.update([0u8; 32]);
        Self(h)
    }

    pub fn digest(&self, msg: &[u8], ctr: u32) -> [u8; 32] {
        let mut h = self.0.clone();
        h.update(msg);
        h.update(ctr.to_le_bytes());
        h.finalize().into()
    }
}

/// The two 128-bit little-endian halves of a digest.
pub fn halves(d: &[u8; 32]) -> [u128; 2] {
    let (lo, hi) = d.split_at(16);
    [
        u128::from_le_bytes(lo.try_into().unwrap()),
        u128::from_le_bytes(hi.try_into().unwrap()),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn midstate_matches_one_shot(
            tag in prop::collection::vec(any::<u8>(), 0..80),
            salt in any::<[u8; 32]>(),
            msg in prop::collection::vec(any::<u8>(), 0..200),
            ctr in any::<u32>(),
        ) {
            let t = Sha256::digest(&tag);
            let one_shot: [u8; 32] = Sha256::new()
                .chain_update(t)
                .chain_update(t)
                .chain_update(salt)
                .chain_update([0u8; 32])
                .chain_update(&msg)
                .chain_update(ctr.to_le_bytes())
                .finalize()
                .into();
            prop_assert_eq!(Salted::new(&tag, &salt).digest(&msg, ctr), one_shot);
        }
    }
}
