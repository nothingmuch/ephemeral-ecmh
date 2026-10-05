//! The encoding `edwards` and `weier` share, and hashing by it: a point is
//! x and sgn0(y), x | sgn0(y) << BITS in `Signed::Bytes` with x packed
//! (`Packed`), and the all-ones x, which no packed x is, stands for O. It
//! needs BITS + 1 bits of `Bytes`; the rest must be zero.

use crate::field::Packed;
use core::fmt::Debug;

/// A field whose points the signed-x codec encodes, and in how many bytes.
/// The codec carries x and the sign in one u128, so BITS < 128 and `Bytes`
/// is at most 16 bytes.
pub trait Signed: Packed {
    /// At least BITS + 1 bits and at most 16 bytes, little-endian.
    type Bytes: Copy + Eq + Debug + Default + AsRef<[u8]> + AsMut<[u8]>;

    /// v's low bytes, as many as fit.
    fn to_bytes(v: u128) -> Self::Bytes {
        let mut b = Self::Bytes::default();
        let n = b.as_ref().len();
        b.as_mut().copy_from_slice(&v.to_le_bytes()[..n]);
        b
    }

    /// Every bit of c, so callers can reject what they don't use.
    fn from_bytes(c: Self::Bytes) -> u128 {
        let mut v = [0; 16];
        v[..c.as_ref().len()].copy_from_slice(c.as_ref());
        u128::from_le_bytes(v)
    }
}

impl Signed for crate::field::fp127::Fp {
    type Bytes = [u8; 16];
}

/// The x bits, and O's x.
fn x_mask<F: Packed>() -> u128 {
    u128::MAX >> (128 - F::BITS)
}

pub(crate) fn mask<F: Signed>() -> u128 {
    u128::MAX >> (127 - F::BITS)
}

pub(crate) fn sign<F: Signed>(c: u128) -> bool {
    c >> F::BITS & 1 == 1
}

/// The low BITS bits of c, reduced.
pub(crate) fn x<F: Packed>(c: u128) -> F {
    F::reduce(c & x_mask::<F>())
}

/// `None` is O.
pub(crate) fn encode<F: Signed>(p: Option<(F, F)>) -> F::Bytes {
    F::to_bytes(match p {
        None => x_mask::<F>(),
        Some((x, y)) => x.pack() | (y.sgn0() as u128) << F::BITS,
    })
}

/// The point c encodes on y^2 = rhs(x), if any, by one square root; in it,
/// `None` is O.
pub(crate) fn decode<F: Signed>(c: u128, rhs: impl FnOnce(F) -> F) -> Option<Option<(F, F)>> {
    if c > mask::<F>() {
        return None;
    }
    let (v, sign) = (c & x_mask::<F>(), sign::<F>(c));
    if v == x_mask::<F>() {
        return (!sign).then_some(None);
    }
    let x = F::unpack(v)?;
    let mut y = rhs(x).sqrt()?;
    if y.is_zero() && sign {
        return None;
    }
    if y.sgn0() != sign {
        y = -y;
    }
    Some(Some((x, y)))
}

/// The x and sign bits of a digest half, if x is packed. Hashing skips
/// the rest: O's x would decode to O, and the others aren't canonical.
pub(crate) fn candidate<F: Signed>(c: u128) -> Option<u128> {
    let c = c & mask::<F>();
    F::unpack(c & x_mask::<F>()).map(|_| c)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::fp127;
    use proptest::prelude::*;

    /// The x, the sign and O's x fit in `Bytes`, which fits in a u128, keeps
    /// exactly the bits that fit and reads them back losslessly.
    fn layout<F: Signed>(v: u128) -> Result<(), TestCaseError> {
        let n = 8 * core::mem::size_of::<F::Bytes>() as u32;
        prop_assert!(F::BITS < n && n <= 128);
        prop_assert!(F::unpack(x_mask::<F>()).is_none());
        let kept = v & (u128::MAX >> (128 - n));
        prop_assert_eq!(F::from_bytes(F::to_bytes(v)), kept);
        Ok(())
    }

    proptest! {
        #[test]
        fn layout_fp127(v in any::<u128>()) {
            layout::<fp127::Fp>(v)?;
        }
    }
}
