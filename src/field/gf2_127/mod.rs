//! GF(2^127) = F_2\[z\]/(z^127 + z^63 + 1): the field of `binary127`.
//!
//! - Modulus: the trinomial z^127 + z^63 + 1, so GF(2^127) is a degree-127
//!   extension of F_2 directly, with no subfield tower in between.
//! - Representation: crrl's `GFb127`, re-exported as `Gf`: two 64-bit
//!   words, reduced lazily, so an element may sit in 128 bits; `to_u128`
//!   normalizes to the canonical 127-bit polynomial.
//! - Backend: crrl picks it at compile time. PMULL (`gfb254_arm64pmull`) on
//!   aarch64 with the `aes` target feature, PCLMULQDQ (`gfb254_x86clmul`)
//!   on x86_64 with `pclmulqdq` and `sse4.1`, else portable 64-bit integer
//!   code (`gfb254_m64`). The field arithmetic is crrl's; this module adds
//!   conversions, alternative half-trace tables, and tests.

pub use crrl::field::GFb127 as Gf;

mod tables;

pub const MASK127: u128 = (1 << 127) - 1;

/// Canonical 127-bit integer form (bit i = coefficient of z^i).
#[inline]
pub fn to_u128(x: Gf) -> u128 {
    u128::from_le_bytes(x.encode())
}

/// Bits above 126 are dropped (not reduced), so the map is onto and exact.
#[inline]
pub fn from_u128(v: u128) -> Gf {
    let v = v & MASK127;
    Gf::w64le(v as u64, (v >> 64) as u64)
}

#[inline]
pub fn is_zero(x: Gf) -> bool {
    x.iszero() != 0
}

#[inline]
pub fn eq(a: Gf, b: Gf) -> bool {
    a.equals(b) != 0
}

static HT8: [[u128; 256]; 16] = super::window_tables(&tables::HALFTRACE);
static HT4: [[u128; 16]; 32] = super::window_tables(&tables::HALFTRACE);

/// The half-trace by byte tables: 16 lookups into 64 KiB, trading cache
/// space for fewer lookups than crrl's per-bit tables.
/// `benches/compare.rs` measures this tradeoff against crrl and nibble tables.
#[inline(always)]
pub fn halftrace8(x: Gf) -> Gf {
    let v = to_u128(x);
    let mut h = 0;
    for (j, t) in HT8.iter().enumerate() {
        h ^= t[(v >> (8 * j)) as u8 as usize];
    }
    from_u128(h)
}

/// The half-trace by nibble tables: 32 lookups into 8 KiB, for when the
/// byte tables would crowd the cache.
#[inline(always)]
pub fn halftrace4(x: Gf) -> Gf {
    let v = to_u128(x);
    let mut h = 0;
    for (j, t) in HT4.iter().enumerate() {
        h ^= t[(v >> (4 * j)) as usize & 15];
    }
    from_u128(h)
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use proptest::prelude::*;

    /// Bit-serial reference multiply, independent of crrl's backends.
    pub fn mul_ref(a: u128, b: u128) -> u128 {
        let (mut lo, mut hi) = (0u128, 0u128);
        for i in 0..127 {
            if b >> i & 1 == 1 {
                lo ^= a << i;
                if i > 0 {
                    hi ^= a >> (128 - i);
                }
            }
        }
        // reduce bits 253..127 with z^127 = z^63 + 1
        for i in (127..254).rev() {
            let bit = if i >= 128 {
                hi >> (i - 128) & 1
            } else {
                lo >> 127 & 1
            };
            if bit == 1 {
                let j = i - 127;
                if i >= 128 {
                    hi ^= 1 << (i - 128);
                } else {
                    lo ^= 1 << 127;
                }
                for k in [j, j + 63] {
                    if k >= 128 {
                        hi ^= 1 << (k - 128);
                    } else {
                        lo ^= 1 << k;
                    }
                }
            }
        }
        lo
    }

    pub fn fe() -> impl Strategy<Value = u128> {
        any::<u128>().prop_map(|v| v & MASK127)
    }

    proptest! {
        #[test]
        fn roundtrip(a in fe()) {
            prop_assert_eq!(to_u128(from_u128(a)), a);
        }

        #[test]
        fn mul_matches_reference(a in fe(), b in fe()) {
            prop_assert_eq!(to_u128(from_u128(a) * from_u128(b)), mul_ref(a, b));
        }

        #[test]
        fn square_matches_mul(a in fe()) {
            let x = from_u128(a);
            prop_assert_eq!(to_u128(x.square()), mul_ref(a, a));
        }

        #[test]
        fn invert(a in fe()) {
            let x = from_u128(a);
            let y = x.invert();
            if a == 0 {
                prop_assert!(is_zero(y));
            } else {
                prop_assert_eq!(mul_ref(a, to_u128(y)), 1);
            }
        }

        #[test]
        fn sqrt(a in fe()) {
            let s = from_u128(a).sqrt();
            prop_assert_eq!(to_u128(s.square()), a);
        }

        #[test]
        fn trace_is_bit0(a in fe()) {
            // Tr(z^i) = 0 for 0 < i < 127, so Tr(v) = v_0.
            prop_assert_eq!(from_u128(a).trace() as u128, a & 1);
        }

        #[test]
        fn halftrace_solves_quadratic(a in fe()) {
            // H(v)^2 + H(v) = v + Tr(v)  (m odd)
            let v = from_u128(a);
            let h = v.halftrace();
            prop_assert_eq!(to_u128(h.square() + h), a ^ (a & 1));
        }

        #[test]
        fn halftrace_tables_match_crrl(a in any::<u128>()) {
            // any 128-bit representative, as crrl leaves them
            let v = Gf::w64le(a as u64, (a >> 64) as u64);
            let h = to_u128(v.halftrace());
            prop_assert_eq!(to_u128(halftrace8(v)), h);
            prop_assert_eq!(to_u128(halftrace4(v)), h);
        }

        #[test]
        fn distributive(a in fe(), b in fe(), c in fe()) {
            let (x, y, w) = (from_u128(a), from_u128(b), from_u128(c));
            prop_assert!(eq(x * (y + w), x * y + x * w));
        }
    }
}
