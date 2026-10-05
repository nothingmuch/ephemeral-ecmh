//! Portable backend: carry-less products from integer multiplications,
//! folds by shifts. Always built in tests, so it is checked on every target.

use super::tables::RHO;

#[derive(Clone, Copy)]
pub struct Gf([u64; 2]);

impl Gf {
    /// a0 + a1 u, low word first; any 64-bit representatives.
    pub const fn w64le(x0: u64, x1: u64) -> Self {
        Self([x0, x1])
    }

    /// The two 64-bit representatives, not reduced.
    #[inline(always)]
    pub fn limbs(self) -> [u64; 2] {
        self.0
    }
}

#[inline(always)]
pub fn add(a: Gf, b: Gf) -> Gf {
    Gf([a.0[0] ^ b.0[0], a.0[1] ^ b.0[1]])
}

/// 32 x 32 -> 64 with data bits 4 apart: at most 8 of them meet in a
/// position, so integer carries stay inside the 3 spare bits.
fn clmul32(x: u32, y: u32) -> u64 {
    const M: [u64; 4] = [
        0x1111_1111_1111_1111,
        0x2222_2222_2222_2222,
        0x4444_4444_4444_4444,
        0x8888_8888_8888_8888,
    ];
    let (xs, ys) = (M.map(|m| x as u64 & m), M.map(|m| y as u64 & m));
    let mut z = 0;
    for (i, m) in M.iter().enumerate() {
        let t = (0..4).fold(0, |t, j| t ^ (xs[j] * ys[(i + 4 - j) % 4]));
        z |= t & m;
    }
    z
}

/// 64 x 64 -> 128, one Karatsuba level over `clmul32`.
pub fn clmul(a: u64, b: u64) -> u128 {
    let (a0, a1, b0, b1) = (a as u32, (a >> 32) as u32, b as u32, (b >> 32) as u32);
    let (lo, hi) = (clmul32(a0, b0), clmul32(a1, b1));
    let mid = clmul32(a0 ^ a1, b0 ^ b1) ^ lo ^ hi;
    lo as u128 ^ (mid as u128) << 32 ^ (hi as u128) << 64
}

/// a^2 as a polynomial: bit i moves to bit 2i.
pub fn spread(a: u64) -> u128 {
    let s = |x: u32| {
        let mut x = x as u64;
        x = (x | x << 16) & 0x0000_FFFF_0000_FFFF;
        x = (x | x << 8) & 0x00FF_00FF_00FF_00FF;
        x = (x | x << 4) & 0x0F0F_0F0F_0F0F_0F0F;
        x = (x | x << 2) & 0x3333_3333_3333_3333;
        (x | x << 1) & 0x5555_5555_5555_5555
    };
    s(a as u32) as u128 | (s((a >> 32) as u32) as u128) << 64
}

// Keep the hard-coded reduction shifts consistent with the generated polynomial.
const _: () = assert!(RHO == 1 << 26 | 1 << 18 | 1 << 8 | 1 << 3);

/// lo + hi z^64 mod z^3 f: T = floor(hi rho / z^64), then
/// lo + ((hi + T) rho mod z^64) (see `gf2_61`).
#[inline(always)]
fn reduce(d: u128) -> u64 {
    let (lo, hi) = (d as u64, (d >> 64) as u64);
    let t = hi >> 38 ^ hi >> 46 ^ hi >> 56 ^ hi >> 61;
    let g = hi ^ t;
    lo ^ g << 26 ^ g << 18 ^ g << 8 ^ g << 3
}

#[inline(always)]
pub fn mul61(a: u64, b: u64) -> u64 {
    reduce(clmul(a, b))
}

#[inline(always)]
pub fn square61(a: u64) -> u64 {
    reduce(spread(a))
}

/// Karatsuba over the three unreduced products, then one fold per word.
pub fn mul(a: Gf, b: Gf) -> Gf {
    let ([a0, a1], [b0, b1]) = (a.0, b.0);
    let (p00, p11) = (clmul(a0, b0), clmul(a1, b1));
    let pm = clmul(a0 ^ a1, b0 ^ b1);
    Gf([reduce(p00 ^ p11), reduce(pm ^ p00)])
}

pub fn square(a: Gf) -> Gf {
    let [a0, a1] = a.0;
    Gf([square61(a0 ^ a1), square61(a1)])
}
