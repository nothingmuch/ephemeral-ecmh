//! Portable backend: carry-less products from integer multiplications.
//! Always built in tests, so it is checked on every target.

use super::tables::RHO;

#[derive(Clone, Copy)]
pub struct Gf([u64; 2]);

impl Gf {
    /// Low limb first; any 128 bits (reduced lazily).
    pub const fn w64le(x0: u64, x1: u64) -> Self {
        Self([x0, x1])
    }

    /// The 128-bit representative, not reduced.
    #[inline(always)]
    pub fn limbs(self) -> [u64; 2] {
        self.0
    }
}

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

/// d0 + d1 z^128 mod z^19 f: fold d1 onto rho twice (see the parent module).
fn reduce(d0: u128, d1: u128) -> Gf {
    let c = clmul((d1 >> 64) as u64, RHO);
    let e = clmul(d1 as u64 ^ (c >> 64) as u64, RHO);
    let v = d0 ^ e ^ c << 64;
    Gf([v as u64, (v >> 64) as u64])
}

pub fn mul(a: Gf, b: Gf) -> Gf {
    let ([a0, a1], [b0, b1]) = (a.0, b.0);
    let c1 = clmul(a0, b1) ^ clmul(a1, b0);
    reduce(clmul(a0, b0) ^ c1 << 64, clmul(a1, b1) ^ c1 >> 64)
}

pub fn square(a: Gf) -> Gf {
    reduce(clmul(a.0[0], a.0[0]), clmul(a.0[1], a.0[1]))
}
