//! Z_q / 2^64 = (Z/2^64)\[t\]/(F): the degree-m unramified extension of the
//! 2-adic integers, one machine word per coefficient.
//!
//! F lifts the field's sparse modulus with 0/1 coefficients, so reduction is
//! a few subtractions per coefficient. The price is Frobenius: it is t -> t^2
//! only for the (dense) Teichmüller modulus, so here it is a dense m x m
//! matrix, built once per field.

use core::marker::PhantomData;
use core::ops::{Add, Mul, Neg, Sub};

/// Coefficient capacity; m <= 127 so that a residue mod 2 fits a u128.
pub const CAP: usize = 128;

/// A binary field F_2\[z\]/(z^M + sum of z^k over LOW). The polynomial
/// must be irreducible: Frobenius::new panics if it finds otherwise.
pub trait Modulus: Copy + Eq + core::fmt::Debug + 'static {
    const M: usize;
    /// The exponents below M, each once.
    const LOW: &'static [usize];
}

/// crate::field::gf2_127: GF(2^127) = F_2\[z\]/(z^127 + z^63 + 1).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Gf127;

impl Modulus for Gf127 {
    const M: usize = 127;
    const LOW: &'static [usize] = &[63, 0];
}

/// crate::field::gf2_109: GF(2^109) = F_2\[z\]/(z^109 + z^5 + z^4 + z^2 + 1).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Gf109;

impl Modulus for Gf109 {
    const M: usize = 109;
    const LOW: &'static [usize] = &[5, 4, 2, 0];
}

/// Multiplications, squarings and Frobenius matrix products on this
/// thread, so tests can pin the cost of a count.
#[cfg(test)]
pub(crate) mod ops {
    use std::cell::Cell;

    pub const MUL: usize = 0;
    pub const SQR: usize = 1;
    pub const FROB: usize = 2;

    thread_local!(static N: Cell<[u64; 3]> = const { Cell::new([0; 3]) });

    pub fn bump(i: usize) {
        N.with(|n| {
            let mut v = n.get();
            v[i] += 1;
            n.set(v)
        })
    }

    pub fn get() -> [u64; 3] {
        N.with(Cell::get)
    }
}

/// Coefficients at and above M are zero.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Zq<F: Modulus> {
    pub c: [u64; CAP],
    f: PhantomData<F>,
}

impl<F: Modulus> Zq<F> {
    pub const ZERO: Self = Self::scalar(0);
    pub const ONE: Self = Self::scalar(1);

    pub const fn scalar(v: u64) -> Self {
        let mut c = [0; CAP];
        c[0] = v;
        Self { c, f: PhantomData }
    }

    /// The lift with 0/1 coefficients of a field element (bit i = z^i).
    pub fn from_bits(v: u128) -> Self {
        assert!(F::M < 128 && v >> F::M == 0);
        let mut c = [0; CAP];
        for (i, x) in c.iter_mut().enumerate().take(F::M) {
            *x = (v >> i & 1) as u64;
        }
        Self { c, f: PhantomData }
    }

    /// Reduction mod 2.
    pub fn bits(&self) -> u128 {
        (0..F::M).fold(0, |acc, i| acc | ((self.c[i] & 1) as u128) << i)
    }

    fn map(&self, f: impl Fn(u64) -> u64) -> Self {
        let mut c = [0; CAP];
        for (y, &x) in c.iter_mut().zip(&self.c[..F::M]) {
            *y = f(x);
        }
        Self { c, f: PhantomData }
    }

    pub fn scale(&self, k: u64) -> Self {
        self.map(|x| x.wrapping_mul(k))
    }

    /// Times 2^n.
    pub fn shl(&self, n: u32) -> Self {
        self.map(|x| x << n)
    }

    /// Exact division of a multiple of 2^n: the quotient is right mod 2^(64-n).
    pub fn shr(&self, n: u32) -> Self {
        self.map(|x| x >> n)
    }

    pub fn square(&self) -> Self {
        #[cfg(test)]
        ops::bump(ops::SQR);
        let m = F::M;
        let a = &self.c;
        let mut p = [0u64; 2 * CAP];
        for i in 0..m {
            for j in i + 1..m {
                p[i + j] = p[i + j].wrapping_add(a[i].wrapping_mul(a[j]));
            }
        }
        for (i, &x) in a.iter().enumerate().take(m) {
            p[2 * i] = (p[2 * i] << 1).wrapping_add(x.wrapping_mul(x));
            p[2 * i + 1] <<= 1;
        }
        reduce(p)
    }

    pub fn pow(&self, mut e: usize) -> Self {
        let (mut r, mut b) = (Self::ONE, *self);
        while e > 0 {
            if e & 1 == 1 {
                r = r * b;
            }
            b = b.square();
            e >>= 1;
        }
        r
    }

    /// Inverse of a unit: extended Euclid on its nonzero reduction in GF(2^m),
    /// then Newton lifting.
    pub fn inv(&self) -> Self {
        let mut y = Self::from_bits(gf2_inv(self.bits(), modulus_bits::<F>()));
        for _ in 0..6 {
            y = y * (Self::scalar(2) - *self * y);
        }
        y
    }
}

/// Fold t^i for i >= M down with t^M = -sum t^k, from the top, so that
/// every fold lands below the coefficient being folded.
fn reduce<F: Modulus>(mut p: [u64; 2 * CAP]) -> Zq<F> {
    let m = F::M;
    for i in (m..2 * m - 1).rev() {
        let h = p[i];
        for &k in F::LOW {
            p[i - m + k] = p[i - m + k].wrapping_sub(h);
        }
    }
    let mut c = [0; CAP];
    c[..m].copy_from_slice(&p[..m]);
    Zq { c, f: PhantomData }
}

impl<F: Modulus> Mul for Zq<F> {
    type Output = Self;
    /// Schoolbook: m^2 word multiplications.
    fn mul(self, o: Self) -> Self {
        #[cfg(test)]
        ops::bump(ops::MUL);
        let m = F::M;
        let mut p = [0u64; 2 * CAP];
        for i in 0..m {
            let ai = self.c[i];
            for j in 0..m {
                p[i + j] = p[i + j].wrapping_add(ai.wrapping_mul(o.c[j]));
            }
        }
        reduce(p)
    }
}

impl<F: Modulus> Add for Zq<F> {
    type Output = Self;
    fn add(mut self, o: Self) -> Self {
        for i in 0..F::M {
            self.c[i] = self.c[i].wrapping_add(o.c[i]);
        }
        self
    }
}

impl<F: Modulus> Sub for Zq<F> {
    type Output = Self;
    fn sub(mut self, o: Self) -> Self {
        for i in 0..F::M {
            self.c[i] = self.c[i].wrapping_sub(o.c[i]);
        }
        self
    }
}

impl<F: Modulus> Neg for Zq<F> {
    type Output = Self;
    fn neg(self) -> Self {
        self.map(u64::wrapping_neg)
    }
}

/// The reduced modulus z^M + sum z^k as bits; M <= 127.
fn modulus_bits<F: Modulus>() -> u128 {
    F::LOW.iter().fold(1 << F::M, |acc, &k| acc | 1 << k)
}

/// a^-1 mod (2, f) for a != 0 of lower degree than f, by Euclid on bit vectors.
fn gf2_inv(a: u128, f: u128) -> u128 {
    assert!(a != 0);
    let deg = |v: u128| 127 - v.leading_zeros() as i32;
    let (mut u, mut v, mut g1, mut g2) = (a, f, 1u128, 0u128);
    while u != 1 {
        // u = 0 once gcd(a, f) != 1 has been reached
        assert!(u != 0, "not a unit: is the modulus irreducible?");
        let mut j = deg(u) - deg(v);
        if j < 0 {
            (u, v, g1, g2) = (v, u, g2, g1);
            j = -j;
        }
        u ^= v << j;
        g1 ^= g2 << j;
    }
    g1
}

/// sigma and sigma^-1 as matrices (column i = the image of t^i), and the
/// trace form.
pub struct Frobenius<F: Modulus> {
    fwd: Vec<Zq<F>>,
    inv: Vec<Zq<F>>,
    /// Tr(t^i), the power sums of the roots of F.
    power_sums: [u64; CAP],
}

impl<F: Modulus> Default for Frobenius<F> {
    fn default() -> Self {
        Self::new()
    }
}

impl<F: Modulus> Frobenius<F> {
    /// About 750 multiplication-sized operations: once per field, not per
    /// curve.
    pub fn new() -> Self {
        let m = F::M;
        let t = Zq::<F>::from_bits(2);
        // sigma(t) is the root of F above t^2; Newton, 1 bit -> 64 bits.
        let mut theta = t.square();
        for _ in 0..6 {
            let (f, df) = eval_modulus(&theta);
            theta = theta - f * df.inv();
        }
        let fwd = powers(theta);
        let mut s = t;
        for _ in 1..m {
            s = apply(&fwd, &s);
        }
        // Reducible but squarefree moduli reach this point; for prime m, sigma
        // then has order below m.
        assert!(
            apply(&fwd, &s) == t,
            "sigma^m != 1: modulus not irreducible"
        );
        let inv = powers(s);
        // Newton's identities: p_i = -(sum_{j<i} c_j p_{i-j}) - i c_i, with
        // F = t^M + sum_j c_j t^(M-j).
        let c = |j: usize| F::LOW.contains(&(m - j)) as u64;
        let mut power_sums = [0u64; CAP];
        power_sums[0] = m as u64;
        for i in 1..m {
            let mut s = (i as u64).wrapping_mul(c(i));
            for j in 1..i {
                s = s.wrapping_add(c(j).wrapping_mul(power_sums[i - j]));
            }
            power_sums[i] = s.wrapping_neg();
        }
        Self {
            fwd,
            inv,
            power_sums,
        }
    }

    /// m^2 word multiplications, like a schoolbook product.
    pub fn sigma(&self, x: &Zq<F>) -> Zq<F> {
        apply(&self.fwd, x)
    }

    pub fn sigma_inv(&self, x: &Zq<F>) -> Zq<F> {
        apply(&self.inv, x)
    }

    pub fn trace(&self, x: &Zq<F>) -> u64 {
        (0..F::M).fold(0u64, |acc, i| {
            acc.wrapping_add(x.c[i].wrapping_mul(self.power_sums[i]))
        })
    }
}

/// (F(r), F'(r))
fn eval_modulus<F: Modulus>(r: &Zq<F>) -> (Zq<F>, Zq<F>) {
    let m = F::M;
    let mut f = r.pow(m);
    let mut df = r.pow(m - 1).scale(m as u64);
    for &k in F::LOW {
        f = f + r.pow(k);
        if k > 0 {
            df = df + r.pow(k - 1).scale(k as u64);
        }
    }
    (f, df)
}

fn powers<F: Modulus>(r: Zq<F>) -> Vec<Zq<F>> {
    let mut v = vec![Zq::ONE];
    for i in 1..F::M {
        v.push(v[i - 1] * r);
    }
    v
}

fn apply<F: Modulus>(cols: &[Zq<F>], x: &Zq<F>) -> Zq<F> {
    #[cfg(test)]
    ops::bump(ops::FROB);
    let m = F::M;
    let mut out = Zq::ZERO;
    for (col, &xi) in cols.iter().zip(&x.c[..m]) {
        for j in 0..m {
            out.c[j] = out.c[j].wrapping_add(xi.wrapping_mul(col.c[j]));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::gf2_127::{MASK127, from_u128, tests::mul_ref};
    use proptest::prelude::*;
    use std::sync::OnceLock;

    type Z = Zq<Gf127>;

    fn frob() -> &'static Frobenius<Gf127> {
        static F: OnceLock<Frobenius<Gf127>> = OnceLock::new();
        F.get_or_init(Frobenius::new)
    }

    fn elt() -> impl Strategy<Value = Z> {
        any::<[u64; 127]>().prop_map(|v| {
            let mut c = [0; CAP];
            c[..127].copy_from_slice(&v);
            Zq { c, f: PhantomData }
        })
    }

    /// Nonzero mod 2, hence a unit in the unramified extension.
    fn unit() -> impl Strategy<Value = Z> {
        (elt(), 1..=MASK127).prop_map(|(x, r)| x - Z::from_bits(x.bits()) + Z::from_bits(r))
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(32))]

        #[test]
        fn reduces_to_gf2(a in elt(), b in elt()) {
            prop_assert_eq!((a * b).bits(), mul_ref(a.bits(), b.bits()));
            prop_assert_eq!(a.square(), a * a);
            prop_assert_eq!((a + b).bits(), a.bits() ^ b.bits());
        }

        #[test]
        fn ring_laws(a in elt(), b in elt(), c in elt()) {
            prop_assert_eq!(a * b, b * a);
            prop_assert_eq!((a * b) * c, a * (b * c));
            prop_assert_eq!(a * (b + c), a * b + a * c);
            prop_assert_eq!(a - b + b, a);
            prop_assert_eq!(a + -a, Z::ZERO);
            prop_assert_eq!(a.shl(3), a.scale(8));
            prop_assert_eq!(a.shl(5).shr(5).shl(5), a.shl(5));
        }

        #[test]
        fn inverse(a in unit()) {
            prop_assert_eq!(a * a.inv(), Z::ONE);
        }

        #[test]
        fn frobenius(a in elt(), b in elt()) {
            let f = frob();
            prop_assert_eq!(f.sigma(&(a * b)), f.sigma(&a) * f.sigma(&b));
            prop_assert_eq!(f.sigma(&(a + b)), f.sigma(&a) + f.sigma(&b));
            prop_assert_eq!(f.sigma(&a).bits(), a.square().bits());
            prop_assert_eq!(f.sigma_inv(&f.sigma(&a)), a);
            prop_assert_eq!(f.trace(&f.sigma(&a)), f.trace(&a));
        }

        #[test]
        fn trace_reduces_to_gf2(a in elt()) {
            prop_assert_eq!(frob().trace(&a) & 1, from_u128(a.bits()).trace() as u64);
        }
    }

    /// z^127 + z^63 + z = z (z^63 + z^31 + 1)^2: a repeated-factor modulus
    /// for testing rejection of nonunits without looping.
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    struct Square;

    impl Modulus for Square {
        const M: usize = 127;
        const LOW: &'static [usize] = &[63, 1];
    }

    /// z^127 + z^62 + 1, squarefree with factors of degree 2, 27, 28, 32, 38.
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    struct Squarefree;

    impl Modulus for Squarefree {
        const M: usize = 127;
        const LOW: &'static [usize] = &[62, 0];
    }

    #[test]
    #[should_panic(expected = "not a unit")]
    fn modulus_with_square_factor_panics() {
        Frobenius::<Square>::new();
    }

    #[test]
    #[should_panic(expected = "sigma^m != 1")]
    fn squarefree_reducible_modulus_panics() {
        Frobenius::<Squarefree>::new();
    }

    #[test]
    fn trace_of_one_and_sigma_has_order_m() {
        let f = frob();
        assert_eq!(f.trace(&Z::ONE), 127);
        let t = Z::from_bits(2);
        let mut s = t;
        for _ in 0..127 {
            s = f.sigma(&s);
        }
        assert_eq!(s, t);
    }
}
