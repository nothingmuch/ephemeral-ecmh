//! Dense univariate polynomials over GF(2^127), F_p, p = 2^127 - 1, and the
//! quadratic fields GF(p^2), for the torsion sieve (`sieve`): products,
//! division, gcd, arithmetic modulo a fixed polynomial, and the roots of a
//! polynomial in the field.
//!
//! Everything is schoolbook. The sieve's moduli are division polynomials of
//! degree d = (l^2 - 1)/2, and its time goes to d-by-d products modulo them
//! (`Modulus::square`). Schoolbook multiplication is quadratic in d; the
//! concrete crossover for Karatsuba products with Barrett reduction over
//! F_p has not been measured. Over GF(2^127) squaring mod m is linear and
//! costs d^2/2 M (see `sieve` for the cost model).

use crate::field::gf2_127::{Gf, from_u128};
use core::ops::{Add, Mul, Neg, Sub};

/// A coefficient field: `field::Field`, and what root finding needs.
pub trait Field: crate::field::Field + Sub<Output = Self> + Neg<Output = Self> {
    /// Squaring is additive in characteristic 2: a polynomial squares
    /// coefficient by coefficient, and a -> a^2 mod m is semilinear.
    const CHAR2: bool;
    /// q, the number of elements.
    const Q: u128;
    /// n mod the characteristic.
    fn small(n: u64) -> Self;
    /// The k-th shift root finding tries: z^k in GF(2^127), a basis, so
    /// some k < 127 separates any two roots; k in F_p; j + c i in GF(p^2),
    /// k = 3j + c. A line of shifts j + c i, j in F_p, fails to separate
    /// roots r and s whose norms N(r + j + c i) and N(s + j + c i) are
    /// both squares (Im r = Im s = -c) or equal (Re r = Re s, Im r + Im s
    /// = -2c), so for at most two c: conjugates for c = 0, for one.
    fn shift(k: u32) -> Self;
}

impl Field for Gf {
    const CHAR2: bool = true;
    const Q: u128 = 1 << 127;
    fn small(n: u64) -> Self {
        if n & 1 == 1 { Gf::ONE } else { Gf::ZERO }
    }
    fn shift(k: u32) -> Self {
        assert!(k < 127, "z^0..z^126 separate any two roots");
        from_u128(1 << k)
    }
}

/// Coefficients from x^0 up, with no trailing zeros (0 is empty).
#[derive(Clone, Debug)]
pub struct Poly<F>(Vec<F>);

impl<F: Field> PartialEq for Poly<F> {
    fn eq(&self, o: &Self) -> bool {
        self.0.len() == o.0.len() && self.0.iter().zip(&o.0).all(|(&a, &b)| (a - b).is_zero())
    }
}

impl<F: Field> Poly<F> {
    pub fn new(mut c: Vec<F>) -> Self {
        while c.last().is_some_and(|a| a.is_zero()) {
            c.pop();
        }
        Self(c)
    }

    pub fn zero() -> Self {
        Self(Vec::new())
    }

    pub fn constant(c: F) -> Self {
        Self::new(vec![c])
    }

    /// x^k
    pub fn monomial(k: usize) -> Self {
        let mut c = vec![F::ZERO; k + 1];
        c[k] = F::ONE;
        Self(c)
    }

    pub fn coeffs(&self) -> &[F] {
        &self.0
    }

    pub fn is_zero(&self) -> bool {
        self.0.is_empty()
    }

    /// `None` for 0.
    pub fn degree(&self) -> Option<usize> {
        self.0.len().checked_sub(1)
    }

    pub fn eval(&self, x: F) -> F {
        self.0.iter().rev().fold(F::ZERO, |acc, &c| acc * x + c)
    }

    pub fn scale(&self, c: F) -> Self {
        Self::new(self.0.iter().map(|&a| a * c).collect())
    }

    pub fn monic(&self) -> Self {
        match self.0.last() {
            Some(&l) => self.scale(l.inv()),
            None => Self::zero(),
        }
    }

    /// d S + d(d-1)/2 M for d coefficients, or d S in characteristic 2.
    pub fn square(&self) -> Self {
        let a = &self.0;
        if a.is_empty() {
            return Self::zero();
        }
        let mut out = vec![F::ZERO; 2 * a.len() - 1];
        if !F::CHAR2 {
            for (i, &x) in a.iter().enumerate() {
                for (o, &y) in out[2 * i + 1..].iter_mut().zip(&a[i + 1..]) {
                    *o = *o + x * y;
                }
            }
            for o in out.iter_mut() {
                *o = *o + *o;
            }
        }
        for (i, &x) in a.iter().enumerate() {
            out[2 * i] = out[2 * i] + x.square();
        }
        Self::new(out)
    }

    /// (q, r) with self = q d + r and deg r < deg d; d must be nonzero.
    pub fn divrem(&self, d: &Self) -> (Self, Self) {
        let n = d.0.len();
        assert!(n > 0, "division by zero");
        if self.0.len() < n {
            return (Self::zero(), self.clone());
        }
        let inv = d.0[n - 1].inv();
        let mut r = self.0.clone();
        let mut q = vec![F::ZERO; r.len() - n + 1];
        for i in (0..q.len()).rev() {
            let c = r[i + n - 1] * inv;
            q[i] = c;
            for (x, &y) in r[i..i + n - 1].iter_mut().zip(&d.0) {
                *x = *x - c * y;
            }
        }
        r.truncate(n - 1);
        (Self::new(q), Self::new(r))
    }

    /// Monic gcd (0 if both are 0). Remainders only matter up to a scalar,
    /// so each step cancels a leading term by cross-multiplying: one
    /// inversion in all, for the final monic, rather than one per step
    /// (the remainders' degrees drop by one per step; the cost of an
    /// inversion relative to M depends on the field and implementation).
    pub fn gcd(&self, o: &Self) -> Self {
        let (mut a, mut b) = (self.0.clone(), o.0.clone());
        while let Some(&lb) = b.last() {
            // a <- lb a - la x^s b until deg a < deg b
            while a.len() >= b.len() {
                let la = a.pop().unwrap();
                let s = a.len() + 1 - b.len();
                for x in &mut a[..s] {
                    *x = *x * lb;
                }
                for (x, &y) in a[s..].iter_mut().zip(&b) {
                    *x = *x * lb - la * y;
                }
                while a.last().is_some_and(|c| c.is_zero()) {
                    a.pop();
                }
            }
            (a, b) = (b, a);
        }
        Self(a).monic()
    }
}

impl<F: Field> Add for &Poly<F> {
    type Output = Poly<F>;
    fn add(self, o: &Poly<F>) -> Poly<F> {
        let (long, short) = if self.0.len() >= o.0.len() {
            (self, o)
        } else {
            (o, self)
        };
        let mut c = long.0.clone();
        for (x, &y) in c.iter_mut().zip(&short.0) {
            *x = *x + y;
        }
        Poly::new(c)
    }
}

impl<F: Field> Neg for &Poly<F> {
    type Output = Poly<F>;
    fn neg(self) -> Poly<F> {
        Poly(self.0.iter().map(|&a| -a).collect())
    }
}

impl<F: Field> Sub for &Poly<F> {
    type Output = Poly<F>;
    fn sub(self, o: &Poly<F>) -> Poly<F> {
        self + &-o
    }
}

impl<F: Field> Mul for &Poly<F> {
    type Output = Poly<F>;
    fn mul(self, o: &Poly<F>) -> Poly<F> {
        let (a, b) = (&self.0, &o.0);
        if a.is_empty() || b.is_empty() {
            return Poly::zero();
        }
        let mut out = vec![F::ZERO; a.len() + b.len() - 1];
        for (i, &x) in a.iter().enumerate() {
            for (o, &y) in out[i..].iter_mut().zip(b) {
                *o = *o + x * y;
            }
        }
        Poly::new(out)
    }
}

/// Arithmetic modulo a fixed m of degree d >= 1, kept monic.
pub struct Modulus<F> {
    m: Vec<F>,
    /// Characteristic 2: x^(2i) mod m for ceil(d/2) <= i < d, each padded
    /// to d coefficients. Squaring is then a sum of rows, d S + d^2/2 M,
    /// where squaring and reducing would cost d S + d^2 M.
    rows: Vec<Vec<F>>,
}

impl<F: Field> Modulus<F> {
    pub fn new(m: &Poly<F>) -> Self {
        assert!(m.degree().is_some_and(|d| d >= 1), "modulus of degree 0");
        let mut md = Self {
            m: m.monic().0,
            rows: Vec::new(),
        };
        if F::CHAR2 {
            let d = md.degree();
            let h = d.div_ceil(2);
            let mut row = md.rem(Poly::monomial(2 * h)).0;
            for _ in h..d {
                row.resize(d, F::ZERO);
                md.rows.push(row.clone());
                row.splice(0..0, [F::ZERO; 2]);
                md.reduce(&mut row);
            }
        }
        md
    }

    pub fn degree(&self) -> usize {
        self.m.len() - 1
    }

    pub fn modulus(&self) -> Poly<F> {
        Poly(self.m.clone())
    }

    /// Long division by the monic m in place: d M per coefficient above d.
    fn reduce(&self, v: &mut Vec<F>) {
        let d = self.degree();
        for i in (d..v.len()).rev() {
            let c = v[i];
            if c.is_zero() {
                continue;
            }
            for (x, &y) in v[i - d..i].iter_mut().zip(&self.m) {
                *x = *x - c * y;
            }
        }
        v.truncate(d);
        while v.last().is_some_and(|a| a.is_zero()) {
            v.pop();
        }
    }

    pub fn rem(&self, a: Poly<F>) -> Poly<F> {
        let mut v = a.0;
        self.reduce(&mut v);
        Poly(v)
    }

    pub fn mul(&self, a: &Poly<F>, b: &Poly<F>) -> Poly<F> {
        self.rem(a * b)
    }

    /// a^2 mod m, for a already reduced.
    pub fn square(&self, a: &Poly<F>) -> Poly<F> {
        if !F::CHAR2 {
            return self.rem(a.square());
        }
        let d = self.degree();
        let h = d.div_ceil(2);
        let mut out = vec![F::ZERO; d];
        for (i, &c) in a.0.iter().enumerate() {
            let s = c.square();
            if i < h {
                out[2 * i] = out[2 * i] + s;
            } else {
                for (o, &r) in out.iter_mut().zip(&self.rows[i - h]) {
                    *o = *o + s * r;
                }
            }
        }
        Poly::new(out)
    }

    /// a^e mod m for e >= 1, a already reduced: a squaring per bit of e
    /// after the first, and a product by a per set bit. `mul_a` is that
    /// product, which for a = x is a shift, d M.
    fn pow(&self, a: &Poly<F>, e: u128, mul_a: impl Fn(&Poly<F>) -> Poly<F>) -> Poly<F> {
        assert!(e >= 1);
        (0..127 - e.leading_zeros())
            .rev()
            .fold(a.clone(), |acc, i| {
                let acc = self.square(&acc);
                if e >> i & 1 == 1 { mul_a(&acc) } else { acc }
            })
    }

    /// x^e mod m, e >= 1.
    pub fn pow_x(&self, e: u128) -> Poly<F> {
        let shift = |a: &Poly<F>| {
            let mut v = Vec::with_capacity(a.0.len() + 1);
            v.push(F::ZERO);
            v.extend_from_slice(&a.0);
            self.rem(Poly(v))
        };
        self.pow(&self.rem(Poly::monomial(1)), e, shift)
    }
}

/// The distinct roots of a squarefree m in F, as the monic gcd(m, x^q - x),
/// or gcd(m, x^(q+1) - x^2), the same roots (0 too), when q + 1 has fewer
/// set bits. Over GF(2^127) and F_p, p = 2^127 - 1, that leaves x^(2^127),
/// squarings alone; over GF(p^2) each set bit adds a product by x, d M
/// against a squaring's d^2.
pub fn field_roots<F: Field>(m: &Poly<F>) -> Poly<F> {
    if m.degree() == Some(0) {
        return Poly::constant(F::ONE);
    }
    let md = Modulus::new(m);
    let (e, k) = match F::Q.checked_add(1) {
        Some(e) if e.count_ones() < F::Q.count_ones() => (e, 2),
        _ => (F::Q, 1),
    };
    m.gcd(&(&md.pow_x(e) - &md.rem(Poly::monomial(k))))
}

/// Takes the value 0 at the roots x0 of m on one side of a split by delta:
/// Tr(delta x0) = 0 over GF(2^127), x0 + delta a nonzero square in an odd
/// field.
fn splitter<F: Field>(md: &Modulus<F>, delta: F) -> Poly<F> {
    if F::CHAR2 {
        // sum of (delta x)^(2^i), i < 127
        let mut t = md.rem(Poly::new(vec![F::ZERO, delta]));
        let mut acc = t.clone();
        for _ in 1..F::Q.trailing_zeros() {
            t = md.square(&t);
            acc = &acc + &t;
        }
        acc
    } else {
        // (x + delta)^((q-1)/2) - 1
        let a = md.rem(Poly::new(vec![delta, F::ONE]));
        let r = md.pow(&a, (F::Q - 1) / 2, |r| md.mul(r, &a));
        &r - &Poly::constant(F::ONE)
    }
}

/// Two proper factors of m, which has degree >= 2 and distinct roots in F.
fn split<F: Field>(m: &Poly<F>) -> (Poly<F>, Poly<F>) {
    let md = Modulus::new(m);
    let d = md.degree();
    for k in 0.. {
        let g = m.gcd(&splitter(&md, F::shift(k)));
        if g.degree().is_some_and(|e| 0 < e && e < d) {
            let h = m.divrem(&g).0;
            return (g, h);
        }
    }
    unreachable!()
}

/// Depth-first over the roots of m, which must have distinct roots all in
/// F (as from `field_roots`), smaller factors first: the first root that
/// `test` accepts.
pub fn find_root<F: Field, W>(m: &Poly<F>, test: &mut impl FnMut(F) -> Option<W>) -> Option<W> {
    match m.degree() {
        None | Some(0) => None,
        Some(1) => test(-(m.0[0] * m.0[1].inv())),
        Some(_) => {
            let (a, b) = split(m);
            let (a, b) = if a.0.len() <= b.0.len() {
                (a, b)
            } else {
                (b, a)
            };
            find_root(&a, test).or_else(|| find_root(&b, test))
        }
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;
    use crate::field::batch::Invert;
    use proptest::prelude::*;

    fn gf() -> impl Strategy<Value = Gf> {
        crate::field::gf2_127::tests::fe().prop_map(from_u128)
    }

    fn poly<F: Field + core::fmt::Debug>(
        e: impl Strategy<Value = F>,
        max: usize,
    ) -> impl Strategy<Value = Poly<F>> {
        prop::collection::vec(e, 0..max).prop_map(Poly::new)
    }

    /// Monic, of degree 1..max.
    fn modulus<F: Field + core::fmt::Debug>(
        e: impl Strategy<Value = F>,
        max: usize,
    ) -> impl Strategy<Value = Poly<F>> {
        prop::collection::vec(e, 1..max).prop_map(|mut c| {
            c.push(F::ONE);
            Poly::new(c)
        })
    }

    fn linear<F: Field>(r: F) -> Poly<F> {
        Poly::new(vec![-r, F::ONE])
    }

    fn product<F: Field>(rs: &[F]) -> Poly<F> {
        rs.iter()
            .fold(Poly::constant(F::ONE), |acc, &r| &acc * &linear(r))
    }

    /// Distinct roots from a list with repeats.
    fn dedup<F: Field>(rs: &[F]) -> Vec<F> {
        let mut out: Vec<F> = Vec::new();
        for &r in rs {
            if !out.iter().any(|&s| (r - s).is_zero()) {
                out.push(r);
            }
        }
        out
    }

    fn ring_laws<F: Field + core::fmt::Debug>(a: &Poly<F>, b: &Poly<F>, c: &Poly<F>) {
        assert_eq!(&(a * b) * c, a * &(b * c));
        assert_eq!(a * &(b + c), &(a * b) + &(a * c));
        assert_eq!(a * b, b * a);
        assert_eq!(a.square(), a * a);
        assert_eq!(&(a - b) + b, *a);
        let x = F::shift(5);
        assert!(((a * b).eval(x) - a.eval(x) * b.eval(x)).is_zero());
    }

    fn division<F: Field + core::fmt::Debug>(a: &Poly<F>, d: &Poly<F>) {
        if d.is_zero() {
            return;
        }
        let (q, r) = a.divrem(d);
        assert_eq!(&(&q * d) + &r, *a);
        assert!(r.0.len() < d.0.len());
        assert_eq!(q.0.len(), (a.0.len() + 1).saturating_sub(d.0.len()));
    }

    fn modular<F: Field + core::fmt::Debug>(m: &Poly<F>, a: &Poly<F>, b: &Poly<F>) {
        let md = Modulus::new(m);
        let (a, b) = (md.rem(a.clone()), md.rem(b.clone()));
        assert_eq!(a, a.divrem(m).1);
        assert_eq!(md.mul(&a, &b), (&a * &b).divrem(m).1);
        assert_eq!(md.square(&a), a.square().divrem(m).1);
        // a^11 = ((a^2)^2 a)^2 a, and x^11
        let mut want = a.clone();
        for bit in [false, true, true] {
            want = (&want * &want).divrem(m).1;
            if bit {
                want = (&want * &a).divrem(m).1;
            }
        }
        assert_eq!(md.pow(&a, 11, |r| md.mul(r, &a)), want);
        assert_eq!(md.pow_x(11), Poly::monomial(11).divrem(m).1);
    }

    /// gcd(c (x - r0)...(x - rk), (x - r0)(x - t)) = x - r0 if t is no r_i.
    fn gcd_of_products<F: Field + core::fmt::Debug>(rs: &[F], t: F, c: F) {
        if c.is_zero() || rs.iter().any(|&r| (r - t).is_zero()) {
            return;
        }
        let a = product(rs).scale(c);
        let b = &linear(rs[0]) * &linear(t);
        assert_eq!(a.gcd(&b), linear(rs[0]));
        assert_eq!(b.gcd(&a), linear(rs[0]));
        assert_eq!(a.gcd(&Poly::zero()), a.monic());
    }

    /// Roots rs times a factor with no roots in F: field_roots and
    /// find_root recover exactly the distinct rs.
    fn roots<F: Field + core::fmt::Debug>(rs: &[F], irreducible: &Poly<F>) {
        let want = dedup(rs);
        let m = &product(&want) * irreducible;
        let g = field_roots(&m);
        assert_eq!(g, product(&want));
        let mut seen = Vec::new();
        let none: Option<()> = find_root(&g, &mut |r| {
            seen.push(r);
            None
        });
        assert!(none.is_none());
        assert_eq!(seen.len(), want.len());
        assert!(
            want.iter()
                .all(|&r| seen.iter().any(|&s| (r - s).is_zero()))
        );
        if let Some(&r0) = want.last() {
            let hit = find_root(&g, &mut |r| (r - r0).is_zero().then_some(r));
            assert!(hit.is_some_and(|r| (r - r0).is_zero()));
        }
    }

    proptest! {
        #[test]
        fn gf2_127_ring_laws(a in poly(gf(), 10), b in poly(gf(), 10), c in poly(gf(), 10)) {
            ring_laws(&a, &b, &c);
        }

        #[test]
        fn gf2_127_division(a in poly(gf(), 16), d in poly(gf(), 8)) {
            division(&a, &d);
        }

        #[test]
        fn gf2_127_modular(m in modulus(gf(), 12), a in poly(gf(), 16), b in poly(gf(), 16)) {
            modular(&m, &a, &b);
        }

        #[test]
        fn gf2_127_gcd_of_products(rs in prop::collection::vec(gf(), 1..6), t in gf(), c in gf()) {
            gcd_of_products(&rs, t, c);
        }

        #[test]
        fn gf2_127_roots(rs in prop::collection::vec(gf(), 0..8), c in gf()) {
            // x^2 + x + c has no root in GF(2^127) iff Tr(c) = 1
            let c = if c.trace() == 1 { c } else { c + Gf::ONE };
            roots(&rs, &Poly::new(vec![c, Gf::ONE, Gf::ONE]));
        }
    }

    #[test]
    fn roots_of_split_moduli() {
        // every root in F, 0 among them (over F_p, x^(p+1) - x^2 has it too)
        let gf: Vec<Gf> = (0..6).map(|k| from_u128(k * 0x1234_5678_9abc)).collect();
        roots(&gf, &Poly::constant(Gf::ONE));
    }

    /// Field multiplications, squarings and inversions, the units of the
    /// sieve's cost model (`sieve`). Additions are not counted; the model
    /// treats them as negligible beside these.
    #[derive(Clone, Copy, Debug, Default, PartialEq)]
    pub(crate) struct Ops {
        pub m: u64,
        pub s: u64,
        pub i: u64,
    }

    thread_local! {
        static OPS: core::cell::Cell<Ops> = Default::default();
    }

    fn tally(op: fn(&mut Ops) -> &mut u64) {
        OPS.with(|c| {
            let mut o = c.get();
            *op(&mut o) += 1;
            c.set(o);
        });
    }

    /// The ops f performs on `Counted` elements, on this thread.
    pub(crate) fn count<R>(f: impl FnOnce() -> R) -> (R, Ops) {
        let a = OPS.with(|c| c.get());
        let r = f();
        let b = OPS.with(|c| c.get());
        let ops = Ops {
            m: b.m - a.m,
            s: b.s - a.s,
            i: b.i - a.i,
        };
        (r, ops)
    }

    /// F, with its multiplications, squarings and inversions counted.
    #[derive(Clone, Copy, Debug)]
    pub(crate) struct Counted<F>(pub F);

    impl<F: Field> Add for Counted<F> {
        type Output = Self;
        fn add(self, o: Self) -> Self {
            Counted(self.0 + o.0)
        }
    }

    impl<F: Field> Sub for Counted<F> {
        type Output = Self;
        fn sub(self, o: Self) -> Self {
            Counted(self.0 - o.0)
        }
    }

    impl<F: Field> Neg for Counted<F> {
        type Output = Self;
        fn neg(self) -> Self {
            Counted(-self.0)
        }
    }

    impl<F: Field> Mul for Counted<F> {
        type Output = Self;
        fn mul(self, o: Self) -> Self {
            tally(|t| &mut t.m);
            Counted(self.0 * o.0)
        }
    }

    impl<F: Field> Invert for Counted<F> {
        const ONE: Self = Counted(F::ONE);
        fn inv(self) -> Self {
            tally(|t| &mut t.i);
            Counted(self.0.inv())
        }
    }

    impl<F: Field> crate::field::Field for Counted<F> {
        const ZERO: Self = Counted(F::ZERO);
        fn is_zero(self) -> bool {
            self.0.is_zero()
        }
        fn equals(self, o: Self) -> bool {
            self.0.equals(o.0)
        }
        fn square(self) -> Self {
            tally(|t| &mut t.s);
            Counted(self.0.square())
        }
    }

    impl<F: Field> Field for Counted<F> {
        const CHAR2: bool = F::CHAR2;
        const Q: u128 = F::Q;
        fn small(n: u64) -> Self {
            Counted(F::small(n))
        }
        fn shift(k: u32) -> Self {
            Counted(F::shift(k))
        }
    }
}
