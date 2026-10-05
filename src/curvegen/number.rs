//! Integer arithmetic used by curve-order certificate checks.
//!
//! Montgomery residues accelerate repeated modular products; primality
//! acceptance uses the same fixed-base probable-prime test in every selector.

/// The embedding degree of an accepted r must exceed this: r divides no
/// q^k - 1 with k <= 2^20, so the pairing reductions of Menezes–Okamoto–
/// Vanstone and Frey–Rück lead into an extension GF(q^k) of more than
/// 2^20 times the field's bits. One bound for every family.
pub const EMBEDDING_MIN: u32 = 1 << 20;

/// a*b mod m for any m > 0. Prime-order r can exceed 2^127, so the
/// doubling and adding must not overflow: x + y mod m as x - (m - y).
fn mulmod(a: u128, b: u128, m: u128) -> u128 {
    let add = |x: u128, y: u128| if x >= m - y { x - (m - y) } else { x + y };
    let (a, mut r) = (a % m, 0u128);
    for i in (0..128).rev() {
        r = add(r, r);
        if b >> i & 1 == 1 {
            r = add(r, a);
        }
    }
    r
}

/// Arithmetic mod an odd m in Montgomery form, R = 2^128. A product is
/// eight 64-bit multiplies (the 256-bit product, then the reduction's),
/// against `mulmod`'s 128 dependent add-and-compare steps. Certificate
/// checks spend ~4,500 products in `is_prime` and about 2^11 in
/// `embedding_degree_ok`.
pub(crate) struct Mont {
    m: u128,
    /// -m^-1 mod R
    m_neg_inv: u128,
    /// R^2 mod m, to enter the form
    r2: u128,
    /// R mod m, the form's 1
    pub(crate) one: u128,
}

/// The 256-bit product as (low, high) halves.
pub(crate) fn mul_wide(a: u128, b: u128) -> (u128, u128) {
    let (a0, a1) = (a as u64 as u128, a >> 64);
    let (b0, b1) = (b as u64 as u128, b >> 64);
    let (mid, c1) = (a0 * b1).overflowing_add(a1 * b0);
    let (lo, c2) = (a0 * b0).overflowing_add(mid << 64);
    let hi = a1 * b1 + (mid >> 64) + ((c1 as u128) << 64) + c2 as u128;
    (lo, hi)
}

impl Mont {
    pub(crate) fn new(m: u128) -> Self {
        assert!(m & 1 == 1);
        // Newton's iteration doubles the correct low bits: 3 from m itself
        // (m * m = 1 mod 8), then 6, 12, ..., 192
        let mut inv = m;
        for _ in 0..6 {
            inv = inv.wrapping_mul(2u128.wrapping_sub(m.wrapping_mul(inv)));
        }
        let one = m.wrapping_neg() % m;
        Mont {
            m,
            m_neg_inv: inv.wrapping_neg(),
            r2: mulmod(one, one, m),
            one,
        }
    }

    /// (lo + hi R) / R mod m, for lo + hi R < m R. The sum before the
    /// final subtraction is below 2m, which can exceed 2^128.
    fn redc(&self, lo: u128, hi: u128) -> u128 {
        let (_, uhi) = mul_wide(lo.wrapping_mul(self.m_neg_inv), self.m);
        // lo + (u m mod R) is 0 or R: it carries iff lo != 0
        let (t, o1) = hi.overflowing_add(uhi);
        let (t, o2) = t.overflowing_add((lo != 0) as u128);
        if o1 || o2 || t >= self.m {
            t.wrapping_sub(self.m)
        } else {
            t
        }
    }

    pub(crate) fn mul(&self, a: u128, b: u128) -> u128 {
        let (lo, hi) = mul_wide(a, b);
        self.redc(lo, hi)
    }

    pub(crate) fn to(&self, a: u128) -> u128 {
        self.mul(a % self.m, self.r2)
    }

    fn pow(&self, mut a: u128, mut e: u128) -> u128 {
        let mut r = self.one;
        while e > 0 {
            if e & 1 == 1 {
                r = self.mul(r, a);
            }
            a = self.mul(a, a);
            e >>= 1;
        }
        r
    }
}

/// Miller–Rabin to the first 24 prime bases: a probable-prime test.
///
/// It's proven exact only below ψ₁₃ ≈ 3.3·10²⁴ (≈ 2⁸¹), where the bases up
/// to 41 already are (Sorenson and Webster, 2015). The group orders r it
/// sees are about 2¹²⁶. Composites that pass every fixed base can be
/// constructed. The certificate's r is the prover's; the verifier binds it
/// to the admissible window and to rP = O for the hashed point
/// (`AdmissibleR`, `select::accept`), not to primality, and runs this same
/// test. An accepted certificate therefore shows that r is a probable
/// prime, and what the certificate concludes from r prime holds on that
/// condition.
pub fn is_prime(n: u128) -> bool {
    const BASES: [u128; 24] = [
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89,
    ];
    if n < 2 {
        return false;
    }
    for p in BASES {
        if n.is_multiple_of(p) {
            return n == p;
        }
    }
    let s = (n - 1).trailing_zeros();
    let d = (n - 1) >> s;
    // n is odd here: 2 is among the bases
    let mont = Mont::new(n);
    let minus_one = n - mont.one;
    'base: for a in BASES {
        let mut x = mont.pow(mont.to(a), d);
        if x == mont.one || x == minus_one {
            continue;
        }
        for _ in 1..s {
            x = mont.mul(x, x);
            if x == minus_one {
                continue 'base;
            }
        }
        return false;
    }
    true
}

/// q^k != 1 mod r for every 1 <= k <= bound (bound <= 2^32), by baby-step
/// giant-step with s = ceil(sqrt bound): the baby steps q^i (0 <= i < s),
/// then the giant steps q^(js) (1 <= j <= s). The order k0 of q, if at
/// most s^2, is js - i at the first j where q^(js) = q^i: j = ceil(k0/s)
/// matches, and an earlier match would be a multiple of k0 below it.
/// Certificates call this with an odd probable prime r and a unit q; even
/// or tiny r take the plain loop, and a nonunit q (gcd(q, r) > 1) has no
/// power equal to 1, so the predicate holds. Both cases are exact for any
/// r > 1.
fn embedding_degree_exceeds(q_mod_r: u128, r: u128, bound: u32) -> bool {
    if r > 1 && gcd(q_mod_r, r) != 1 {
        return true;
    }
    if r < 3 || r & 1 == 0 {
        let mut qk = 1;
        for _ in 0..bound {
            qk = mulmod(qk, q_mod_r, r);
            if qk == 1 {
                return false;
            }
        }
        return true;
    }
    let s = (bound as f64).sqrt().ceil() as u64;
    let mont = Mont::new(r);
    let q = mont.to(q_mod_r);
    // Canonical Montgomery residues let the lookup compare integers directly.
    let mut baby = Vec::with_capacity(s as usize);
    let mut x = mont.one;
    for i in 0..s {
        if i > 0 && x == mont.one {
            // the order is i < s, and the baby steps would repeat
            return i > bound.into();
        }
        baby.push((x, i));
        x = mont.mul(x, q);
    }
    baby.sort_unstable();
    let (giant, mut y) = (x, mont.one);
    for j in 1..=s {
        y = mont.mul(y, giant);
        if let Ok(pos) = baby.binary_search_by_key(&y, |&(v, _)| v) {
            return j * s - baby[pos].1 > bound.into();
        }
    }
    true
}

/// The embedding degree of r over F_q exceeds `EMBEDDING_MIN`: q^k != 1 mod
/// r for every 1 <= k <= `EMBEDDING_MIN`, by `embedding_degree_exceeds`.
pub fn embedding_degree_ok(q_mod_r: u128, r: u128) -> bool {
    embedding_degree_exceeds(q_mod_r, r, EMBEDDING_MIN)
}

/// Stein's binary gcd, by shifts and subtractions: no u128 division, which
/// has no instruction and would cost a library call per Euclidean step.
/// gcd(a, 0) = gcd(0, a) = a.
fn gcd(mut a: u128, mut b: u128) -> u128 {
    if a == 0 || b == 0 {
        return a | b;
    }
    let shift = (a | b).trailing_zeros();
    a >>= a.trailing_zeros();
    loop {
        b >>= b.trailing_zeros();
        if a > b {
            (a, b) = (b, a);
        }
        b -= a;
        if b == 0 {
            return a << shift;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    /// The smallest k with q^k = 1 mod r, by brute force.
    fn order(q: u128, r: u128) -> u32 {
        let (mut x, mut k) = (q % r, 1);
        while x != 1 {
            x = x * (q % r) % r;
            k += 1;
        }
        k
    }

    #[test]
    fn embedding_degree_of_a_nonunit_is_unbounded() {
        // 3 is a nonunit mod 9 without being 0 mod 9; 0 mod 3 and 1 are the ends
        for (q, r) in [(0, 3), (3, 9), (6, 9), (0, 2), (5, 1)] {
            assert!(embedding_degree_exceeds(q, r, 100), "q {q} r {r}");
            assert!(embedding_degree_ok(q, r), "q {q} r {r}");
        }
        assert!(!embedding_degree_exceeds(4, 9, 100)); // 4^3 = 64 = 1 mod 9
        assert!(!embedding_degree_exceeds(1, 2, 1));
    }

    #[test]
    fn embedding_degree_bsgs_matches_brute_force() {
        // small primes, so the brute-force order is quick
        for r in [3u128, 5, 7, 101, 257, 65537, 1_000_003, 999_999_937] {
            for q in [0u128, 2, 3, 10, 12345, 1 << 40] {
                if q % r == 0 {
                    // no power of a nonunit is 1
                    assert!(embedding_degree_exceeds(q % r, r, 100), "q {q} r {r}");
                    continue;
                }
                let k = order(q, r);
                for bound in [1, k - 1, k, k + 1, 1000, EMBEDDING_MIN] {
                    if bound == 0 {
                        continue;
                    }
                    assert_eq!(
                        embedding_degree_exceeds(q % r, r, bound),
                        k > bound,
                        "q {q} r {r} k {k} bound {bound}"
                    );
                }
            }
        }
    }

    fn gcd_euclid(mut a: u128, mut b: u128) -> u128 {
        while b != 0 {
            (a, b) = (b, a % b);
        }
        a
    }

    proptest! {
        #[test]
        fn gcd_matches_euclid(a in any::<u128>(), b in any::<u128>(), k in 0u32..100) {
            prop_assert_eq!(gcd(a, b), gcd_euclid(a, b));
            // a shared power of two, zeros among them
            let (a2, b2) = (a >> k << k, b >> k << k);
            prop_assert_eq!(gcd(a2, b2), gcd_euclid(a2, b2));
            prop_assert_eq!(gcd(a, 0), a);
            prop_assert_eq!(gcd(0, b), b);
        }
    }

    fn is_prime_naive(n: u128) -> bool {
        n >= 2
            && (2..)
                .take_while(|d| d * d <= n)
                .all(|d| !n.is_multiple_of(d))
    }

    proptest! {
        #[test]
        fn mulmod_matches_u128(a in any::<u64>(), b in any::<u64>(), m in 1u128..u64::MAX as u128) {
            prop_assert_eq!(mulmod(a as u128, b as u128, m), (a as u128 * b as u128) % m);
        }

        #[test]
        fn mulmod_large_is_commutative_and_distributive(
            a in any::<u128>(), b in any::<u128>(), c in any::<u128>(), m in (1u128 << 126)..
        ) {
            prop_assert_eq!(mulmod(a, b, m), mulmod(b, a, m));
            let addm = |x: u128, y: u128| (x % m).checked_add(y % m).map_or_else(
                || (x % m).wrapping_add(y % m).wrapping_sub(m), |s| s % m);
            prop_assert_eq!(mulmod(a, addm(b, c), m), addm(mulmod(a, b, m), mulmod(a, c, m)));
        }

        #[test]
        fn mont_mul_matches_mulmod(a in any::<u128>(), b in any::<u128>(), m in any::<u128>()) {
            let (m, top) = (m | 1, u128::MAX - (m & 0xfe));
            for m in [m, top] {
                let mont = Mont::new(m);
                // leaving the form is a product with 1
                let ab = mont.mul(mont.mul(mont.to(a), mont.to(b)), 1);
                prop_assert_eq!(ab, mulmod(a, b, m));
            }
        }

        #[test]
        fn embedding_degree_matches_mulmod(q in any::<u128>(), r in any::<u128>(), k in 1u128..100) {
            // an r of small order k, and a random one, at a bound the plain loop affords
            let bound = 100;
            let slow = |q: u128, r: u128| {
                let mut qk = 1;
                !(0..bound).any(|_| { qk = mulmod(qk, q, r); qk == 1 })
            };
            let small = (1u128 << k) - 1;
            for (q, r) in [(q % r.max(1), r.max(1)), (2, small), (2, small | 1 << 127)] {
                prop_assert_eq!(embedding_degree_exceeds(q, r, bound), slow(q, r));
            }
        }

        #[test]
        fn is_prime_matches_trial_division(n in 0u128..2_000_000) {
            prop_assert_eq!(is_prime(n), is_prime_naive(n));
        }

        #[test]
        fn semiprimes_of_large_primes_are_composite(i in 0usize..4, j in 0usize..4) {
            const PS: [u128; 4] = [18446744073709551557, 2305843009213693951, 2147483647, 4294967291];
            prop_assert!(is_prime(PS[i]) && !is_prime(PS[i] * PS[j]));
        }
    }

    #[test]
    fn known_primes() {
        assert!(is_prime((1u128 << 89) - 1));
        assert!(!is_prime((1u128 << 67) - 1)); // Cole's 193707721 * 761838257287
        // Carmichael numbers
        for n in [
            561u128,
            41041,
            825265,
            321197185,
            5394826801,
            232250619601,
            9746347772161,
        ] {
            assert!(!is_prime(n));
        }
    }
}
