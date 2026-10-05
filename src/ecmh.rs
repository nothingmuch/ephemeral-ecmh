//! ECMH over a curve implementing the required group capabilities:
//! H(multiset) = sum of H(item).
//!
//! `Ecmh` is the streaming form (one hashed point into an accumulator per
//! insert/remove, as a RIBLT cell does). `digest_batch` is the bulk form:
//! hash and sum through the group's batch operations. Each asks of the
//! curve only the `group` capabilities it uses.

use crate::group::{Encode, Group, HashToCurve, Negate, SumBatch};
use crate::hash::Salted;

pub const TAG_ITEM: &[u8] = b"ephemeral-ecmh/item";

/// Streaming multiset hash under one salt on one curve.
pub struct Ecmh<G: Group> {
    pub group: G,
    pub salt: Salted,
    pub acc: G::Point,
}

impl<G: HashToCurve + Negate + Encode> Ecmh<G> {
    pub fn new(group: G, salt: &[u8; 32]) -> Self {
        Self {
            group,
            salt: Salted::new(TAG_ITEM, salt),
            acc: group.identity(),
        }
    }

    pub fn insert(&mut self, item: &[u8]) {
        let p = self.group.hash(&self.salt, item);
        self.acc = self.group.add_affine(&self.acc, &p);
    }

    pub fn remove(&mut self, item: &[u8]) {
        let p = self.group.neg(&self.group.hash(&self.salt, item));
        self.acc = self.group.add_affine(&self.acc, &p);
    }

    pub fn digest(&self) -> G::Encoding {
        self.group.encode_point(&self.acc)
    }
}

/// Bulk multiset hash using the group's batch operations.
pub fn digest_batch<G: HashToCurve + SumBatch + Encode>(
    group: G,
    salt: &[u8; 32],
    items: &[&[u8]],
) -> G::Encoding {
    let h = Salted::new(TAG_ITEM, salt);
    group.encode(&group.sum_batch(&group.hash_batch(&h, items)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::curve::binary::lambda;
    use crate::curve::binary127;
    use crate::curve::edwards127;
    use crate::curve::weier127;
    use crate::group::{Accumulate, Decode};
    use proptest::prelude::*;

    type Items = Vec<Vec<u8>>;

    fn items() -> impl Strategy<Value = Items> {
        prop::collection::vec(prop::collection::vec(any::<u8>(), 0..12), 0..12)
    }

    fn streaming<G: HashToCurve + Negate + Encode>(
        g: G,
        salt: &[u8; 32],
        xs: &Items,
    ) -> G::Encoding {
        let mut e = Ecmh::new(g, salt);
        xs.iter().for_each(|x| e.insert(x));
        e.digest()
    }

    fn batch<G: HashToCurve + SumBatch + Encode>(g: G, salt: &[u8; 32], xs: &Items) -> G::Encoding {
        let refs: Vec<&[u8]> = xs.iter().map(|x| x.as_slice()).collect();
        digest_batch(g, salt, &refs)
    }

    /// Batch encode/decode agree with one-at-a-time on running sums (O included)
    /// and on every single-bit flip of O's and the first sum's encodings;
    /// whatever decodes re-encodes to the same bytes.
    fn check_codec<G: HashToCurve + Negate + Decode>(
        g: G,
        salt: [u8; 32],
        xs: &Items,
    ) -> Result<(), TestCaseError> {
        let h = Salted::new(TAG_ITEM, &salt);
        let pts: Vec<G::Point> = std::iter::once(g.identity())
            .chain(xs.iter().scan(g.identity(), |p, x| {
                *p = g.add_affine(p, &g.hash(&h, x));
                Some(*p)
            }))
            .collect();
        let es = g.encode_batch(&pts);
        let one: Vec<_> = pts.iter().map(|p| g.encode(&g.to_affine(p))).collect();
        prop_assert_eq!(&es, &one);
        let direct: Vec<_> = pts.iter().map(|p| g.encode_point(p)).collect();
        prop_assert_eq!(&direct, &one);
        let mut all = es.clone();
        for e in es.iter().take(2) {
            for bit in 0..8 * e.as_ref().len() {
                let mut f = *e;
                f.as_mut()[bit / 8] ^= 1 << (bit % 8);
                all.push(f);
            }
        }
        prop_assert!(g.decode_batch(&[]).is_empty());
        let re = |a: Option<G::Affine>| a.map(|a| g.encode(&a));
        let batch: Vec<_> = g.decode_batch(&all).into_iter().map(re).collect();
        let one: Vec<_> = all.iter().map(|e| re(g.decode(e))).collect();
        prop_assert_eq!(&batch, &one);
        for (e, d) in all.iter().zip(&batch) {
            prop_assert!(d.is_none_or(|d| d == *e));
        }
        for (e, d) in es.iter().zip(&batch) {
            prop_assert_eq!(*d, Some(*e));
        }
        Ok(())
    }

    fn check_laws<G: HashToCurve + Negate + SumBatch + Decode>(
        g: G,
        salt: [u8; 32],
        xs: Items,
        ys: Items,
        perm: prop::sample::Index,
    ) -> Result<(), TestCaseError> {
        check_codec(g, salt, &xs)?;
        let empty = streaming(g, &salt, &vec![]);
        prop_assert_eq!(streaming(g, &salt, &xs), batch(g, &salt, &xs));
        let mut shuffled = xs.clone();
        if !shuffled.is_empty() {
            let k = perm.index(shuffled.len());
            shuffled.rotate_left(k);
            shuffled.reverse();
        }
        prop_assert_eq!(batch(g, &salt, &shuffled), batch(g, &salt, &xs));
        let mut e = Ecmh::new(g, &salt);
        xs.iter().chain(&ys).for_each(|x| e.insert(x));
        xs.iter().for_each(|x| e.remove(x));
        prop_assert_eq!(e.digest(), streaming(g, &salt, &ys));
        prop_assert_eq!(g.is_identity(&e.acc), ys.is_empty());
        ys.iter().for_each(|x| e.remove(x));
        prop_assert_eq!(e.digest(), empty);
        prop_assert!(g.is_identity(&e.acc));
        let h = Salted::new(TAG_ITEM, &salt);
        let pts: Vec<_> = xs.iter().map(|x| g.hash(&h, x)).collect();
        let one: Vec<_> = pts.iter().map(|p| g.prepare(p)).collect();
        // equals_addend on unnormalized accumulators: accepts the item, rejects
        // its negation and O, matches affine equality after another item is
        // added, and accepts again once it is removed.
        for (a, b) in one.iter().zip(one.iter().skip(1)) {
            let p = g.add(&g.identity(), a);
            prop_assert!(g.equals_addend(&p, a));
            prop_assert!(!g.equals_addend(&p, &g.neg_addend(a)));
            prop_assert!(!g.equals_addend(&g.identity(), a));
            let pb = g.add(&p, b);
            prop_assert_eq!(
                g.equals_addend(&pb, a),
                g.encode(&g.to_affine(&pb)) == g.encode(&g.to_affine(&p))
            );
            prop_assert!(g.equals_addend(&g.add(&pb, &g.neg_addend(b)), a));
        }
        for adds in [one, g.prepare_batch(&pts)] {
            let acc = adds.iter().fold(g.identity(), |acc, a| g.add(&acc, a));
            prop_assert_eq!(g.encode(&g.to_affine(&acc)), streaming(g, &salt, &xs));
            let acc = adds
                .iter()
                .fold(acc, |acc, a| g.add(&acc, &g.neg_addend(a)));
            prop_assert_eq!(g.encode(&g.to_affine(&acc)), empty);
        }
        // The digest depends on the salt (a collision is possible but negligible).
        if !xs.is_empty() {
            let mut other = salt;
            other[0] ^= 1;
            prop_assert_ne!(batch(g, &salt, &xs), batch(g, &other, &xs));
        }
        Ok(())
    }

    /// The λ group is `binary127`'s under another accumulator: every digest
    /// and every running sum of prepared addends, signed, must agree.
    fn check_lambda_matches_binary127(
        c: binary127::Curve,
        salt: [u8; 32],
        xs: Items,
        signs: Vec<bool>,
    ) -> Result<(), TestCaseError> {
        let l = lambda::Curve(c);
        prop_assert_eq!(streaming(l, &salt, &xs), streaming(c, &salt, &xs));
        prop_assert_eq!(batch(l, &salt, &xs), batch(c, &salt, &xs));
        let h = Salted::new(TAG_ITEM, &salt);
        let pts: Vec<_> = xs.iter().map(|x| c.hash(&h, x)).collect();
        let (la, ca) = (l.prepare_batch(&pts), c.prepare_batch(&pts));
        let (mut lp, mut cp) = (l.identity(), c.identity());
        // Revisit items in reverse order with either sign to exercise doubling
        // and cancellation to O as well as additions of distinct points.
        let order = (0..xs.len()).chain((0..xs.len()).rev());
        for (i, &sub) in order.zip(signs.iter().cycle()) {
            let (a, b) = if sub {
                (l.neg_addend(&la[i]), c.neg_addend(&ca[i]))
            } else {
                (la[i], ca[i])
            };
            lp = l.add(&lp, &a);
            cp = c.add(&cp, &b);
            prop_assert_eq!(l.to_affine(&lp), c.to_affine(&cp));
            prop_assert_eq!(l.is_identity(&lp), c.is_identity(&cp));
        }
        Ok(())
    }

    /// `check_lambda_matches_binary127` for any two accumulators over the
    /// same curve, points and encodings (a Weierstrass family's projective
    /// and Jacobian accumulators): digests, and signed running sums of
    /// prepared addends, agree.
    fn check_same_digests<G, H>(
        g: G,
        h: H,
        salt: [u8; 32],
        xs: Items,
        signs: Vec<bool>,
    ) -> Result<(), TestCaseError>
    where
        G: HashToCurve + Negate + SumBatch + Encode,
        H: Group<Affine = G::Affine>
            + HashToCurve
            + Negate
            + SumBatch
            + Encode<Encoding = G::Encoding>,
    {
        prop_assert_eq!(streaming(g, &salt, &xs), streaming(h, &salt, &xs));
        prop_assert_eq!(batch(g, &salt, &xs), batch(h, &salt, &xs));
        let s = Salted::new(TAG_ITEM, &salt);
        let pts: Vec<_> = xs.iter().map(|x| g.hash(&s, x)).collect();
        let (ga, ha) = (g.prepare_batch(&pts), h.prepare_batch(&pts));
        let (mut gp, mut hp) = (g.identity(), h.identity());
        // Revisit items in reverse order with either sign to exercise doubling
        // and cancellation to O as well as additions of distinct points.
        let order = (0..xs.len()).chain((0..xs.len()).rev());
        for (i, &sub) in order.zip(signs.iter().cycle()) {
            let (a, b) = if sub {
                (g.neg_addend(&ga[i]), h.neg_addend(&ha[i]))
            } else {
                (ga[i], ha[i])
            };
            gp = g.add(&gp, &a);
            hp = h.add(&hp, &b);
            prop_assert_eq!(g.encode(&g.to_affine(&gp)), h.encode(&h.to_affine(&hp)));
            prop_assert_eq!(g.is_identity(&gp), h.is_identity(&hp));
        }
        Ok(())
    }

    proptest! {
        #[test]
        fn binary127_ecmh(c in binary127::tests::curve(), salt in any::<[u8; 32]>(), xs in items(), ys in items(), perm in any::<prop::sample::Index>()) {
            check_laws(c, salt, xs, ys, perm)?;
        }

        #[test]
        fn binary127_lambda_ecmh(c in binary127::tests::curve(), salt in any::<[u8; 32]>(), xs in items(), ys in items(), perm in any::<prop::sample::Index>()) {
            check_laws(lambda::Curve(c), salt, xs, ys, perm)?;
        }

        #[test]
        fn binary127_lambda_matches_binary127(c in binary127::tests::curve(), salt in any::<[u8; 32]>(), xs in items(), signs in prop::collection::vec(any::<bool>(), 1..8)) {
            check_lambda_matches_binary127(c, salt, xs, signs)?;
        }

        #[test]
        fn edwards127_ecmh(c in edwards127::tests::curve(), salt in any::<[u8; 32]>(), xs in items(), ys in items(), perm in any::<prop::sample::Index>()) {
            check_laws(c, salt, xs, ys, perm)?;
        }

        #[test]
        fn weier127_ecmh(c in weier127::tests::odd_curve(), salt in any::<[u8; 32]>(), xs in items(), ys in items(), perm in any::<prop::sample::Index>()) {
            check_laws(c, salt, xs, ys, perm)?;
        }

        #[test]
        fn weier127_jacobian_ecmh(c in weier127::tests::odd_curve(), salt in any::<[u8; 32]>(), xs in items(), ys in items(), perm in any::<prop::sample::Index>()) {
            check_laws(crate::curve::weier::jacobian::Curve(c), salt, xs, ys, perm)?;
        }

        #[test]
        fn weier127_jacobian_matches_weier127(c in weier127::tests::odd_curve(), salt in any::<[u8; 32]>(), xs in items(), signs in prop::collection::vec(any::<bool>(), 1..8)) {
            check_same_digests(c, crate::curve::weier::jacobian::Curve(c), salt, xs, signs)?;
        }

        #[test]
        fn binary109_ecmh(c in crate::curve::binary109::tests::curve(), salt in any::<[u8; 32]>(), xs in items(), ys in items(), perm in any::<prop::sample::Index>()) {
            check_laws(c, salt, xs, ys, perm)?;
        }

        #[test]
        fn binary109_lambda_ecmh(c in crate::curve::binary109::tests::curve(), salt in any::<[u8; 32]>(), xs in items(), ys in items(), perm in any::<prop::sample::Index>()) {
            check_laws(lambda::Curve(c), salt, xs, ys, perm)?;
        }

        #[test]
        fn binary109_lambda_matches_binary109(c in crate::curve::binary109::tests::curve(), salt in any::<[u8; 32]>(), xs in items(), signs in prop::collection::vec(any::<bool>(), 1..8)) {
            check_same_digests(c, lambda::Curve(c), salt, xs, signs)?;
        }

        #[test]
        fn binary127_w_ecmh(c in binary127::tests::curve(), salt in any::<[u8; 32]>(), xs in items(), ys in items(), perm in any::<prop::sample::Index>()) {
            check_laws(crate::curve::binary::wcodec::Curve::new(c), salt, xs, ys, perm)?;
        }

        #[test]
        fn binary109_w_ecmh(c in crate::curve::binary109::tests::curve(), salt in any::<[u8; 32]>(), xs in items(), ys in items(), perm in any::<prop::sample::Index>()) {
            check_laws(crate::curve::binary::wcodec::Curve::new(c), salt, xs, ys, perm)?;
        }

        #[test]
        fn binary127_unscaled_ecmh(c in binary127::tests::curve(), salt in any::<[u8; 32]>(), xs in items(), ys in items(), perm in any::<prop::sample::Index>()) {
            check_laws(crate::curve::binary::unscaled::Curve::new(c), salt, xs, ys, perm)?;
        }

        #[test]
        fn binary109_unscaled_ecmh(c in crate::curve::binary109::tests::curve(), salt in any::<[u8; 32]>(), xs in items(), ys in items(), perm in any::<prop::sample::Index>()) {
            check_laws(crate::curve::binary::unscaled::Curve::new(c), salt, xs, ys, perm)?;
        }
    }
}
