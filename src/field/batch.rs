//! Batch inversion (Montgomery's trick) for any of the fields.
//!
//! The usual prefix-chain implementation is limited by multiplication latency.
//! A product tree exposes independent multiplications at each level.

use crate::field::gf2_127::Gf;
use std::ops::Mul;

pub trait Invert: Copy + Mul<Output = Self> {
    const ONE: Self;
    fn inv(self) -> Self;
}

impl Invert for Gf {
    const ONE: Self = Gf::ONE;
    fn inv(self) -> Self {
        self.invert()
    }
}

/// Invert nonzero elements in place: 3(n - 1) M + 1 I for n > 0.
///
/// Montgomery's trick as a product tree rather than a prefix chain. Each
/// level's multiplications are independent, so the core overlaps as many as
/// its multiplier allows, with no lane count to tune per target; only the
/// top few levels (about log2(latency / throughput)) run at latency.
pub fn invert<F: Invert>(v: &mut [F]) {
    if v.len() <= 1 {
        v.iter_mut().for_each(|x| *x = x.inv());
        return;
    }
    // Levels 1.. of the tree, concatenated in fewer than n + log2(n) elements;
    // level 0 is v. levels[k][i] = levels[k-1][2i] levels[k-1][2i+1], with
    // an odd tail carried up unchanged.
    let mut tree = Vec::with_capacity(v.len());
    up(v, &mut tree);
    let mut levels = vec![(0, tree.len())];
    while let Some(&(s, l)) = levels.last().filter(|&&(_, l)| l > 1) {
        for i in 0..l / 2 {
            let x = tree[s + 2 * i] * tree[s + 2 * i + 1];
            tree.push(x);
        }
        if l % 2 == 1 {
            tree.push(tree[s + l - 1]);
        }
        levels.push((s + l, tree.len() - (s + l)));
    }
    // Invert the root, then walk down in place: 1/a = (1/ab) b, 1/b = (1/ab) a.
    let root = tree.len() - 1;
    tree[root] = tree[root].inv();
    for w in levels.windows(2).rev() {
        let ((s, l), (ps, pl)) = (w[0], w[1]);
        let (lower, upper) = tree.split_at_mut(ps);
        down(&mut lower[s..s + l], &upper[..pl]);
    }
    down(v, &tree[..levels[0].1]);
}

fn up<F: Invert>(lvl: &[F], out: &mut Vec<F>) {
    let (pairs, tail) = lvl.as_chunks::<2>();
    out.extend(
        pairs
            .iter()
            .map(|&[a, b]| a * b)
            .chain(tail.first().copied()),
    );
}

fn down<F: Invert>(lvl: &mut [F], parents: &[F]) {
    let (pairs, tail) = lvl.as_chunks_mut::<2>();
    for (c, &pi) in pairs.iter_mut().zip(parents) {
        let [a, b] = *c;
        *c = [pi * b, pi * a];
    }
    if let [t] = tail {
        *t = parents[parents.len() - 1];
    }
}
