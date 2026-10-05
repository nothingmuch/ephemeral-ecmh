//! The affine chord law: one inversion per addition, shared across a
//! batch, and tree sums built on it.

use super::{Affine, Curve};
use crate::field::OddField;
use crate::field::batch::invert as batch_invert;

impl<F: OddField> Curve<F> {
    pub fn add(&self, p: &Affine<F>, q: &Affine<F>) -> Affine<F> {
        let mut out = [Affine::IDENTITY];
        self.add_batch(core::slice::from_ref(p), core::slice::from_ref(q), &mut out);
        out[0]
    }

    /// out\[i\] = a\[i\] + b\[i\] with one shared inversion; chord 5M + 1S.
    pub fn add_batch(&self, a: &[Affine<F>], b: &[Affine<F>], out: &mut [Affine<F>]) {
        let n = out.len();
        assert!(a.len() == n && b.len() == n);
        let mut num = vec![F::ZERO; n];
        let mut den = vec![F::ONE; n];
        let mut live = vec![false; n];
        for i in 0..n {
            let (p, q) = (&a[i], &b[i]);
            if p.is_identity() {
                out[i] = *q;
            } else if q.is_identity() {
                out[i] = *p;
            } else if p.x == q.x {
                if p.y != q.y || p.y.is_zero() {
                    out[i] = Affine::IDENTITY;
                } else {
                    let t = p.x.square() - F::ONE;
                    num[i] = t + t + t;
                    den[i] = p.y + p.y;
                    live[i] = true;
                }
            } else {
                num[i] = q.y - p.y;
                den[i] = q.x - p.x;
                live[i] = true;
            }
        }
        batch_invert(&mut den);
        for i in 0..n {
            if live[i] {
                let (p, q) = (&a[i], &b[i]);
                let l = num[i] * den[i];
                let x = l.square() - p.x - q.x;
                out[i] = Affine {
                    x,
                    y: l * (p.x - x) - p.y,
                };
            }
        }
    }

    pub fn sum_batch(&self, points: &[Affine<F>]) -> Affine<F> {
        let mut level = points.to_vec();
        while level.len() > 1 {
            let half = level.len() / 2;
            let mut next = vec![Affine::IDENTITY; half];
            self.add_batch(&level[..half], &level[half..2 * half], &mut next);
            if level.len() % 2 == 1 {
                next.push(level[level.len() - 1]);
            }
            level = next;
        }
        level.first().copied().unwrap_or(Affine::IDENTITY)
    }
}
