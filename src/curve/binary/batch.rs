//! Batched affine arithmetic: one shared inversion for many additions.

use super::{Affine, Model};
use crate::field::Field;
use crate::field::batch::Invert;
use crate::field::batch::invert as batch_invert;

#[derive(Clone, Copy)]
enum Slot {
    Done,
    Chord,
    Double,
}

/// out\[i\] = a\[i\] + b\[i\], with one shared inversion.
///
/// Chord: 5M + 1S per slot. P = Q, P = -Q and O inputs branch, which is
/// acceptable on public data; those cases are routine (peeling computes
/// P - P).
pub fn add_batch<M: Model>(a: &[Affine<M>], b: &[Affine<M>], out: &mut [Affine<M>]) {
    let n = out.len();
    assert!(a.len() == n && b.len() == n);
    let mut den = vec![M::F::ONE; n];
    let mut slot = vec![Slot::Done; n];
    for i in 0..n {
        let (p, q) = (&a[i], &b[i]);
        if p.is_identity() {
            out[i] = *q;
        } else if q.is_identity() {
            out[i] = *p;
        } else if p.x.equals(q.x) {
            if p.y.equals(q.y) {
                den[i] = p.x;
                slot[i] = Slot::Double;
            } else {
                out[i] = Affine::IDENTITY;
            }
        } else {
            den[i] = p.x + q.x;
            slot[i] = Slot::Chord;
        }
    }
    batch_invert(&mut den);
    for i in 0..n {
        let (p, q) = (&a[i], &b[i]);
        match slot[i] {
            Slot::Done => {}
            Slot::Chord => {
                let l = (p.y + q.y) * den[i];
                let x = l.square() + l + p.x + q.x + M::A;
                out[i] = Affine {
                    x,
                    y: l * (p.x + x) + x + p.y,
                };
            }
            Slot::Double => {
                let l = p.x + p.y * den[i];
                let x = l.square() + l + M::A;
                out[i] = Affine {
                    x,
                    y: p.x.square() + l * x + x,
                };
            }
        }
    }
}

/// Sum of many points by pairwise tree levels, one inversion per level.
pub fn sum_batch<M: Model>(points: &[Affine<M>]) -> Affine<M> {
    let mut level = points.to_vec();
    while level.len() > 1 {
        let half = level.len() / 2;
        let mut next = vec![Affine::IDENTITY; half];
        add_batch(&level[..half], &level[half..2 * half], &mut next);
        if level.len() % 2 == 1 {
            next.push(level[level.len() - 1]);
        }
        level = next;
    }
    level.first().copied().unwrap_or(Affine::IDENTITY)
}
