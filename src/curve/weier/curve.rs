//! The curve for one namespace's b, and the odd-order curves that are groups.

use super::Affine;
use crate::field::OddField;

/// y^2 = x^3 - 3x + b; `new` admits only nonsingular curves with b != 0,
/// which the methods assume of a literal. Its formulas are those of any
/// such curve: certificate verification computes on candidates of every
/// order, and RCB's are exceptional where P - Q has order 2 (`select`'s
/// `Arithmetic` states what a verdict then means). The group is
/// `OddCurve`.
#[derive(Clone, Copy, Debug)]
pub struct Curve<F> {
    pub b: F,
}

impl<F: OddField> Curve<F> {
    /// `None` for singular curves (4a^3 + 27b^2 = 0, i.e. b = ±2, and every
    /// b in characteristic 3) or b = 0, which would put (0, 0) on the curve
    /// and collide with `Affine::IDENTITY`.
    pub fn new(b: F) -> Option<Self> {
        let two = F::ONE + F::ONE;
        let three = two + F::ONE;
        (!three.is_zero() && b != F::ZERO && b != two && b != -two).then_some(Self { b })
    }

    pub fn rhs(&self, x: F) -> F {
        x * (x.square() - (F::ONE + F::ONE + F::ONE)) + self.b
    }

    /// Whether x is an x-coordinate.
    pub fn x_on_curve(&self, x: F) -> bool {
        self.rhs(x).is_square()
    }

    pub fn is_on_curve(&self, p: &Affine<F>) -> bool {
        p.is_identity() || p.y.square() == self.rhs(p.x)
    }
}

/// A curve of odd order, the `group` traits' Weierstrass group: a root x0
/// of x^3 - 3x + b is the point (x0, 0) of order 2, so the cubic has no
/// root exactly when #E(F) is odd. There RCB's formulas are complete; the
/// chord law is complete on every curve.
#[derive(Clone, Copy, Debug)]
pub struct OddCurve<F>(Curve<F>);

impl<F: OddField> OddCurve<F> {
    /// `None` unless x^3 - 3x + b is irreducible. That also refuses what
    /// `Curve::new` does: b = 0 (the root 0), b = ±2 (a double root) and
    /// characteristic 3 (x^3 + b has the root -b^(1/3)). The test is
    /// gcd(x^q - x, x^3 - 3x + b): about log2 q squarings modulo the cubic,
    /// once per curve.
    pub fn new(c: Curve<F>) -> Option<Self> {
        let three = F::ONE + F::ONE + F::ONE;
        cubic_irreducible(-three, c.b).then_some(Self(c))
    }

    pub fn curve(&self) -> Curve<F> {
        self.0
    }
}

/// A cubic is irreducible iff gcd(x^q - x, x³ + ax + b) = 1.
// Fixed-degree arithmetic keeps this check generic over the curve's field.
pub(super) fn cubic_irreducible<F: OddField>(a: F, b: F) -> bool {
    let multiply = |left: [F; 3], right: [F; 3]| {
        let mut product = [F::ZERO; 5];
        for (i, x) in left.into_iter().enumerate() {
            for (j, y) in right.into_iter().enumerate() {
                product[i + j] = product[i + j] + x * y;
            }
        }
        for k in (3..5).rev() {
            product[k - 3] = product[k - 3] - product[k] * b;
            product[k - 2] = product[k - 2] - product[k] * a;
        }
        [product[0], product[1], product[2]]
    };
    let mut power = [F::ZERO, F::ONE, F::ZERO];
    let mut result = [F::ONE, F::ZERO, F::ZERO];
    let mut q = F::ORDER;
    while q != 0 {
        if q & 1 != 0 {
            result = multiply(result, power);
        }
        q >>= 1;
        if q != 0 {
            power = multiply(power, power);
        }
    }
    result[1] = result[1] - F::ONE;
    let mut left = vec![b, a, F::ZERO, F::ONE];
    let mut right = result.to_vec();
    trim(&mut right);
    while let Some(leading) = right.last() {
        let inverse = leading.inv();
        while left.len() >= right.len() {
            let factor = *left.last().unwrap() * inverse;
            let offset = left.len() - right.len();
            for (x, &y) in left[offset..].iter_mut().zip(&right) {
                *x = *x - factor * y;
            }
            trim(&mut left);
        }
        (left, right) = (right, left);
    }
    left.len() == 1
}

fn trim<F: OddField>(coefficients: &mut Vec<F>) {
    while coefficients.last().is_some_and(|c| c.is_zero()) {
        coefficients.pop();
    }
}
