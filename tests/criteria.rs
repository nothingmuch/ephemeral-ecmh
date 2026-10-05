//! Direct selection needs an order oracle, not group encodings or certificates.
use ephemeral_ecmh::curvegen::criteria::{self, AdmissibleR, Count, Criteria, Filter, OrderError};

const R: u128 = (1 << 20) + 7;
struct OrderOnly;
impl Criteria for OrderOnly {
    type Curve = u32;
    const TAG: &'static [u8] = b"order-only";
    const R: AdmissibleR = AdmissibleR::new(R + 5, 6, 1);
    fn candidate(_: &[u8; 32], j: u32) -> Option<u32> {
        (j != 0).then_some(j)
    }
}
struct Counts(Vec<u32>);
impl Count<u32> for Counts {
    fn order(&mut self, c: &u32) -> u128 {
        self.0.push(*c);
        if *c == 3 { R } else { R + 1 }
    }
}
struct SkipOne;
impl Filter<u32> for SkipOne {
    fn rejects(&mut self, c: &u32) -> bool {
        *c == 1
    }
}
#[test]
fn selection_needs_neither_arithmetic_nor_a_witness() {
    let mut count = Counts(vec![]);
    let found = criteria::find::<OrderOnly>(&[0; 32], &mut count, &mut SkipOne);
    assert_eq!((found.index, found.r, found.curve), (3, R, 3));
    assert_eq!(count.0, [2, 3]);
}
#[test]
fn order_rejections_are_values_not_certificate_panics() {
    assert_eq!(OrderOnly::check_order(R), Ok(R));
    assert_eq!(OrderOnly::check_order(0), Err(OrderError::Hasse));
    assert_eq!(OrderOnly::check_order(R + 1), Err(OrderError::RNotPrime));
    let anomalous = AdmissibleR::new(R, 6, 1);
    assert_eq!(anomalous.check_order(R), Err(OrderError::Anomalous));
    let mov = AdmissibleR::new(R + 1, 6, 1);
    assert_eq!(mov.check_order(R), Err(OrderError::Mov));
    let h = AdmissibleR::new(2 * R, 6, 2);
    assert_eq!(h.check_order(2 * R + 1), Err(OrderError::Cofactor));
}
