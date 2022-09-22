use geometry::{Number, Point2};

use crate::filters::Filter;

/// Simple filter applying weights equally to all samples
/// Very easy to compute, but will introduce aliasing
pub struct Box {}

impl Box {
    /// Create a new box filter
    pub fn new() -> Box {
        Box {}
    }
}

impl<T: Number> Filter<T> for Box {
    fn eval(&self, _: Point2<T>) -> T {
        T::ONE
    }
}
