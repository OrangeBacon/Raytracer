use geometry::{Number, Point2, Vector2};

use crate::filters::Filter;

/// Simple filter applying weights equally to all samples
/// Very easy to compute, but will introduce aliasing
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Default)]
pub struct Box<T: Number> {
    radius: Vector2<T>,
}

impl<T: Number> Box<T> {
    /// Create a new box filter
    pub fn new(radius: Vector2<T>) -> Box<T> {
        Box { radius }
    }
}

impl<T: Number> Filter<T> for Box<T> {
    fn eval(&self, _: Point2<T>) -> T {
        T::ONE
    }

    fn radius(&self) -> Vector2<T> {
        self.radius
    }
}
