use geometry::{Number, Vector2};

use crate::filters::Filter;

pub struct Triangle<T: Number> {
    radius: Vector2<T>,
}

impl<T: Number> Triangle<T> {
    /// Create a new triangle filter with a given radius
    pub fn new(radius: Vector2<T>) -> Triangle<T> {
        Triangle { radius }
    }
}

impl<T: Number> Filter<T> for Triangle<T> {
    fn eval(&self, point: geometry::Point2<T>) -> T {
        T::ZERO.max(self.radius.x - point.x.abs()) * T::ZERO.max(self.radius.y - point.y.abs())
    }
}
