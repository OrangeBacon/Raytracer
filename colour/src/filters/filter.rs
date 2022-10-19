use geometry::{Number, Point2, Vector2};

/// Interface to a filter using in image reconstruction
pub trait Filter<T: Number> {
    fn eval(&self, point: Point2<T>) -> T;

    fn radius(&self) -> Vector2<T>;
}
