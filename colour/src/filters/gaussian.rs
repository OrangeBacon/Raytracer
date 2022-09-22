use geometry::{Number, Point2, Vector2};

use crate::filters::Filter;

pub struct Gaussian<T: Number> {
    alpha: T,
    exp_x: T,
    exp_y: T,
}

impl<T: Number> Gaussian<T> {
    /// create a new gaussian filter with a given radius and falloff rate (alpha)
    pub fn new(radius: Vector2<T>, alpha: T) -> Gaussian<T> {
        Gaussian {
            alpha,
            exp_x: (-alpha * radius.x * radius.x).exp(),
            exp_y: (-alpha * radius.y * radius.y).exp(),
        }
    }

    /// Compute the 1D gaussian
    fn gaussian(&self, d: T, exp_v: T) -> T {
        T::ZERO.max((-self.alpha * d * d).exp() - exp_v)
    }
}

impl<T: Number> Filter<T> for Gaussian<T> {
    fn eval(&self, point: Point2<T>) -> T {
        self.gaussian(point.x, self.exp_x) * self.gaussian(point.y, self.exp_y)
    }
}
