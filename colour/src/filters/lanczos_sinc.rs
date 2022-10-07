use geometry::{Number, Point2, Vector2};

use crate::filters::Filter;

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LanczosSinc<T: Number> {
    radius: Vector2<T>,
    tau: T,
}

impl<T: Number> LanczosSinc<T> {
    /// Create a new windowed lanczos sinc filter
    pub fn new(radius: Vector2<T>, tau: T) -> LanczosSinc<T> {
        LanczosSinc { radius, tau }
    }

    fn sinc(&self, x: T) -> T {
        let x = x.abs();
        if x < T::cast(1e-5) {
            return T::ONE;
        }

        (T::PI * x).sin() / (T::PI * x)
    }

    fn windowed_sinc(&self, x: T, radius: T) -> T {
        let x = x.abs();
        if x > radius {
            return T::ZERO;
        }

        let lanczos = self.sinc(x / self.tau);
        self.sinc(x) * lanczos
    }
}

impl<T: Number> Filter<T> for LanczosSinc<T> {
    fn eval(&self, point: Point2<T>) -> T {
        self.windowed_sinc(point.x, self.radius.x) * self.windowed_sinc(point.y, self.radius.y)
    }
}
