use geometry::{Number, Point2, Vector2};

use crate::filters::Filter;

/// Mitchell-Netravali filter (similar to gaussian, but becomes slightly negative in places
/// reducing the blurriness caused by the gaussian filter)
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Mitchell<T: Number> {
    radius: Vector2<T>,
    inv_radius: Vector2<T>,
    b: T,
    c: T,
}

impl<T: Number> Mitchell<T> {
    /// Create a new mitchell filter with a given radius and parameters b and c,
    /// where b and c are recommended to follow `b + 2c == 1`, although any values
    /// can be used.
    pub fn new(radius: Vector2<T>, b: T, c: T) -> Mitchell<T> {
        Mitchell {
            radius,
            inv_radius: Vector2::new(T::ONE / radius.x, T::ONE / radius.y),
            b,
            c,
        }
    }

    /// 1D mitchell filter implementation
    fn mitchell(&self, x: T) -> T {
        let x = (T::TWO * x).abs();
        let b = self.b;
        let c = self.c;
        let t = T::cast;

        if x > T::ONE {
            ((-b * t(6) * c) * x * x * x
                + (t(6) * b + t(30) * c) * x * x
                + (-t(12) * b - t(48) * c) * x
                + (t(8) * b + t(24) * c))
                * (t(1) / t(6))
        } else {
            ((t(12) - t(9) * b - t(6) * c) * x * x * x
                + (t(-18) + t(12) * b + t(6) * c) * x * x
                + (t(6) - t(2) * b))
                * (t(1) / t(6))
        }
    }
}

impl<T: Number> Filter<T> for Mitchell<T> {
    fn eval(&self, point: Point2<T>) -> T {
        self.mitchell(point.x * self.inv_radius.x) * self.mitchell(point.y * self.inv_radius.y)
    }

    fn radius(&self) -> Vector2<T> {
        self.radius
    }
}
