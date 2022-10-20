use crate::{ConstZero, Number, Point2, Vector2};

/// Concentric mapping from the unit square to the unit circle
pub fn concentric_sample_disk<T: Number>(point: Point2<T>) -> Point2<T> {
    // map uniform number to [-1, 1]^2
    let offset = point * T::TWO - Vector2::new(T::ONE, T::ONE);

    if offset.x == T::ZERO && offset.y == T::ZERO {
        return Point2::ZERO;
    }

    let (theta, r) = if offset.x.abs() > offset.y.abs() {
        (T::PI / T::cast(4) * (offset.y / offset.x), offset.x)
    } else {
        (
            T::PI / T::TWO - T::PI / T::cast(4) * (offset.y / offset.x),
            offset.y,
        )
    };

    Point2::new(theta.cos(), theta.sin()) * r
}
