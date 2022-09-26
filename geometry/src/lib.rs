mod animated_transform;
mod bounds2;
mod bounds3;
mod float;
mod interval;
mod matrix4x4;
mod normal3;
mod number;
mod point2;
mod point3;
mod quaternion;
mod ray;
mod rng;
mod transform;
mod vector2;
mod vector3;

#[cfg(test)]
mod tests;

#[cfg(feature = "double")]
pub type Float = f64;

#[cfg(not(feature = "double"))]
pub type Float = f32;

pub use number::ConstZero;
pub use number::Integer;
pub use number::Number;

pub use float::*;

pub use vector2::Vector2;
pub use vector2::Vector2f;
pub use vector2::Vector2i;

pub use vector3::Vector3;
pub use vector3::Vector3f;
pub use vector3::Vector3i;

pub use point2::Point2;
pub use point2::Point2f;
pub use point2::Point2i;

pub use point3::Point3;
pub use point3::Point3f;
pub use point3::Point3i;

pub use normal3::Normal3;
pub use normal3::Normal3f;
pub use normal3::Normal3i;

pub use ray::Ray;
pub use ray::RayDifferential;

pub use bounds3::Bounds3;
pub use bounds3::Bounds3f;
pub use bounds3::Bounds3i;

pub use bounds2::Bounds2;
pub use bounds2::Bounds2f;
pub use bounds2::Bounds2i;

pub use matrix4x4::Matrix4x4;

pub use transform::Applicable;
pub use transform::Transform;

pub use quaternion::Quaternion;

pub use animated_transform::AnimatedTransform;

pub use interval::Interval;

pub use rng::Rng;

/// Linearly interpolate between two floats
pub fn lerp<T: Number>(t: T, a: T, b: T) -> T {
    (T::ONE - t) * a + t * b
}

/// Solve a quadratic equation ax^2 + bx + c = 0
/// If no real solutions are found, returns None.  If both solutions are the
/// same, that solution will be both return values.
pub fn quadratic<T: Number>(a: T, b: T, c: T) -> Option<(T, T)> {
    let discriminant = b * b - T::cast(4) * a * c;
    if discriminant < T::ZERO {
        return None;
    }
    let root = discriminant.sqrt();

    let q = if b < T::ZERO {
        -T::HALF * (b - root)
    } else {
        -T::HALF * (b + root)
    };

    let t0 = q / a;
    let t1 = c / q;

    Some((t0.min(t1) as _, t0.max(t1) as _))
}

/// Solve equation Ax = B where A is a 2x2 matrix and B is a 2 element vector
/// returns x as (x1, x2) if solution found, otherwise None
pub fn solve_2x2_system<T: Number>(a: [[T; 2]; 2], b: [T; 2]) -> Option<(T, T)> {
    let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    if det.abs() < T::cast(1e-10) {
        return None;
    }

    let x0 = (a[1][1] * b[0] - a[0][1] * b[1]) / det;
    let x1 = (a[0][0] * b[1] - a[1][0] * b[0]) / det;

    if x0.is_nan() || x1.is_nan() {
        return None;
    }

    Some((x0, x1))
}

pub fn offset_ray_origin<T: Number>(
    point: Point3<T>,
    error: Vector3<T>,
    normal: Normal3<T>,
    w: Vector3<T>,
) -> Point3<T> {
    let d = normal.to_vector().abs().dot(error);
    let mut offset = normal.to_vector() * d;
    if normal.to_vector().dot(w) < T::ZERO {
        offset = -offset;
    }
    let mut po = point + offset;

    for i in 0..3 {
        if offset[i] > T::ZERO {
            po[i] = next_float_up(po[i]);
        } else if offset[i] < T::ZERO {
            po[i] = next_float_down(po[i]);
        }
    }

    po
}

/// Integer base 2 logarithm of a number
pub fn log2u32(num: u32) -> u32 {
    31 - num.leading_zeros()
}
pub fn log2usize(num: usize) -> usize {
    std::mem::size_of::<usize>() * 8 - 1 - num.leading_zeros() as usize
}

/// Round a number up to the nearest power of 2.  Returns the input if it is
/// a power of two.
pub fn round_up_pow_2(mut x: usize) -> usize {
    x -= 1;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    x |= x >> 32;
    x + 1
}

/// Is the input number a power of two
pub fn is_pow_2(x: usize) -> bool {
    x != 0 && (x & (x - 1)) == 0
}
