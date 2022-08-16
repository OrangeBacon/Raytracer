mod animated_transform;
mod bounds2;
mod bounds3;
mod interval;
mod matrix4x4;
mod normal3;
mod number;
mod point2;
mod point3;
mod quaternion;
mod ray;
mod transform;
mod vector2;
mod vector3;

mod interaction;
#[cfg(test)]
mod tests;
pub mod float;

#[cfg(feature = "double")]
mod types {
    pub type Float = f64;
    pub type FloatBits = u64;
}

#[cfg(not(feature = "double"))]
mod types {
    pub type Float = f32;
    pub type FloatBits = u32;
}

pub use types::Float;
pub use types::FloatBits;

pub use number::ConstZero;
pub use number::Integer;
pub use number::Number;

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

pub use transform::Transform;

pub use quaternion::Quaternion;

pub use animated_transform::AnimatedTransform;

pub use interval::Interval;

pub use interaction::Interaction;
pub use interaction::SurfaceInteractable;
pub use interaction::SurfaceInteraction;

/// Linearly interpolate between two floats
pub fn lerp(t: Float, a: Float, b: Float) -> Float {
    (1.0 - t) * a + t * b
}

/// Solve a quadratic equation ax^2 + bx + c = 0
/// If no real solutions are found, returns None.  If both solutions are the
/// same, that solution will be both return values.
pub fn quadratic(a: Float, b: Float, c: Float) -> Option<(Float, Float)> {
    let a = a as f64;
    let b = b as f64;
    let c = c as f64;
    let discriminant = b * b - 4.0 * a * c;
    if discriminant < 0.0 {
        return None;
    }
    let root = discriminant.sqrt();

    let q = if b < 0.0 {
        -0.5 * (b - root)
    } else {
        -0.5 * (b + root)
    };

    let t0 = q / a;
    let t1 = c / q;

    Some((t0.min(t1) as _, t0.max(t1) as _))
}

/// Solve equation Ax = B where A is a 2x2 matrix and B is a 2 element vector
/// returns x as (x1, x2) if solution found, otherwise None
pub fn solve_2x2_system(a: [[Float; 2]; 2], b: [Float; 2]) -> Option<(Float, Float)> {
    let det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    if det.abs() < 1e-10 {
        return None;
    }

    let x0 = (a[1][1] * b[0] - a[0][1] * b[1]) / det;
    let x1 = (a[0][0] * b[1] - a[1][0] * b[0]) / det;

    if x0.is_nan() || x1.is_nan() {
        return None;
    }

    Some((x0, x1))
}
