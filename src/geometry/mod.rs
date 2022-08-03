mod number;
mod point2;
mod point3;
mod vector2;
mod vector3;

#[cfg(feature = "double")]
pub type Float = f64;

#[cfg(not(feature = "double"))]
pub type Float = f32;

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
