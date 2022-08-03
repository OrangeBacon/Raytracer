mod bounds3;
mod normal3;
mod number;
mod point2;
mod point3;
mod ray;
mod vector2;
mod vector3;
mod bounds2;

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