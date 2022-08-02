mod number;
pub mod vector2;
pub mod vector3;

#[cfg(feature = "double")]
pub type Float = f64;

#[cfg(not(feature = "double"))]
pub type Float = f32;
