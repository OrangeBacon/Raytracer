use std::fmt::Debug;

use geometry::{Bounds3, Number, Ray};

use crate::SurfaceInteraction;

/// Link between geometry and shading code
pub trait Primitive<T: Number>: Debug {
    /// Bound for the primitive in world space
    fn world_bound(&self) -> Bounds3<T>;

    /// Get an intersection between the primitive and a given ray
    fn intersect(&self, ray: Ray<(), T>) -> Option<SurfaceInteraction<(), T>>;

    /// Does the ray intersect the given primitive, faster version of
    /// `Self::intersect(..).is_some()`
    fn does_intersect(&self, ray: Ray<(), T>) -> bool;

    /// Get the emissive properties of the primitive
    fn area_light(&self);

    /// Get the material instance assigned to the primitive
    fn material(&self);

    /// Create representations of the light scattering properties of the material
    fn scattering_functions(&self);
}
