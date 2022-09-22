mod cone;
mod curve;
mod cylinder;
mod disk;
mod hyperboloid;
mod interaction;
mod paraboloid;
mod quadric;
mod sphere;
mod triangle;
mod triangle_mesh;

use std::ops::{Deref, DerefMut};

pub use cone::Cone;
pub use curve::Curve;
pub use curve::CurveCommon;
pub use curve::CurveType;
pub use cylinder::Cylinder;
pub use disk::Disk;
pub use hyperboloid::Hyperboloid;
pub use interaction::*;
pub use paraboloid::Paraboloid;
pub use sphere::Sphere;
pub use triangle::Triangle;
pub use triangle_mesh::TriangleMesh;

use geometry::{Bounds3, Number, Ray, Transform};

/// Data that should be stored by all shapes.
/// TODO: Transform cache + references, don't copy 4 mat4x4 with each shape
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default, Hash)]
pub struct ShapeData<T: Number> {
    pub object_to_world: Transform<T>,
    pub world_to_object: Transform<T>,
    pub reverse_orientation: bool,
    pub transform_swaps_handedness: bool,
}

/// General interface for all shapes
pub trait Shape<T: Number>: Deref<Target = ShapeData<T>> + DerefMut + SurfaceInteractable {
    /// Get the bounding box of the shape in object space
    fn object_bound(&self) -> Bounds3<T>;

    /// Get the bounding box of the shape in world space
    fn world_bound(&self) -> Bounds3<T> {
        self.object_bound() * self.object_to_world
    }

    /// Find the intersection of self and a ray
    fn intersect(
        &self,
        ray: Ray<(), T>,
        test_alpha: bool,
    ) -> Option<(T, SurfaceInteraction<&dyn SurfaceInteractable, (), T>)>;

    /// Find whether self and a ray intersect
    fn does_intersect(&self, ray: Ray<(), T>, test_alpha: bool) -> bool {
        self.intersect(ray, test_alpha).is_some()
    }

    /// Get the surface area of a shape
    fn area(&self) -> T;
}

impl<T: Number> ShapeData<T> {
    /// Create a new shape data from transform matrices
    fn new(
        object_to_world: Transform<T>,
        world_to_object: Transform<T>,
        reverse_orientation: bool,
    ) -> Self {
        Self {
            object_to_world,
            world_to_object,
            reverse_orientation,
            transform_swaps_handedness: object_to_world.swaps_handedness(),
        }
    }
}
