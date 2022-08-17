use geometry::{Bounds3f, Float, Ray, SurfaceInteractable, SurfaceInteraction};

/// General interface for all shapes
pub trait Shape: SurfaceInteractable {
    /// Get the bounding box of the shape in object space
    fn object_bound(&self) -> Bounds3f;

    /// Get the bounding box of the shape in world space
    fn world_bound(&self) -> Bounds3f;

    /// Find the intersection of self and a ray
    fn intersect(
        &self,
        ray: Ray,
        test_alpha: bool,
    ) -> Option<(Float, SurfaceInteraction<&dyn SurfaceInteractable>)>;

    /// Find whether self and a ray intersect
    fn does_intersect(&self, ray: Ray, test_alpha: bool) -> bool {
        self.intersect(ray, test_alpha).is_some()
    }

    /// Get the surface area of a shape
    fn area(&self) -> Float;
}
