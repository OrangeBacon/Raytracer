use std::ops::{Deref, DerefMut};

use geometry::{
    float::EFloat, Bounds3f, ConstZero, Float, Point3f, SurfaceInteractable, SurfaceInteraction,
    Transform, Vector2f,
};

use crate::{Shape, ShapeData};

/// 3D sphere, represented using spherical coordinates.
/// x^2 + y^2 + z^2 = r^2;
/// x = r sin theta cos phi;
/// y = r sin theta sin phi;
/// z = r cos theta;
/// only renders the portion of the sphere between z_min and z_max
pub struct Sphere {
    data: ShapeData,
    radius: Float,
    z_min: Float,
    z_max: Float,
    theta_min: Float,
    theta_max: Float,
    phi_max: Float,
}

impl Sphere {
    /// Create a new sphere instance
    pub fn new(
        object_to_world: Transform,
        world_to_object: Transform,
        reverse_orientation: bool,
        radius: Float,
        z: Vector2f,
        phi_max: Float,
    ) -> Self {
        let z_min = z.min_component().clamp(-radius, radius);
        let z_max = z.max_component().clamp(-radius, radius);
        Self {
            data: ShapeData {
                object_to_world,
                world_to_object,
                reverse_orientation,
                transform_swaps_handedness: object_to_world.swaps_handedness(),
            },
            radius,
            z_min,
            z_max,
            theta_min: (z_min / radius).clamp(-1.0, 1.0).acos(),
            theta_max: (z_max / radius).clamp(-1.0, 1.0).acos(),
            phi_max: phi_max.clamp(0.0, 360.0).to_radians(),
        }
    }
}

impl Shape for Sphere {
    fn object_bound(&self) -> Bounds3f {
        Bounds3f::new(
            Point3f::new(-self.radius, -self.radius, self.z_min),
            Point3f::new(self.radius, self.radius, self.z_max),
        )
    }

    fn intersect(
        &self,
        ray: geometry::Ray,
        test_alpha: bool,
    ) -> Option<(Float, SurfaceInteraction<&dyn SurfaceInteractable>)> {
        todo!()
    }

    fn area(&self) -> Float {
        self.phi_max * self.radius * (self.z_max - self.z_min)
    }
}

impl SurfaceInteractable for Sphere {
    fn reverses_orientation(&self) -> bool {
        self.reverse_orientation ^ self.transform_swaps_handedness
    }
}

impl Deref for Sphere {
    type Target = ShapeData;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl DerefMut for Sphere {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}
