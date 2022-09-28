use std::ops::{Deref, DerefMut};

use geometry::{Bounds3, ConstZero, Normal3, Number, Point2, Point3, Ray, Transform, Vector3};

use crate::{PartialDerivatives, Shape, ShapeData, SurfaceInteraction};

/// A disk with inner and outer radii centred on the z axis
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default, Hash)]
pub struct Disk<T: Number> {
    data: ShapeData<T>,
    height: T,
    radius: T,
    inner_radius: T,
    phi_max: T,
}

impl<T: Number> Disk<T> {
    /// Create a new disk instance
    pub fn new(
        object_to_world: Transform<T>,
        world_to_object: Transform<T>,
        reverse_orientation: bool,
        height: T,
        radius: T,
        inner_radius: T,
        phi_max: T,
    ) -> Self {
        Self {
            data: ShapeData::new(object_to_world, world_to_object, reverse_orientation),
            height,
            radius,
            inner_radius,
            phi_max: phi_max.clamp(T::ZERO, T::cast(360)).to_radians(),
        }
    }
}

impl<T: Number> Shape<T> for Disk<T> {
    fn data(&self) -> &ShapeData<T> {
        &self.data
    }

    fn data_mut(&mut self) -> &mut ShapeData<T> {
        &mut self.data
    }

    fn object_bound(&self) -> Bounds3<T> {
        Bounds3::new(
            Point3::new(-self.radius, -self.radius, self.height),
            Point3::new(self.radius, self.radius, self.height),
        )
    }

    fn intersect(&self, ray: Ray<T>, _test_alpha: bool) -> Option<(T, SurfaceInteraction<T>)> {
        let (ray, _) = self.world_to_object.apply_err(ray);

        let t_shape_hit = (self.height - ray.origin.z) / ray.direction.z;
        if t_shape_hit <= T::ZERO || t_shape_hit >= ray.t_max {
            return None;
        }

        // no intersections parallel to the disk
        if ray.direction.z == T::ZERO {
            return None;
        }

        let mut p_hit = ray.at(t_shape_hit);
        let dist2 = p_hit.x * p_hit.x + p_hit.y * p_hit.y;
        if dist2 > self.radius * self.radius || dist2 < self.inner_radius * self.inner_radius {
            return None;
        }

        let mut phi = p_hit.y.atan2(p_hit.x);
        if phi < T::ZERO {
            phi += T::TWO * T::PI;
        }
        if phi > self.phi_max {
            return None;
        }

        let u = phi / self.phi_max;
        let r_hit = dist2.sqrt();
        let one_minus_v = (r_hit - self.inner_radius) / (self.radius - self.inner_radius);
        let v = T::ONE - one_minus_v;
        let dpdu = Vector3::new(-self.phi_max * p_hit.y, self.phi_max * p_hit.x, T::ZERO);
        let dpdv =
            Vector3::new(p_hit.x, p_hit.y, T::ZERO) * (self.inner_radius - self.radius) / r_hit;

        p_hit.z = self.height;

        let intersection = SurfaceInteraction::new(
            p_hit,
            Vector3::ZERO,
            Point2::new(u, v),
            -ray.direction,
            PartialDerivatives {
                dpdu,
                dpdv,
                dndu: Normal3::ZERO,
                dndv: Normal3::ZERO,
            },
            ray.time,
        )
        .with_shape(*self)
            * self.object_to_world;

        Some((t_shape_hit, intersection))
    }

    fn does_intersect(&self, ray: Ray<T>, _test_alpha: bool) -> bool {
        let (ray, _) = self.world_to_object.apply_err(ray);

        let t_shape_hit = (self.height - ray.origin.z) / ray.direction.z;
        if t_shape_hit <= T::ZERO || t_shape_hit >= ray.t_max {
            return false;
        }

        // no intersections parallel to the disk
        if ray.direction.z == T::ZERO {
            return false;
        }

        let p_hit = ray.at(t_shape_hit);
        let dist2 = p_hit.x * p_hit.x + p_hit.y * p_hit.y;
        if dist2 > self.radius * self.radius || dist2 < self.inner_radius * self.inner_radius {
            return false;
        }

        let mut phi = p_hit.y.atan2(p_hit.x);
        if phi < T::ZERO {
            phi += T::TWO * T::PI;
        }
        if phi > self.phi_max {
            return false;
        }

        true
    }

    fn area(&self) -> T {
        self.phi_max
            * (T::HALF)
            * (self.radius * self.radius - self.inner_radius * self.inner_radius)
    }
}

impl<T: Number> Deref for Disk<T> {
    type Target = ShapeData<T>;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl<T: Number> DerefMut for Disk<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}
