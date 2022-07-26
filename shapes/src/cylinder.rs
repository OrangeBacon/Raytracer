use std::ops::{Deref, DerefMut};

use geometry::{
    gamma, Bounds3, ConstZero, Number, Point2, Point3, Ray, Transform, Vector2, Vector3,
};

use crate::{PartialDerivatives, Shape, ShapeData, SurfaceInteraction};

/// A cylinder centred on the z axis
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default, Hash)]
pub struct Cylinder<T: Number> {
    data: ShapeData<T>,
    radius: T,
    z_min: T,
    z_max: T,
    phi_max: T,
}

impl<T: Number> Cylinder<T> {
    /// Create a new cylinder instance
    pub fn new(
        object_to_world: Transform<T>,
        world_to_object: Transform<T>,
        reverse_orientation: bool,
        radius: T,
        z: Vector2<T>,
        phi_max: T,
    ) -> Self {
        Self {
            data: ShapeData::new(object_to_world, world_to_object, reverse_orientation),
            radius,
            z_min: z.min_component(),
            z_max: z.max_component(),
            phi_max: phi_max.clamp(T::ZERO, T::cast(360)).to_radians(),
        }
    }
}

impl<T: Number> Shape<T> for Cylinder<T> {
    fn data(&self) -> &ShapeData<T> {
        &self.data
    }

    fn data_mut(&mut self) -> &mut ShapeData<T> {
        &mut self.data
    }

    fn object_bound(&self) -> Bounds3<T> {
        Bounds3::new(
            Point3::new(-self.radius, -self.radius, self.z_min),
            Point3::new(self.radius, self.radius, self.z_max),
        )
    }

    fn intersect(&self, ray: Ray<T>, _test_alpha: bool) -> Option<(T, SurfaceInteraction<T>)> {
        // ignore z-axis for the quadric hit test
        let ray = {
            let mut ray = ray;
            ray.origin.z = T::ZERO;
            ray.direction.z = T::ZERO;
            ray
        };
        let (mut t_shape_hit, t1) = self.quadric_coefficients(self.radius, ray)?;

        let mut p_hit = ray.at(t_shape_hit.value());
        let hit_rad = (p_hit.x * p_hit.x + p_hit.y * p_hit.y).sqrt();
        p_hit.x *= self.radius / hit_rad;
        p_hit.y *= self.radius / hit_rad;
        let mut phi = p_hit.y.atan2(p_hit.x);
        if phi < T::ZERO {
            phi += T::TWO * T::PI;
        }

        if p_hit.z < self.z_min || p_hit.z > self.z_max || phi > self.phi_max {
            if t_shape_hit == t1 {
                return None;
            }
            t_shape_hit = t1;
            if t1.upper_bound() > ray.t_max {
                return None;
            }
            p_hit = ray.at(t_shape_hit.value());
            let hit_rad = (p_hit.x * p_hit.x + p_hit.y * p_hit.y).sqrt();
            p_hit.x *= self.radius / hit_rad;
            p_hit.y *= self.radius / hit_rad;
            phi = p_hit.y.atan2(p_hit.x);
            if phi < T::ZERO {
                phi += T::TWO * T::PI;
            }
            if p_hit.z < self.z_min || p_hit.z > self.z_max || phi > self.phi_max {
                return None;
            }
        }

        let u = phi / self.phi_max;
        let v = (p_hit.z - self.z_min) / (self.z_max - self.z_min);
        let dpdu = Vector3::new(-self.phi_max * p_hit.y, self.phi_max * p_hit.x, T::ZERO);
        let dpdv = Vector3::new(T::ZERO, T::ZERO, self.z_max - self.z_min);
        let d2pdu2 = Vector3::new(p_hit.x, p_hit.y, T::ZERO) * self.phi_max * -self.phi_max;

        let (dndu, dndv) = self.derive_normal(dpdu, dpdv, d2pdu2, Vector3::ZERO, Vector3::ZERO);

        let p_error = Vector3::new(p_hit.x, p_hit.y, T::ZERO).abs() * gamma::<T>(3);

        let intersection = SurfaceInteraction::new(
            p_hit,
            p_error,
            Point2::new(u, v),
            -ray.direction,
            PartialDerivatives {
                dpdu,
                dpdv,
                dndu,
                dndv,
            },
            ray.time,
        )
        .with_shape(*self)
            * self.object_to_world;

        Some((t_shape_hit.value(), intersection))
    }

    fn does_intersect(&self, ray: Ray<T>, _test_alpha: bool) -> bool {
        // ignore z-axis for the quadric hit test
        let ray = {
            let mut ray = ray;
            ray.origin.z = T::ZERO;
            ray.direction.z = T::ZERO;
            ray
        };
        let (mut t_shape_hit, t1) = match self.quadric_coefficients(self.radius, ray) {
            Some(a) => a,
            None => return false,
        };

        let mut p_hit = ray.at(t_shape_hit.value());
        let mut phi = p_hit.y.atan2(p_hit.x);
        if phi < T::ZERO {
            phi += T::TWO * T::PI;
        }

        if p_hit.z < self.z_min || p_hit.z > self.z_max || phi > self.phi_max {
            if t_shape_hit == t1 {
                return false;
            }
            t_shape_hit = t1;
            if t1.upper_bound() > ray.t_max {
                return false;
            }
            p_hit = ray.at(t_shape_hit.value());
            phi = p_hit.y.atan2(p_hit.x);
            if phi < T::ZERO {
                phi += T::TWO * T::PI;
            }
            if p_hit.z < self.z_min || p_hit.z > self.z_max || phi > self.phi_max {
                return false;
            }
        }

        true
    }

    fn area(&self) -> T {
        (self.z_max - self.z_min) * self.radius * self.phi_max
    }
}

impl<T: Number> Deref for Cylinder<T> {
    type Target = ShapeData<T>;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl<T: Number> DerefMut for Cylinder<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}
