use std::ops::{Deref, DerefMut};

use geometry::{
    Bounds3, ConstZero, Number, PartialDerivatives, Point2, Point3, Ray, SurfaceInteractable,
    SurfaceInteraction, Transform, Vector2, Vector3,
};

use crate::{Shape, ShapeData};

/// 3D sphere, represented using spherical coordinates.
/// x^2 + y^2 + z^2 = r^2;
/// x = r sin theta cos phi;
/// y = r sin theta sin phi;
/// z = r cos theta;
/// only renders the portion of the sphere between z_min and z_max
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Default)]
pub struct Sphere<T: Number> {
    data: ShapeData<T>,
    radius: T,
    z_min: T,
    z_max: T,
    theta_min: T,
    theta_max: T,
    phi_max: T,
}

impl<T: Number> Sphere<T> {
    /// Create a new sphere instance
    pub fn new(
        object_to_world: Transform<T>,
        world_to_object: Transform<T>,
        reverse_orientation: bool,
        radius: T,
        z: Vector2<T>,
        phi_max: T,
    ) -> Self {
        let z_min = z.min_component().clamp(-radius, radius);
        let z_max = z.max_component().clamp(-radius, radius);
        Self {
            data: ShapeData::new(object_to_world, world_to_object, reverse_orientation),
            radius,
            z_min,
            z_max,
            theta_min: (z_min / radius).clamp(-T::ONE, T::ONE).acos(),
            theta_max: (z_max / radius).clamp(-T::ONE, T::ONE).acos(),
            phi_max: phi_max.clamp(T::ZERO, T::cast(360)).to_radians(),
        }
    }
}

impl<T: Number> Shape<T> for Sphere<T> {
    fn object_bound(&self) -> Bounds3<T> {
        Bounds3::new(
            Point3::new(-self.radius, -self.radius, self.z_min),
            Point3::new(self.radius, self.radius, self.z_max),
        )
    }

    fn intersect(
        &self,
        ray: Ray<(), T>,
        _test_alpha: bool,
    ) -> Option<(T, SurfaceInteraction<&dyn SurfaceInteractable, (), T>)> {
        let (mut t_shape_hit, t1) = self.quadric_coefficients(self.radius, ray)?;

        // compute sphere hit position
        let mut p_hit = ray.at(t_shape_hit.value());
        if p_hit.x == T::ZERO && p_hit.y == T::ZERO {
            p_hit.x = T::cast(1e-5) * self.radius;
        }

        let mut phi = p_hit.y.atan2(p_hit.x);
        if phi < T::ZERO {
            phi += T::TWO * T::PI;
        }

        // test sphere intersection against clipping parameters
        if (self.z_min > -self.radius && p_hit.z < self.z_min)
            || (self.z_max < self.radius && p_hit.z > self.z_max)
            || phi > self.phi_max
        {
            if t_shape_hit == t1 || t1.upper_bound() > ray.t_max {
                return None;
            }
            t_shape_hit = t1;

            p_hit = ray.at(t_shape_hit.value());
            if p_hit.x == T::ZERO && p_hit.y == T::ZERO {
                p_hit.x = T::cast(1e-5) * self.radius;
            }

            let mut phi = p_hit.y.atan2(p_hit.x);
            if phi < T::ZERO {
                phi += T::TWO * T::PI;
            }

            if (self.z_min > -self.radius && p_hit.z < self.z_min)
                || (self.z_max < self.radius && p_hit.z > self.z_max)
                || phi > self.phi_max
            {
                return None;
            }
        }

        // Parametric representation of sphere hit
        let u = phi / self.phi_max;
        let theta = (p_hit.z / self.radius).clamp(-T::ONE, T::ONE).acos();
        let v = (theta - self.theta_min) / (self.theta_max - self.theta_min);

        let z_radius = (p_hit.x * p_hit.x + p_hit.y * p_hit.y).sqrt();
        let inv_z_radius = T::ONE / z_radius;
        let cos_phi = p_hit.x * inv_z_radius;
        let sin_phi = p_hit.y * inv_z_radius;
        let dpdu = Vector3::new(-self.phi_max * p_hit.y, self.phi_max * p_hit.x, T::ZERO);
        let dpdv = Vector3::new(
            p_hit.z * cos_phi,
            p_hit.z * sin_phi,
            -self.radius * theta.sin(),
        ) * (self.theta_max - self.theta_min);

        // Normal vector derivatives

        let d2pdu2 = Vector3::new(p_hit.x, p_hit.y, T::ZERO) * -self.phi_max * self.phi_max;
        let d2pduv = Vector3::new(-sin_phi, cos_phi, T::ZERO)
            * (self.theta_max - self.theta_min)
            * p_hit.z
            * self.phi_max;
        let d2pdv2 =
            p_hit.to_vec() * -(self.theta_max - self.theta_min) * (self.theta_max * self.theta_min);

        let (dndu, dndv) = self.derive_normal(dpdu, dpdv, d2pdu2, d2pduv, d2pdv2);

        let intersection = SurfaceInteraction::new(
            p_hit,
            Vector3::ZERO,
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
        .with_shape(self as &dyn SurfaceInteractable)
            * self.object_to_world;

        Some((t_shape_hit.value(), intersection))
    }

    fn does_intersect(&self, ray: Ray<(), T>, _test_alpha: bool) -> bool {
        // transform ray to object space
        let (mut t_shape_hit, t1) = match self.quadric_coefficients(self.radius, ray) {
            Some(a) => a,
            None => return false,
        };

        // compute sphere hit position
        let mut p_hit = ray.at(t_shape_hit.value());
        if p_hit.x == T::ZERO && p_hit.y == T::ZERO {
            p_hit.x = T::cast(1e-5) * self.radius;
        }

        let mut phi = p_hit.y.atan2(p_hit.x);
        if phi < T::ZERO {
            phi += T::TWO * T::PI;
        }

        // test sphere intersection against clipping parameters
        if (self.z_min > -self.radius && p_hit.z < self.z_min)
            || (self.z_max < self.radius && p_hit.z > self.z_max)
            || phi > self.phi_max
        {
            if t_shape_hit == t1 || t1.upper_bound() > ray.t_max {
                return false;
            }
            t_shape_hit = t1;

            p_hit = ray.at(t_shape_hit.value());
            if p_hit.x == T::ZERO && p_hit.y == T::ZERO {
                p_hit.x = T::cast(1e-5) * self.radius;
            }

            let mut phi = p_hit.y.atan2(p_hit.x);
            if phi < T::ZERO {
                phi += T::TWO * T::PI;
            }

            if (self.z_min > -self.radius && p_hit.z < self.z_min)
                || (self.z_max < self.radius && p_hit.z > self.z_max)
                || phi > self.phi_max
            {
                return false;
            }
        }

        true
    }

    fn area(&self) -> T {
        self.phi_max * self.radius * (self.z_max - self.z_min)
    }
}

impl<T: Number> SurfaceInteractable for Sphere<T> {
    fn reverses_orientation(&self) -> bool {
        self.reverse_orientation ^ self.transform_swaps_handedness
    }
}

impl<T: Number> Deref for Sphere<T> {
    type Target = ShapeData<T>;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl<T: Number> DerefMut for Sphere<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}
