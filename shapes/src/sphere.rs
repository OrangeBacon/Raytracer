use std::ops::{Deref, DerefMut};

use geometry::{
    Bounds3, ConstZero, EFloat, Normal3, Number, PartialDerivatives, Point2, Point3, Ray,
    SurfaceInteractable, SurfaceInteraction, Transform, Vector2, Vector3,
};

use crate::{Shape, ShapeData};

/// 3D sphere, represented using spherical coordinates.
/// x^2 + y^2 + z^2 = r^2;
/// x = r sin theta cos phi;
/// y = r sin theta sin phi;
/// z = r cos theta;
/// only renders the portion of the sphere between z_min and z_max
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
            data: ShapeData {
                object_to_world,
                world_to_object,
                reverse_orientation,
                transform_swaps_handedness: object_to_world.swaps_handedness(),
            },
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
        // transform ray to object space
        let (ray, (o_err, d_err)) = self.world_to_object.apply_err(ray);

        let ox = EFloat::new_with_err(ray.origin.x, o_err.x);
        let oy = EFloat::new_with_err(ray.origin.y, o_err.y);
        let oz = EFloat::new_with_err(ray.origin.z, o_err.z);
        let dx = EFloat::new_with_err(ray.direction.x, d_err.x);
        let dy = EFloat::new_with_err(ray.direction.y, d_err.y);
        let dz = EFloat::new_with_err(ray.direction.z, d_err.z);

        let a = dx * dx + dy * dy + dz * dz;
        let b = (dx * ox + dy * oy + dz * oz) * T::TWO;
        let c = ox * ox + oy * oy + oz * oz - EFloat::new(self.radius) * EFloat::new(self.radius);

        let (t0, t1) = EFloat::quadratic(a, b, c)?;

        if t0.upper_bound() > ray.t_max || t1.lower_bound() <= T::ZERO {
            return None;
        }

        let mut t_shape_hit = if t0.lower_bound() <= T::ZERO {
            if t1.upper_bound() > ray.t_max {
                return None;
            }
            t1
        } else {
            t0
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

        let e1 = dpdu.dot(dpdu);
        let f1 = dpdu.dot(dpdv);
        let g1 = dpdv.dot(dpdv);
        let normal = dpdu.cross(dpdv).normalise();
        let e2 = normal.dot(d2pdu2);
        let f2 = normal.dot(d2pduv);
        let g2 = normal.dot(d2pdv2);

        let inv_egf2 = T::ONE / (e1 * g1 - f1 * f1);
        let dndu = Normal3::from_vector(
            dpdu * (f2 * f1 - e2 * g1) * inv_egf2 + dpdv * (e2 * f1 - f2 * e1) * inv_egf2,
        );
        let dndv = Normal3::from_vector(
            dpdu * (g2 * f1 - f2 * g1) * inv_egf2 + dpdv * (f2 * f1 - g2 * e1) * inv_egf2,
        );

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
        let (ray, (o_err, d_err)) = self.world_to_object.apply_err(ray);

        let ox = EFloat::new_with_err(ray.origin.x, o_err.x);
        let oy = EFloat::new_with_err(ray.origin.y, o_err.y);
        let oz = EFloat::new_with_err(ray.origin.z, o_err.z);
        let dx = EFloat::new_with_err(ray.direction.x, d_err.x);
        let dy = EFloat::new_with_err(ray.direction.y, d_err.y);
        let dz = EFloat::new_with_err(ray.direction.z, d_err.z);

        let a = dx * dx + dy * dy + dz * dz;
        let b = (dx * ox + dy * oy + dz * oz) * T::TWO;
        let c = ox * ox + oy * oy + oz * oz - EFloat::new(self.radius) * EFloat::new(self.radius);

        let (t0, t1) = match EFloat::quadratic(a, b, c) {
            Some(v) => v,
            None => return false,
        };

        if t0.upper_bound() > ray.t_max || t1.lower_bound() <= T::ZERO {
            return false;
        }

        let mut t_shape_hit = if t0.lower_bound() <= T::ZERO {
            if t1.upper_bound() > ray.t_max {
                return false;
            }
            t1
        } else {
            t0
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
