use std::{
    f32::consts::PI,
    ops::{Deref, DerefMut},
};

use geometry::{
    Bounds3f, ConstZero, EFloat, Float, Normal3f, PartialDerivatives, Point2f, Point3f, Ray,
    SurfaceInteractable, SurfaceInteraction, Transform, Vector2f, Vector3f,
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
        ray: Ray<(), Float>,
        _test_alpha: bool,
    ) -> Option<(Float, SurfaceInteraction<&dyn SurfaceInteractable>)> {
        // transform ray to object space
        let (ray, (o_err, d_err)) = self.world_to_object.apply_err(ray);

        let ox = EFloat::new_with_err(ray.origin.x, o_err.x);
        let oy = EFloat::new_with_err(ray.origin.y, o_err.y);
        let oz = EFloat::new_with_err(ray.origin.z, o_err.z);
        let dx = EFloat::new_with_err(ray.direction.x, d_err.x);
        let dy = EFloat::new_with_err(ray.direction.y, d_err.y);
        let dz = EFloat::new_with_err(ray.direction.z, d_err.z);

        let a = dx * dx + dy * dy + dz * dz;
        let b = 2.0 * (dx * ox + dy * oy + dz * oz);
        let c = ox * ox + oy * oy + oz * oz - EFloat::new(self.radius) * EFloat::new(self.radius);

        let (t0, t1) = EFloat::quadratic(a, b, c)?;

        if t0.upper_bound() > ray.t_max || t1.lower_bound() <= 0.0 {
            return None;
        }

        let mut t_shape_hit = if t0.lower_bound() <= 0.0 {
            if t1.upper_bound() > ray.t_max {
                return None;
            }
            t1
        } else {
            t0
        };

        // compute sphere hit position
        let mut p_hit = ray.at(t_shape_hit.value());
        if p_hit.x == 0.0 && p_hit.y == 0.0 {
            p_hit.x = 1e-5 * self.radius;
        }

        let mut phi = p_hit.y.atan2(p_hit.x);
        if phi < 0.0 {
            phi += 2.0 * PI
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
            if p_hit.x == 0.0 && p_hit.y == 0.0 {
                p_hit.x = 1e-5 * self.radius;
            }

            let mut phi = p_hit.y.atan2(p_hit.x);
            if phi < 0.0 {
                phi += 2.0 * PI
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
        let theta = (p_hit.z / self.radius).clamp(-1.0, 1.0).acos();
        let v = (theta - self.theta_min) / (self.theta_max - self.theta_min);

        let z_radius = (p_hit.x * p_hit.x + p_hit.y * p_hit.y).sqrt();
        let inv_z_radius = 1.0 / z_radius;
        let cos_phi = p_hit.x * inv_z_radius;
        let sin_phi = p_hit.y * inv_z_radius;
        let dpdu = Vector3f::new(-self.phi_max * p_hit.y, self.phi_max * p_hit.x, 0.0);
        let dpdv = (self.theta_max - self.theta_min)
            * Vector3f::new(
                p_hit.z * cos_phi,
                p_hit.z * sin_phi,
                -self.radius * theta.sin(),
            );

        // Normal vector derivatives

        let d2pdu2 = -self.phi_max * self.phi_max * Vector3f::new(p_hit.x, p_hit.y, 0.0);
        let d2pduv = (self.theta_max - self.theta_min)
            * p_hit.z
            * self.phi_max
            * Vector3f::new(-sin_phi, cos_phi, 0.0);
        let d2pdv2 =
            -(self.theta_max - self.theta_min) * (self.theta_max * self.theta_min) * p_hit.to_vec();

        let e1 = dpdu.dot(dpdu);
        let f1 = dpdu.dot(dpdv);
        let g1 = dpdv.dot(dpdv);
        let normal = dpdu.cross(dpdv).normalise();
        let e2 = normal.dot(d2pdu2);
        let f2 = normal.dot(d2pduv);
        let g2 = normal.dot(d2pdv2);

        let inv_egf2 = 1.0 / (e1 * g1 - f1 * f1);
        let dndu = Normal3f::from_vector(
            (f2 * f1 - e2 * g1) * inv_egf2 * dpdu + (e2 * f1 - f2 * e1) * inv_egf2 * dpdv,
        );
        let dndv = Normal3f::from_vector(
            (g2 * f1 - f2 * g1) * inv_egf2 * dpdu + (f2 * f1 - g2 * e1) * inv_egf2 * dpdv,
        );

        let intersection = SurfaceInteraction::new(
            p_hit,
            Vector3f::ZERO,
            Point2f::new(u, v),
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

    fn does_intersect(&self, ray: Ray<(), Float>, _test_alpha: bool) -> bool {
        // transform ray to object space
        let (ray, (o_err, d_err)) = self.world_to_object.apply_err(ray);

        let ox = EFloat::new_with_err(ray.origin.x, o_err.x);
        let oy = EFloat::new_with_err(ray.origin.y, o_err.y);
        let oz = EFloat::new_with_err(ray.origin.z, o_err.z);
        let dx = EFloat::new_with_err(ray.direction.x, d_err.x);
        let dy = EFloat::new_with_err(ray.direction.y, d_err.y);
        let dz = EFloat::new_with_err(ray.direction.z, d_err.z);

        let a = dx * dx + dy * dy + dz * dz;
        let b = 2.0 * (dx * ox + dy * oy + dz * oz);
        let c = ox * ox + oy * oy + oz * oz - EFloat::new(self.radius) * EFloat::new(self.radius);

        let (t0, t1) = match EFloat::quadratic(a, b, c) {
            Some(v) => v,
            None => return false,
        };

        if t0.upper_bound() > ray.t_max || t1.lower_bound() <= 0.0 {
            return false;
        }

        let mut t_shape_hit = if t0.lower_bound() <= 0.0 {
            if t1.upper_bound() > ray.t_max {
                return false;
            }
            t1
        } else {
            t0
        };

        // compute sphere hit position
        let mut p_hit = ray.at(t_shape_hit.value());
        if p_hit.x == 0.0 && p_hit.y == 0.0 {
            p_hit.x = 1e-5 * self.radius;
        }

        let mut phi = p_hit.y.atan2(p_hit.x);
        if phi < 0.0 {
            phi += 2.0 * PI
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
            if p_hit.x == 0.0 && p_hit.y == 0.0 {
                p_hit.x = 1e-5 * self.radius;
            }

            let mut phi = p_hit.y.atan2(p_hit.x);
            if phi < 0.0 {
                phi += 2.0 * PI
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
