use std::ops::{Deref, DerefMut};

use geometry::{Bounds3, EFloat, Number, Point2, Point3, Ray, Transform, Vector2, Vector3};

use crate::{PartialDerivatives, Shape, ShapeData, SurfaceInteractable, SurfaceInteraction};

/// A paraboloid centred on the z axis
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default, Hash)]
pub struct Paraboloid<T: Number> {
    data: ShapeData<T>,
    radius: T,
    z_min: T,
    z_max: T,
    phi_max: T,
}

impl<T: Number> Paraboloid<T> {
    /// Create a new paraboloid instance
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

impl<T: Number> Shape<T> for Paraboloid<T> {
    fn object_bound(&self) -> geometry::Bounds3<T> {
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
        let (ray, (o_err, d_err)) = self.world_to_object.apply_err(ray);

        let ox = EFloat::new_with_err(ray.origin.x, o_err.x);
        let oy = EFloat::new_with_err(ray.origin.y, o_err.y);
        let oz = EFloat::new_with_err(ray.origin.z, o_err.z);
        let dx = EFloat::new_with_err(ray.direction.x, d_err.x);
        let dy = EFloat::new_with_err(ray.direction.y, d_err.y);
        let dz = EFloat::new_with_err(ray.direction.z, d_err.z);

        let k = EFloat::new(self.z_max) / (EFloat::new(self.radius) * EFloat::new(self.radius));

        let a = k * (dx * dx + dy * dy);
        let b = (dx * ox + dy * oy) * T::TWO * k - dz;
        let c = k * (ox * ox + oy * oy) - oz;

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

        let mut p_hit = ray.at(t_shape_hit.value());
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
        let dpdv = Vector3::new(
            p_hit.x / (T::TWO * p_hit.z),
            p_hit.y / (T::TWO * p_hit.z),
            T::ONE,
        ) * (self.z_max - self.z_min);

        let d2pdu2 = Vector3::new(p_hit.x, p_hit.y, T::ZERO) * -self.phi_max * self.phi_max;
        let d2pduv = Vector3::new(
            p_hit.x / (T::TWO * p_hit.z),
            p_hit.y / (T::TWO * p_hit.z),
            T::ZERO,
        ) * self.phi_max
            * (self.z_max - self.z_min);
        let d2pdv2 = Vector3::new(
            p_hit.x / (T::cast(4) * p_hit.z * p_hit.z),
            p_hit.y / (T::cast(4) * p_hit.z * p_hit.z),
            T::ZERO,
        ) * -(self.z_max - self.z_min)
            * (self.z_max - self.z_min);

        let (dndu, dndv) = self.derive_normal(dpdu, dpdv, d2pdu2, d2pduv, d2pdv2);

        let px = ox + t_shape_hit * dx;
        let py = oy + t_shape_hit * dy;
        let pz = oz + t_shape_hit * dz;
        let p_err = Vector3::new(
            px.absolute_error(),
            py.absolute_error(),
            pz.absolute_error(),
        );

        let intersection = SurfaceInteraction::new(
            p_hit,
            p_err,
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
        let (ray, (o_err, d_err)) = self.world_to_object.apply_err(ray);

        let ox = EFloat::new_with_err(ray.origin.x, o_err.x);
        let oy = EFloat::new_with_err(ray.origin.y, o_err.y);
        let oz = EFloat::new_with_err(ray.origin.z, o_err.z);
        let dx = EFloat::new_with_err(ray.direction.x, d_err.x);
        let dy = EFloat::new_with_err(ray.direction.y, d_err.y);
        let dz = EFloat::new_with_err(ray.direction.z, d_err.z);

        let k = EFloat::new(self.z_max) / (EFloat::new(self.radius) * EFloat::new(self.radius));

        let a = k * (dx * dx + dy * dy);
        let b = (dx * ox + dy * oy) * T::TWO * k - dz;
        let c = k * (ox * ox + oy * oy) - oz;

        let (t0, t1) = match EFloat::quadratic(a, b, c) {
            Some(a) => a,
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
        let radius_2 = self.radius * self.radius;
        let k = T::cast(4) * self.z_max / radius_2;
        (radius_2 * radius_2 * self.phi_max / (T::cast(12) * self.z_max * self.z_max))
            * ((k * self.z_max + T::ONE).pow(T::cast(1.5))
                - (k * self.z_min + T::ONE).pow(T::cast(1.5)))
    }
}

impl<T: Number> SurfaceInteractable for Paraboloid<T> {
    fn reverses_orientation(&self) -> bool {
        self.reverse_orientation ^ self.transform_swaps_handedness
    }
}

impl<T: Number> Deref for Paraboloid<T> {
    type Target = ShapeData<T>;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl<T: Number> DerefMut for Paraboloid<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}
