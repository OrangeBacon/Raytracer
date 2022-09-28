use std::ops::{Deref, DerefMut};

use geometry::{Bounds3, ConstZero, EFloat, Number, Point2, Point3, Ray, Transform, Vector3};

use crate::{PartialDerivatives, Shape, ShapeData, SurfaceInteraction};

/// A Cone centred on the z axis
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default, Hash)]
pub struct Cone<T: Number> {
    data: ShapeData<T>,
    radius: T,
    height: T,
    phi_max: T,
}

impl<T: Number> Cone<T> {
    ///  Create a new cone instance
    pub fn new(
        object_to_world: Transform<T>,
        world_to_object: Transform<T>,
        reverse_orientation: bool,
        radius: T,
        height: T,
        phi_max: T,
    ) -> Self {
        Self {
            data: ShapeData::new(object_to_world, world_to_object, reverse_orientation),
            radius,
            height,
            phi_max: phi_max.clamp(T::ZERO, T::cast(360)).to_radians(),
        }
    }
}

impl<T: Number> Shape<T> for Cone<T> {
    fn data(&self) -> &ShapeData<T> {
        &self.data
    }

    fn data_mut(&mut self) -> &mut ShapeData<T> {
        &mut self.data
    }

    fn object_bound(&self) -> Bounds3<T> {
        Bounds3::new(
            Point3::new(-self.radius, -self.radius, T::ZERO),
            Point3::new(self.radius, self.radius, self.height),
        )
    }

    fn intersect(&self, ray: Ray<T>, _test_alpha: bool) -> Option<(T, SurfaceInteraction<T>)> {
        let (ray, (o_err, d_err)) = self.world_to_object.apply_err(ray);

        let ox = EFloat::new_with_err(ray.origin.x, o_err.x);
        let oy = EFloat::new_with_err(ray.origin.y, o_err.y);
        let oz = EFloat::new_with_err(ray.origin.z, o_err.z);
        let dx = EFloat::new_with_err(ray.direction.x, d_err.x);
        let dy = EFloat::new_with_err(ray.direction.y, d_err.y);
        let dz = EFloat::new_with_err(ray.direction.z, d_err.z);

        let k = EFloat::new(self.radius) / EFloat::new(self.height);
        let k = k * k;

        let a = dx * dx + dy * dy - k * dz * dz;
        let b = (dx * ox + dy * oy - k * dz * (oz - self.height)) * T::TWO;
        let c = ox * ox + oy * oy - k * (oz - self.height) * (oz - self.height);

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

        if p_hit.z < T::ZERO || p_hit.z > self.height || phi > self.phi_max {
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
            if p_hit.z < T::ZERO || p_hit.z > self.height || phi > self.phi_max {
                return None;
            }
        }

        let u = phi / self.phi_max;
        let v = p_hit.z / self.height;

        let dpdu = Vector3::new(-self.phi_max * p_hit.y, self.phi_max * p_hit.x, T::ZERO);
        let dpdv = Vector3::new(
            -p_hit.x / (T::ONE - v),
            -p_hit.y / (T::ONE - v),
            self.height,
        );

        let d2pdu2 = Vector3::new(p_hit.x, p_hit.y, T::ZERO) * -self.phi_max * self.phi_max;
        let d2pduv = Vector3::new(p_hit.y, -p_hit.x, T::ZERO) * self.phi_max / (T::ONE - v);

        let (dndu, dndv) = self.derive_normal(dpdu, dpdv, d2pdu2, d2pduv, Vector3::ZERO);

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
        .with_shape(*self)
            * self.object_to_world;

        Some((t_shape_hit.value(), intersection))
    }

    fn does_intersect(&self, ray: Ray<T>, _test_alpha: bool) -> bool {
        let (ray, (o_err, d_err)) = self.world_to_object.apply_err(ray);

        let ox = EFloat::new_with_err(ray.origin.x, o_err.x);
        let oy = EFloat::new_with_err(ray.origin.y, o_err.y);
        let oz = EFloat::new_with_err(ray.origin.z, o_err.z);
        let dx = EFloat::new_with_err(ray.direction.x, d_err.x);
        let dy = EFloat::new_with_err(ray.direction.y, d_err.y);
        let dz = EFloat::new_with_err(ray.direction.z, d_err.z);

        let k = EFloat::new(self.radius) / EFloat::new(self.height);
        let k = k * k;

        let a = dx * dx + dy * dy - k * dz * dz;
        let b = (dx * ox + dy * oy - k * dz * (oz - self.height)) * T::TWO;
        let c = ox * ox + oy * oy - k * (oz - self.height) * (oz - self.height);

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

        if p_hit.z < T::ZERO || p_hit.z > self.height || phi > self.phi_max {
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
            if p_hit.z < T::ZERO || p_hit.z > self.height || phi > self.phi_max {
                return false;
            }
        }

        true
    }

    fn area(&self) -> T {
        self.radius
            * ((self.height * self.height) + (self.radius * self.radius)).sqrt()
            * self.phi_max
            / T::TWO
    }
}

impl<T: Number> Deref for Cone<T> {
    type Target = ShapeData<T>;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl<T: Number> DerefMut for Cone<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}
