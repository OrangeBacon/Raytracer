use std::ops::{Deref, DerefMut};

use geometry::{
    Bounds3, ConstZero, EFloat, Number, PartialDerivatives, Point2, Point3, Ray,
    SurfaceInteractable, SurfaceInteraction, Transform, Vector3,
};

use crate::{Shape, ShapeData};

/// A hyperboloid
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Default)]
pub struct Hyperboloid<T: Number> {
    data: ShapeData<T>,
    point1: Point3<T>,
    point2: Point3<T>,
    z_min: T,
    z_max: T,
    phi_max: T,
    r_max: T,
    ah: T,
    ch: T,
}

impl<T: Number> Hyperboloid<T> {
    /// Create a new paraboloid instance
    pub fn new(
        object_to_world: Transform<T>,
        world_to_object: Transform<T>,
        reverse_orientation: bool,
        point1: Point3<T>,
        point2: Point3<T>,
        tm: T,
    ) -> Self {
        let radius_1 = point1.as_point2().to_vec().length();
        let radius_2 = point2.as_point2().to_vec().length();

        let mut this = Self {
            data: ShapeData::new(object_to_world, world_to_object, reverse_orientation),
            point1,
            point2,
            phi_max: tm.clamp(T::ZERO, T::cast(360)).to_radians(),
            z_min: point1.z.min(point2.z),
            z_max: point1.z.max(point2.z),
            r_max: radius_1.max(radius_2),
            ..Default::default()
        };

        if point2.z == T::ZERO {
            std::mem::swap(&mut this.point1, &mut this.point2);
        }

        let mut pp = this.point1;
        loop {
            pp += (this.point2 - this.point1) * T::TWO;
            let xy1 = pp.as_point2().to_vec().length_square();
            let xy2 = this.point2.as_point2().to_vec().length_square();
            this.ah = (T::ONE / xy1 - (pp.z - pp.z) / (xy1 * this.point2.z * this.point2.z))
                / (T::ONE - (xy2 * pp.z * pp.z) / (this.point2.z * this.point2.z));
            this.ch = (this.ah * xy2 - T::ONE) / (this.point2.z * this.point2.z);

            if this.ah.is_finite() && !this.ah.is_nan() {
                break;
            }
        }

        this
    }
}

impl<T: Number> Shape<T> for Hyperboloid<T> {
    fn object_bound(&self) -> geometry::Bounds3<T> {
        Bounds3::new(
            Point3::new(-self.r_max, -self.r_max, self.z_min),
            Point3::new(self.r_max, self.r_max, self.z_max),
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

        let a = dx * dx * self.ah + dy * dy * self.ah - dz * dz * self.ch;
        let b = (dx * ox * self.ah + dy * oy * self.ah - dz * oz * self.ch) * T::TWO;
        let c = ox * ox * self.ah + oy * oy * self.ah - oz * oz * self.ch - T::ONE;

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
        let mut v = (p_hit.z - self.point1.z) / (self.point2.z - self.point1.z);
        let pr = self.point1 * (T::ONE - v) + self.point2 * v;
        let mut phi = (pr.x * p_hit.y - p_hit.x * pr.y).atan2(p_hit.x * pr.x + p_hit.y * pr.y);
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
            v = (p_hit.z - self.point1.z) / (self.point2.z - self.point1.z);
            let pr = self.point1 * (T::ONE - v) + self.point2 * v;
            phi = (pr.x * p_hit.y - p_hit.x * pr.y).atan2(p_hit.x * pr.x + p_hit.y * pr.y);
            if phi < T::ZERO {
                phi += T::TWO * T::PI;
            }
            if p_hit.z < self.z_min || p_hit.z > self.z_max || phi > self.phi_max {
                return None;
            }
        }

        let u = phi / self.phi_max;

        let cos = phi.cos();
        let sin = phi.sin();
        let dpdu = Vector3::new(-self.phi_max * p_hit.y, self.phi_max * p_hit.x, T::ZERO);
        let dpdv = Vector3::new(
            (self.point2.x - self.point1.x) * cos - (self.point2.y - self.point1.y) * sin,
            (self.point2.x - self.point1.x) * sin - (self.point2.y - self.point1.y) * cos,
            self.point2.z - self.point1.z,
        ) * (self.z_max - self.z_min);

        let d2pdu2 = Vector3::new(p_hit.x, p_hit.y, T::ZERO) * -self.phi_max * self.phi_max;
        let d2pduv = Vector3::new(-dpdv.y, dpdv.x, T::ZERO) * self.phi_max;
        let d2pdv2 = Vector3::ZERO;

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

        let a = dx * dx * self.ah + dy * dy * self.ah - dz * dz * self.ch;
        let b = (dx * ox * self.ah + dy * oy * self.ah - dz * oz * self.ch) * T::TWO;
        let c = ox * ox * self.ah + oy * oy * self.ah - oz * oz * self.ch - T::ONE;

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
        let mut v = (p_hit.z - self.point1.z) / (self.point2.z - self.point1.z);
        let pr = self.point1 * (T::ONE - v) + self.point2 * v;
        let mut phi = (pr.x * p_hit.y - p_hit.x * pr.y).atan2(p_hit.x * pr.x + p_hit.y * pr.y);
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
            v = (p_hit.z - self.point1.z) / (self.point2.z - self.point1.z);
            let pr = self.point1 * (T::ONE - v) + self.point2 * v;
            phi = (pr.x * p_hit.y - p_hit.x * pr.y).atan2(p_hit.x * pr.x + p_hit.y * pr.y);
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
        let sqr = |a| a * a;
        let quad = |a| sqr(a) * sqr(a);
        let p1 = self.point1;
        let p2 = self.point2;

        self.phi_max / T::cast(6)
            * (T::TWO * quad(p1.x) - T::TWO * p1.x * p1.x * p1.x * p2.x
                + T::TWO * quad(p2.x)
                + T::TWO
                    * (p1.y * p1.y + p1.y * p2.y + p2.y * p2.y)
                    * (sqr(p1.y - p2.y) + sqr(p1.z - p2.z))
                + p2.x
                    * p2.x
                    * (T::cast(5) * p1.y * p1.y + T::TWO * p1.y * p2.y - T::cast(4) * p2.y * p2.y
                        + T::TWO * sqr(p1.z - p2.z))
                + p1.x
                    * p1.x
                    * (-T::cast(4) * p1.y * p1.y
                        + T::TWO * p1.y * p2.y
                        + T::cast(5) * p2.y * p2.y
                        + T::TWO * sqr(p1.z - p2.z))
                - T::TWO
                    * p1.x
                    * p2.x
                    * (p2.x * p2.x - p1.y * p1.y + T::cast(5) * p1.y * p2.y
                        - p2.y * p2.y
                        - p1.z * p1.z
                        + T::TWO * p1.z * p2.z
                        - p2.z * p2.z))
    }
}

impl<T: Number> SurfaceInteractable for Hyperboloid<T> {
    fn reverses_orientation(&self) -> bool {
        self.reverse_orientation ^ self.transform_swaps_handedness
    }
}

impl<T: Number> Deref for Hyperboloid<T> {
    type Target = ShapeData<T>;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl<T: Number> DerefMut for Hyperboloid<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}
