use geometry::{EFloat, Normal3, Number, Ray, Vector3};

use crate::ShapeData;

impl<T: Number> ShapeData<T> {
    /// Calculate t_shape_hit and t0 for the quadratic intersection between the given
    /// ray and description of quadric surface.
    pub fn quadric_coefficients(&self, radius: T, ray: Ray<T>) -> Option<(EFloat<T>, EFloat<T>)> {
        let (ray, (o_err, d_err)) = self.world_to_object.apply_err(ray);

        let ox = EFloat::new_with_err(ray.origin.x, o_err.x);
        let oy = EFloat::new_with_err(ray.origin.y, o_err.y);
        let oz = EFloat::new_with_err(ray.origin.z, o_err.z);
        let dx = EFloat::new_with_err(ray.direction.x, d_err.x);
        let dy = EFloat::new_with_err(ray.direction.y, d_err.y);
        let dz = EFloat::new_with_err(ray.direction.z, d_err.z);

        let a = dx * dx + dy * dy + dz * dz;
        let b = (dx * ox + dy * oy + dz * oz) * T::TWO;
        let c = ox * ox + oy * oy + oz * oz - EFloat::new(radius) * EFloat::new(radius);

        let (t0, t1) = EFloat::quadratic(a, b, c)?;

        if t0.upper_bound() > ray.t_max || t1.lower_bound() <= T::ZERO {
            return None;
        }

        if t0.lower_bound() <= T::ZERO {
            if t1.upper_bound() > ray.t_max {
                return None;
            }
            Some((t1, t1))
        } else {
            Some((t0, t1))
        }
    }

    /// Calculate the partial derivatives of the normal vector, dndu and dndv
    /// using the Weingarten equations
    pub fn derive_normal(
        &self,
        dpdu: Vector3<T>,
        dpdv: Vector3<T>,
        d2pdu2: Vector3<T>,
        d2pduv: Vector3<T>,
        d2pdv2: Vector3<T>,
    ) -> (Normal3<T>, Normal3<T>) {
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

        (dndu, dndv)
    }
}
