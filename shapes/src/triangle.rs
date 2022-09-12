use std::{
    any::TypeId,
    fmt::Debug,
    ops::{Deref, DerefMut},
    sync::Arc,
};

use geometry::{
    gamma, Bounds3, ConstZero, Normal3, Number, PartialDerivatives, Point2, Ray,
    SurfaceInteractable, SurfaceInteraction, Transform, Vector3,
};

use crate::{triangle_mesh::TriangleMesh, Shape, ShapeData};

/// A Single triangle referencing a triangle mesh
#[derive(Debug, Clone, PartialEq, PartialOrd, Default)]
pub struct Triangle<F: Number, T: Debug> {
    data: ShapeData<F>,
    mesh: Arc<TriangleMesh<F, T>>,
    idx: usize,
}

impl<F: Number, T: Debug> Triangle<F, T> {
    /// Create a new triangle from the data stored in a triangle mesh
    pub fn new(
        object_to_world: Transform<F>,
        world_to_object: Transform<F>,
        reverse_orientation: bool,
        mesh: Arc<TriangleMesh<F, T>>,
        tri_number: usize,
    ) -> Self {
        Self {
            data: ShapeData::new(object_to_world, world_to_object, reverse_orientation),
            mesh,
            idx: tri_number * 3,
        }
    }

    /// Create all the triangles represented by a triangle mesh
    pub fn from_mesh(
        mesh: Arc<TriangleMesh<F, T>>,
        world_to_object: Transform<F>,
        reverse_orientation: bool,
    ) -> Vec<Arc<dyn Shape<F>>>
    where
        T: 'static,
    {
        let tri_count = mesh.positions().len() / 3;

        let mut tris: Vec<Arc<dyn Shape<F>>> = Vec::with_capacity(tri_count);
        for i in 0..tris.capacity() {
            let tri = Triangle::new(
                *mesh.object_to_world(),
                world_to_object,
                reverse_orientation,
                Arc::clone(&mesh),
                i,
            );
            tris.push(Arc::new(tri));
        }

        tris
    }

    /// Get the UV coordinates for a given triangle
    fn uvs(&self) -> [Point2<F>; 3] {
        if let Some(uv) = self.mesh.uv() {
            [
                uv[self.mesh.indices()[self.idx]],
                uv[self.mesh.indices()[self.idx + 1]],
                uv[self.mesh.indices()[self.idx + 2]],
            ]
        } else {
            [
                Point2::ZERO,
                Point2::new(F::ONE, F::ZERO),
                Point2::new(F::ZERO, F::ONE),
            ]
        }
    }
}

impl<F: Number, T: Debug> Shape<F> for Triangle<F, T> {
    fn object_bound(&self) -> Bounds3<F> {
        let p0 = self.mesh.positions()[self.mesh.indices()[self.idx]];
        let p1 = self.mesh.positions()[self.mesh.indices()[self.idx + 1]];
        let p2 = self.mesh.positions()[self.mesh.indices()[self.idx + 2]];

        Bounds3::new(p0 * self.world_to_object, p1 * self.world_to_object)
            .union_point(p2 * self.world_to_object)
    }

    fn world_bound(&self) -> Bounds3<F> {
        let p0 = self.mesh.positions()[self.mesh.indices()[self.idx]];
        let p1 = self.mesh.positions()[self.mesh.indices()[self.idx + 1]];
        let p2 = self.mesh.positions()[self.mesh.indices()[self.idx + 2]];

        Bounds3::new(p0, p1).union_point(p2)
    }

    fn intersect(
        &self,
        ray: Ray<(), F>,
        test_alpha: bool,
    ) -> Option<(F, SurfaceInteraction<&dyn SurfaceInteractable, (), F>)> {
        let p0 = self.mesh.positions()[self.mesh.indices()[self.idx]];
        let p1 = self.mesh.positions()[self.mesh.indices()[self.idx + 1]];
        let p2 = self.mesh.positions()[self.mesh.indices()[self.idx + 2]];

        // transform triangle vertices into ray coordinate space
        let p0t = p0 - ray.origin.to_vec();
        let p1t = p1 - ray.origin.to_vec();
        let p2t = p2 - ray.origin.to_vec();

        let kz = ray.direction.abs().max_dimension();
        let kx = if kz + 1 == 3 { kz + 1 } else { 0 };
        let ky = if kx + 1 == 3 { kx + 1 } else { 0 };
        let idx = [kx, ky, kz];

        let d = ray.direction.permute(idx);
        let mut p0t = p0t.permute(idx);
        let mut p1t = p1t.permute(idx);
        let mut p2t = p2t.permute(idx);

        // apply shear transformation
        let sx = -d.x / d.z;
        let sy = -d.y / d.z;
        let sz = F::ONE / d.z;
        p0t.x += sx * p0t.z;
        p0t.y += sy * p0t.z;
        p1t.x += sx * p1t.z;
        p1t.y += sy * p1t.z;
        p2t.x += sx * p2t.z;
        p2t.y += sy * p2t.z;

        // apply edge function coefficients
        let mut e0 = p1t.x * p2t.y - p1t.y * p2t.x;
        let mut e1 = p2t.x * p0t.y - p2t.y * p0t.x;
        let mut e2 = p0t.x * p1t.y - p0t.y * p1t.x;

        // fallback to double precision on triangle edges
        if TypeId::of::<F>() == TypeId::of::<f32>()
            && (e0 == F::ZERO || e1 == F::ZERO || e2 == F::ZERO)
        {
            let p0t = p0t.cast::<f64>();
            let p1t = p1t.cast::<f64>();
            let p2t = p2t.cast::<f64>();

            e0 = F::cast(p1t.x * p2t.y - p1t.y * p2t.x);
            e1 = F::cast(p2t.x * p0t.y - p2t.y * p0t.x);
            e2 = F::cast(p0t.x * p1t.y - p0t.y * p1t.x);
        }

        // triangle edge and determinant tests
        if (e0 < F::ZERO || e1 < F::ZERO || e2 < F::ZERO)
            && (e0 > F::ZERO || e1 > F::ZERO || e2 > F::ZERO)
        {
            return None;
        }
        let det = e0 + e1 + e2;
        if det == F::ZERO {
            return None;
        }

        // compute scaled hit distance to triangle
        p0t.z *= sz;
        p1t.z *= sz;
        p2t.z *= sz;
        let scaled = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
        if (det < F::ZERO && (scaled >= F::ZERO || scaled < ray.t_max * det))
            || (det > F::ZERO && (scaled <= F::ZERO || scaled > ray.t_max * det))
        {
            return None;
        }

        // now we know there is a hit, so compute intersection t and barycentric coordinates
        let inv_det = F::ONE / det;
        let b0 = e0 * inv_det;
        let b1 = e1 * inv_det;
        let b2 = e2 * inv_det;
        let t = scaled * inv_det;

        // ensure triangle is greater than 0
        let max_zt = Vector3::new(p0t.z, p1t.z, p2t.z).abs().max_component();
        let delta_z = max_zt * gamma(3);

        let max_xt = Vector3::new(p0t.x, p1t.x, p2t.x).abs().max_component();
        let max_yt = Vector3::new(p0t.y, p1t.y, p2t.y).abs().max_component();
        let delta_x = (max_xt + max_zt) * gamma(5);
        let delta_y = (max_yt + max_zt) * gamma(5);

        let delta_e = F::TWO * (max_xt * max_yt * gamma(2) + delta_y * max_xt + delta_x * max_yt);
        let max_e = Vector3::new(e0, e1, e2).abs().max_component();
        let delta_t = F::cast(3)
            * (max_e * max_zt * gamma(3) + delta_e * max_zt + delta_z * max_e)
            * inv_det.abs();
        if t <= delta_t {
            return None;
        }

        // compute partial derivatives
        let uv = self.uvs();
        let duv02 = uv[0] - uv[2];
        let duv12 = uv[1] - uv[2];
        let dp02 = p0 - p2;
        let dp12 = p1 - p2;

        let det = duv02.x * duv12.y - duv02.y * duv12.x;
        let (dpdu, dpdv) = if det == F::ZERO {
            (p2 - p0).cross(p1 - p0).normalise().coordinate_system()
        } else {
            let inv_det = F::ONE / det;
            (
                (dp02 * duv12.y - dp12 * duv02.y) * inv_det,
                (dp02 * -duv12.x + dp12 * duv02.x) * inv_det,
            )
        };

        let p_hit = p0 * b0 + p1 * b1 + p2 * b2;
        let uv_hit = uv[0] * b0 + uv[1] * b1 + uv[2] * b2;

        if test_alpha {
            if let Some(alpha) = self.mesh.alpha_mask() {
                let isect_local = SurfaceInteraction::new(
                    p_hit,
                    Vector3::ZERO,
                    uv_hit,
                    Vector3::ZERO,
                    PartialDerivatives {
                        dpdu,
                        dpdv,
                        dndu: Normal3::ZERO,
                        dndv: Normal3::ZERO,
                    },
                    ray.time,
                )
                .with_shape(self as &dyn SurfaceInteractable);
                todo!(
                    "Implement textures here pls: {:?}, {:?}",
                    alpha,
                    isect_local
                )
            }
        }

        let p_error =
            ((p0 * b0).abs() + (p1 * b1).abs() + (p2 * b1).abs()).to_vec() * gamma::<F>(7);

        let mut isect = SurfaceInteraction::new(
            p_hit,
            p_error,
            uv_hit,
            -ray.direction,
            PartialDerivatives {
                dpdu,
                dpdv,
                dndu: Normal3::ZERO,
                dndv: Normal3::ZERO,
            },
            ray.time,
        )
        .with_shape(self as &dyn SurfaceInteractable);

        let normal = dp02.cross(dp12).normalise().to_normal();
        isect.interaction.normal = Some(normal);
        isect.shading.normal = normal;

        if self.mesh.normal().is_some() || self.mesh.tangent().is_some() {
            let ns = if let Some(norm) = self.mesh.normal() {
                (norm[self.mesh.indices()[self.idx]] * b0
                    + norm[self.mesh.indices()[self.idx + 1]] * b1
                    + norm[self.mesh.indices()[self.idx + 2]] * b2)
                    .normalise()
            } else {
                normal
            };

            let mut ss = if let Some(tan) = self.mesh.tangent() {
                (tan[self.mesh.indices()[self.idx]] * b0
                    + tan[self.mesh.indices()[self.idx + 1]] * b1
                    + tan[self.mesh.indices()[self.idx + 2]] * b2)
                    .normalise()
            } else {
                isect.derivatives.dpdu.normalise()
            };

            // compute bitangent vector
            let mut ts = ns.to_vector().cross(ss);
            if ts.length_square() > F::ZERO {
                ts = ts.normalise();
                ss = ts.cross(ns.to_vector());
            } else {
                (ss, ts) = ns.to_vector().coordinate_system();
            }

            let (dndu, dndv) = if let Some(norm) = self.mesh.normal() {
                let duv02 = uv[0] - uv[2];
                let duv12 = uv[1] - uv[2];
                let dn1 =
                    norm[self.mesh.indices()[self.idx]] - norm[self.mesh.indices()[self.idx + 2]];
                let dn2 = norm[self.mesh.indices()[self.idx + 1]]
                    - norm[self.mesh.indices()[self.idx + 2]];

                let det = duv02.x * duv12.y - duv02.y * duv12.x;
                if det == F::ZERO {
                    (Normal3::ZERO, Normal3::ZERO)
                } else {
                    let inv_det = F::ONE / det;
                    (
                        (dn1 * duv12.y - dn2 * duv02.y) * inv_det,
                        (dn1 * -duv12.x + dn2 * duv02.x) * inv_det,
                    )
                }
            } else {
                (Normal3::ZERO, Normal3::ZERO)
            };

            isect.set_shading_geometry(
                PartialDerivatives {
                    dpdu: ss,
                    dpdv: ts,
                    dndu,
                    dndv,
                },
                true,
            );
        }

        if self.mesh.normal().is_some() {
            isect.interaction.normal = Some(normal.face_forward(isect.shading.normal.to_vector()));
        } else if self.reverses_orientation() {
            isect.interaction.normal = Some(-normal);
            isect.shading.normal = -normal;
        }

        Some((t, isect))
    }

    fn does_intersect(&self, ray: Ray<(), F>, test_alpha: bool) -> bool {
        let p0 = self.mesh.positions()[self.mesh.indices()[self.idx]];
        let p1 = self.mesh.positions()[self.mesh.indices()[self.idx + 1]];
        let p2 = self.mesh.positions()[self.mesh.indices()[self.idx + 2]];

        // transform triangle vertices into ray coordinate space
        let p0t = p0 - ray.origin.to_vec();
        let p1t = p1 - ray.origin.to_vec();
        let p2t = p2 - ray.origin.to_vec();

        let kz = ray.direction.abs().max_dimension();
        let kx = if kz + 1 == 3 { kz + 1 } else { 0 };
        let ky = if kx + 1 == 3 { kx + 1 } else { 0 };
        let idx = [kx, ky, kz];

        let d = ray.direction.permute(idx);
        let mut p0t = p0t.permute(idx);
        let mut p1t = p1t.permute(idx);
        let mut p2t = p2t.permute(idx);

        // apply shear transformation
        let sx = -d.x / d.z;
        let sy = -d.y / d.z;
        let sz = F::ONE / d.z;
        p0t.x += sx * p0t.z;
        p0t.y += sy * p0t.z;
        p1t.x += sx * p1t.z;
        p1t.y += sy * p1t.z;
        p2t.x += sx * p2t.z;
        p2t.y += sy * p2t.z;

        // apply edge function coefficients
        let mut e0 = p1t.x * p2t.y - p1t.y * p2t.x;
        let mut e1 = p2t.x * p0t.y - p2t.y * p0t.x;
        let mut e2 = p0t.x * p1t.y - p0t.y * p1t.x;

        // fallback to double precision on triangle edges
        if TypeId::of::<F>() == TypeId::of::<f32>()
            && (e0 == F::ZERO || e1 == F::ZERO || e2 == F::ZERO)
        {
            let p0t = p0t.cast::<f64>();
            let p1t = p1t.cast::<f64>();
            let p2t = p2t.cast::<f64>();

            e0 = F::cast(p1t.x * p2t.y - p1t.y * p2t.x);
            e1 = F::cast(p2t.x * p0t.y - p2t.y * p0t.x);
            e2 = F::cast(p0t.x * p1t.y - p0t.y * p1t.x);
        }

        // triangle edge and determinant tests
        if (e0 < F::ZERO || e1 < F::ZERO || e2 < F::ZERO)
            && (e0 > F::ZERO || e1 > F::ZERO || e2 > F::ZERO)
        {
            return false;
        }
        let det = e0 + e1 + e2;
        if det == F::ZERO {
            return false;
        }

        // compute scaled hit distance to triangle
        p0t.z *= sz;
        p1t.z *= sz;
        p2t.z *= sz;
        let scaled = e0 * p0t.z + e1 * p1t.z + e2 * p2t.z;
        if (det < F::ZERO && (scaled >= F::ZERO || scaled < ray.t_max * det))
            || (det > F::ZERO && (scaled <= F::ZERO || scaled > ray.t_max * det))
        {
            return false;
        }

        // now we know there is a hit, so compute intersection t and barycentric coordinates
        let inv_det = F::ONE / det;
        let b0 = e0 * inv_det;
        let b1 = e1 * inv_det;
        let b2 = e2 * inv_det;

        // compute partial derivatives
        let uv = self.uvs();
        let duv02 = uv[0] - uv[2];
        let duv12 = uv[1] - uv[2];
        let dp02 = p0 - p2;
        let dp12 = p1 - p2;

        let det = duv02.x * duv12.y - duv02.y * duv12.x;
        let (dpdu, dpdv) = if det == F::ZERO {
            (p2 - p0).cross(p1 - p0).normalise().coordinate_system()
        } else {
            let inv_det = F::ONE / det;
            (
                (dp02 * duv12.y - dp12 * duv02.y) * inv_det,
                (dp02 * -duv12.x + dp12 * duv02.x) * inv_det,
            )
        };

        let p_hit = p0 * b0 + p1 * b1 + p2 * b2;
        let uv_hit = uv[0] * b0 + uv[1] * b1 + uv[2] * b2;

        if test_alpha {
            if let Some(alpha) = self.mesh.alpha_mask() {
                let isect_local = SurfaceInteraction::new(
                    p_hit,
                    Vector3::ZERO,
                    uv_hit,
                    Vector3::ZERO,
                    PartialDerivatives {
                        dpdu,
                        dpdv,
                        dndu: Normal3::ZERO,
                        dndv: Normal3::ZERO,
                    },
                    ray.time,
                )
                .with_shape(self as &dyn SurfaceInteractable);
                todo!(
                    "Implement textures here pls: {:?}, {:?}",
                    alpha,
                    isect_local
                )
            }
        }

        true
    }

    fn area(&self) -> F {
        let p0 = self.mesh.positions()[self.mesh.indices()[self.idx]];
        let p1 = self.mesh.positions()[self.mesh.indices()[self.idx + 1]];
        let p2 = self.mesh.positions()[self.mesh.indices()[self.idx + 2]];

        (F::HALF) * (p1 - p0).cross(p2 - p0).length()
    }
}

impl<F: Number, T: Debug> SurfaceInteractable for Triangle<F, T> {
    fn reverses_orientation(&self) -> bool {
        self.reverse_orientation ^ self.transform_swaps_handedness
    }
}

impl<F: Number, T: Debug> Deref for Triangle<F, T> {
    type Target = ShapeData<F>;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl<F: Number, T: Debug> DerefMut for Triangle<F, T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}
