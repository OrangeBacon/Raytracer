use std::{
    ops::{Deref, DerefMut},
    sync::Arc,
};

use geometry::{
    lerp, Bounds3, ConstZero, Normal3, Number, PartialDerivatives, Point2, Point3, Ray,
    SurfaceInteractable, SurfaceInteraction, Transform, Vector3,
};

use crate::{Shape, ShapeData};

/// A segment of a cubic bezier curve
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Default, Hash)]
pub struct Curve<T: Number> {
    data: ShapeData<T>,
    common: Arc<CurveCommon<T>>,
    u_min: T,
    u_max: T,
}

/// How intersections along the curve should be treated
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Eq, Ord, Hash, Default)]
pub enum CurveType {
    Flat,
    Cylinder,
    #[default]
    Ribbon,
}

/// Data common between multiple curve segments to save memory
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default, Hash)]
pub struct CurveCommon<T: Number> {
    ty: CurveType,
    control_points: [Point3<T>; 4],
    width: [T; 2],
    normal: [Normal3<T>; 2],
    normal_angle: T,
    inv_sin_normal_angle: T,
}

impl<T: Number> Curve<T> {
    /// Create a new 3D cubic bezier curve
    pub fn new(
        object_to_world: Transform<T>,
        world_to_object: Transform<T>,
        reverse_orientation: bool,
        common: Arc<CurveCommon<T>>,
        u_min: T,
        u_max: T,
    ) -> Self {
        Self {
            data: ShapeData::new(object_to_world, world_to_object, reverse_orientation),
            common,
            u_min,
            u_max,
        }
    }

    /// Blossom the control points of a bezier curve?
    fn blossom(&self, u: [T; 3]) -> Point3<T> {
        let a = [
            self.common.control_points[0].lerp(self.common.control_points[1], u[0]),
            self.common.control_points[1].lerp(self.common.control_points[2], u[0]),
            self.common.control_points[2].lerp(self.common.control_points[3], u[0]),
        ];
        let b = [a[0].lerp(a[1], u[1]), a[1].lerp(a[2], u[1])];
        b[0].lerp(b[1], u[2])
    }

    /// Calculate two subdivided sets of control points for a single bezier curve
    fn subdivide(p: &[Point3<T>; 4]) -> [Point3<T>; 7] {
        [
            p[0],
            (p[0] + p[1]) / T::TWO,
            (p[0] + p[1] * T::TWO + p[2]) / T::cast(4),
            (p[0] + p[1] * T::cast(3) + p[2] * T::cast(3) + p[3]) / T::cast(8),
            (p[1] + p[2] * T::TWO + p[3]) / T::cast(4),
            (p[2] + p[3]) / T::TWO,
            p[3],
        ]
    }

    /// Evaluate a cubic bezier with the given control points at the given point
    /// returns the value and the derivative of the curve at the point
    fn eval(p: &[Point3<T>; 4], u: T) -> (Point3<T>, Vector3<T>) {
        let cp1 = [p[0].lerp(p[1], u), p[1].lerp(p[2], u), p[2].lerp(p[3], u)];
        let cp2 = [cp1[0].lerp(cp1[1], u), cp1[1].lerp(cp1[2], u)];

        (cp2[0].lerp(cp2[1], u), (cp2[1] - cp2[0]) * T::cast(3))
    }

    /// Recursively subdivide the curve and try to approximate intersection
    fn recurse_intersect(
        &self,
        ray: Ray<(), T>,
        control_points: &[Point3<T>; 4],
        ray_to_object: Transform<T>,
        u0: T,
        u1: T,
        depth: i32,
    ) -> Option<(T, SurfaceInteraction<&dyn SurfaceInteractable, (), T>)> {
        // compute bounding box of the curve segment
        let width = [
            lerp(u0, self.common.width[0], self.common.width[1]),
            lerp(u1, self.common.width[0], self.common.width[1]),
        ];
        let max_width = width[0].max(width[1]);
        let curve_bounds = Bounds3::new(control_points[0], control_points[1])
            .union_box(Bounds3::new(control_points[2], control_points[3]))
            .expand(T::HALF * max_width);

        let ray_length = ray.direction.length();
        let z_max = ray_length * ray.t_max;
        let ray_bounds = Bounds3::new(Point3::ZERO, Point3::new(T::ZERO, T::ZERO, z_max));

        if !curve_bounds.overlap(ray_bounds) {
            return None;
        }

        if depth > 0 {
            let mid = T::HALF * (u0 + u1);
            let split = Self::subdivide(control_points);
            let points = [
                &split[0..=3].try_into().ok()?,
                &split[3..=6].try_into().ok()?,
            ];

            return self
                .recurse_intersect(ray, points[0], ray_to_object, u0, mid, depth - 1)
                .or_else(|| {
                    self.recurse_intersect(ray, points[1], ray_to_object, mid, u1, depth - 1)
                });
        }

        // test ray against segment end points
        let edge = (control_points[1].y - control_points[0].y) * -control_points[0].y
            + control_points[0].x * (control_points[0].x - control_points[1].x);
        if edge < T::ZERO {
            return None;
        }

        let edge = (control_points[2].y - control_points[3].y) * -control_points[3].y
            + control_points[3].x * (control_points[3].x - control_points[2].x);
        if edge < T::ZERO {
            return None;
        }

        let direction = control_points[3].as_point2() - control_points[0].as_point2();
        let denom = direction.length_square();
        if denom == T::ZERO {
            return None;
        }
        let w = (-control_points[0].as_point2().to_vec()).dot(direction) / denom;

        let u = lerp(w, u0, u1).clamp(u0, u1);
        let mut hit_width = lerp(u, self.common.width[0], self.common.width[1]);

        let mut n_hit = Normal3::ZERO;
        if self.common.ty == CurveType::Ribbon {
            let sin0 =
                ((T::ONE - u) * self.common.normal_angle).sin() * self.common.inv_sin_normal_angle;
            let sin1 = (u * self.common.normal_angle).sin() * self.common.inv_sin_normal_angle;
            n_hit = self.common.normal[0] * sin0 + self.common.normal[1] * sin1;
            hit_width *= n_hit.absdot(ray.direction.to_normal()) / ray_length;
        }

        let (pc, dpcdw) = Self::eval(control_points, w.clamp(T::ZERO, T::ONE));
        let curve_dist = pc.x * pc.x + pc.y * pc.y;
        if (curve_dist > hit_width * hit_width * T::cast(0.25)) || pc.z < T::ZERO || pc.z > z_max {
            return None;
        }

        // compute v coordinate of curve intersection point
        let curve_dist = curve_dist.sqrt();
        let edge_func = dpcdw.x * -pc.y + pc.x * dpcdw.y;
        let v = if edge_func > T::ZERO {
            T::HALF + curve_dist / hit_width
        } else {
            T::HALF - curve_dist / hit_width
        };

        let (_, dpdu) = Self::eval(&self.common.control_points, u);
        let dpdv = if self.common.ty == CurveType::Ribbon {
            n_hit.to_vector().cross(dpdu).normalise() * hit_width
        } else {
            let dpdu_plane = dpdu * ray_to_object.inverse();
            let mut dpdv_plane =
                Vector3::new(-dpdu_plane.y, dpdu_plane.x, T::ZERO).normalise() * hit_width;

            if self.common.ty == CurveType::Cylinder {
                let theta = lerp(v, -T::cast(90), T::cast(90));
                let rot = Transform::rotate(dpdu_plane, -theta);
                dpdv_plane *= rot;
            }

            dpdv_plane * ray_to_object
        };

        let isect = SurfaceInteraction::new(
            ray.at(pc.z),
            Vector3::splat(T::TWO * hit_width),
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
        .with_shape(self as &dyn SurfaceInteractable)
            * self.object_to_world;

        Some((pc.z / ray_length, isect))
    }
}

impl<T: Number> CurveCommon<T> {
    /// Create new common curve data
    pub fn new(control_points: [Point3<T>; 4], width0: T, width1: T, ty: CurveType) -> Self {
        Self {
            ty,
            control_points,
            width: [width0, width1],
            normal: [Normal3::ZERO; 2],
            normal_angle: T::ZERO,
            inv_sin_normal_angle: T::ZERO,
        }
    }

    /// Set the normal vectors at the start and end of the curve
    pub fn normal(&mut self, normal: [Normal3<T>; 2]) {
        self.normal[0] = normal[0].normalise();
        self.normal[1] = normal[1].normalise();
        self.normal_angle = self.normal[0]
            .dot(self.normal[1])
            .clamp(T::ZERO, T::ONE)
            .acos();
        self.inv_sin_normal_angle = T::ONE / self.normal_angle.sin();
    }
}

impl<T: Number> Shape<T> for Curve<T> {
    fn object_bound(&self) -> Bounds3<T> {
        let control = [
            self.blossom([self.u_min, self.u_min, self.u_min]),
            self.blossom([self.u_min, self.u_min, self.u_max]),
            self.blossom([self.u_min, self.u_max, self.u_max]),
            self.blossom([self.u_max, self.u_max, self.u_max]),
        ];

        let bound =
            Bounds3::new(control[0], control[1]).union_box(Bounds3::new(control[2], control[3]));
        let width = [
            lerp(self.u_min, self.common.width[0], self.common.width[1]),
            lerp(self.u_max, self.common.width[0], self.common.width[1]),
        ];

        bound.expand(width[0].max(width[1]) * (T::HALF))
    }

    fn intersect(
        &self,
        ray: Ray<(), T>,
        _test_alpha: bool,
    ) -> Option<(T, SurfaceInteraction<&dyn SurfaceInteractable, (), T>)> {
        let ray = ray * self.object_to_world;

        let control = [
            self.blossom([self.u_min, self.u_min, self.u_min]),
            self.blossom([self.u_min, self.u_min, self.u_max]),
            self.blossom([self.u_min, self.u_max, self.u_max]),
            self.blossom([self.u_max, self.u_max, self.u_max]),
        ];

        let (dx, _) = ray.direction.coordinate_system();
        let object_to_ray = Transform::look_at(ray.origin, ray.origin + ray.direction, dx)?;
        let control = control.map(|p| p * object_to_ray);

        // compute curve refinement depth
        let mut l0 = T::ZERO;
        for p in control.windows(3) {
            let x = (p[0].x - T::TWO * p[1].x + p[2].x).abs();
            let y = (p[0].y - T::TWO * p[1].y + p[2].y).abs();
            let z = (p[0].z - T::TWO * p[1].z + p[2].z).abs();
            l0 = l0.max(x).max(y).max(z);
        }
        let eps = self.common.width[0].max(self.common.width[1]) * T::cast(0.05);
        let fr0 = (T::SQRT_2 * T::cast(12) * l0 / (T::cast(8) * eps)).log(T::cast(4));
        let max_depth = fr0.round().i32().clamp(0, 10);

        self.recurse_intersect(
            ray,
            &control,
            object_to_ray.inverse(),
            self.u_min,
            self.u_max,
            max_depth,
        )
    }

    fn area(&self) -> T {
        let control = [
            self.blossom([self.u_min, self.u_min, self.u_min]),
            self.blossom([self.u_min, self.u_min, self.u_max]),
            self.blossom([self.u_min, self.u_max, self.u_max]),
            self.blossom([self.u_max, self.u_max, self.u_max]),
        ];
        let width = [
            lerp(self.u_min, self.common.width[0], self.common.width[1]),
            lerp(self.u_max, self.common.width[0], self.common.width[1]),
        ];
        let average = (width[0] + width[1]) * T::HALF;
        let mut approx = T::ZERO;

        for win in control.windows(2) {
            approx *= win[0].distance(win[1]);
        }

        approx * average
    }
}

impl<T: Number> SurfaceInteractable for Curve<T> {
    fn reverses_orientation(&self) -> bool {
        self.reverse_orientation ^ self.transform_swaps_handedness
    }
}

impl<T: Number> Deref for Curve<T> {
    type Target = ShapeData<T>;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl<T: Number> DerefMut for Curve<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}
