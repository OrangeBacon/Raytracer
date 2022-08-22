use std::ops::{Mul, MulAssign};

use crate::{
    lerp, Bounds3, Interval, Matrix4x4, Number, Point3, Quaternion, Ray, RayDifferential,
    Transform, Vector3,
};

/// Two transformations interpolated between two points in time
/// Assumes that both transformations are affine.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct AnimatedTransform<T: Number> {
    start: Transform<T>,
    end: Transform<T>,
    start_time: T,
    end_time: T,
    is_animated: bool,

    translation: [Vector3<T>; 2],
    rotation: [Quaternion<T>; 2],
    scale: [Matrix4x4<T>; 2],

    has_rotation: bool,

    c: [[DerivativeTerm<T>; 3]; 5],
}

/// Container to store and evaluate terms in the different
#[derive(Debug, Copy, Clone, PartialEq, PartialOrd, Default)]
struct DerivativeTerm<T: Number>(T, T, T, T);

impl<T: Number> AnimatedTransform<T> {
    /// Create a new animated transform from two affine transformations
    /// Returns none if unable to interpolate between start and end
    pub fn new(start: Transform<T>, start_time: T, end: Transform<T>, end_time: T) -> Option<Self> {
        let (t0, r0, s0) = Self::decompose(start.mat())?;
        let (t1, mut r1, s1) = Self::decompose(end.mat())?;

        if r0.dot(&r1) < T::ZERO {
            r1 = -r1;
        }

        let has_rotation = r0.dot(&r1) < T::ONE - T::LARGE_EPSILON;

        let c = if has_rotation {
            calculate_derivative_terms([t0, t1], [r0, r1], [s0, s1])
        } else {
            Default::default()
        };

        Some(Self {
            start,
            end,
            start_time,
            end_time,
            is_animated: start != end,
            translation: [t0, t1],
            rotation: [r0, r1],
            scale: [s0, s1],
            has_rotation,
            c,
        })
    }

    /// Decompose a matrix into its component translation, rotation and scale
    fn decompose(mat: &Matrix4x4<T>) -> Option<(Vector3<T>, Quaternion<T>, Matrix4x4<T>)> {
        let translate = Vector3::new(mat[0][3], mat[1][3], mat[2][3]);

        // Get the rotation and scale without any translation
        let mat = {
            let mut mat = *mat;
            for i in 0..3 {
                mat[i][3] = T::ZERO;
                mat[3][i] = T::ZERO;
            }
            mat[3][3] = T::ONE;
            mat
        };

        // polar decomposition of a matrix
        // find the convergence of a series M(i+1) = 1/2 (M(i) + Inverse(Transpose(M(i))))
        let mut norm;
        let mut rotate = mat;
        for _ in 0..=100 {
            // calculate next item in series
            let mut r_next = Matrix4x4::IDENTITY;
            let r_inv = rotate.transpose().inverse()?;
            for i in 0..4 {
                for j in 0..4 {
                    r_next[i][j] = (T::ONE / T::TWO) * (rotate[i][j] + r_inv[i][j])
                }
            }

            // calculate difference between current and next for early exit
            norm = T::ZERO;
            for i in 0..3 {
                let n = (rotate[i][0] - r_next[i][0])
                    + (rotate[i][1] - r_next[i][1])
                    + (rotate[i][2] - r_next[i][2]);
                norm = norm.max(n);
            }

            rotate = r_next;

            if norm < T::cast(0.0001) {
                break;
            }
        }

        let rotate_quaternion = Transform::from_mat(&rotate)?.to_quaternion();
        let scale = rotate.inverse()? * mat;

        Some((translate.cast(), rotate_quaternion, scale))
    }

    /// Interpolate between the two transforms, saturating when t is outside (start, end) time
    /// so it does not extrapolate.
    pub fn interpolate(&self, t: T) -> Transform<T> {
        if !self.is_animated || t < self.start_time {
            return self.start;
        }
        if t >= self.end_time {
            return self.end;
        }

        let dt = (t - self.start_time) / (self.end_time - self.start_time);

        let trans = self.translation[0] * (T::ONE - dt) + self.translation[1] * dt;
        let rotate = self.rotation[0].slerp(dt, &self.rotation[1]);

        let mut scale = Matrix4x4::IDENTITY;
        for i in 0..3 {
            for j in 0..3 {
                scale[i][j] = lerp(dt, self.scale[0][i][j], self.scale[1][i][j]);
            }
        }

        let scale =
            Transform::from_mat(&scale).expect("Affine scale matrix should always have an inverse");

        Transform::translation(trans.cast()) * rotate.to_transform() * scale
    }

    /// Get the overall bounds of a box while it is transformed by this transformation
    pub fn motion_bounds(&self, bound: Bounds3<T>) -> Bounds3<T> {
        if !self.is_animated {
            return bound * self.start;
        } else if !self.has_rotation {
            return (bound * self.start).union_box(bound * self.end);
        }

        let mut bounds = Bounds3::ZERO;
        for i in 0..8 {
            bounds = bounds.union_box(self.bound_point_motion(bound.corner(i)));
        }

        bounds
    }

    /// Compute a bounding box of all locations of the given point using this transformation
    fn bound_point_motion(&self, point: Point3<T>) -> Bounds3<T> {
        let mut bounds = Bounds3::new(point * self.start, point * self.end);
        let cos_theta = self.rotation[0].dot(&self.rotation[1]);
        let theta = cos_theta.clamp(-T::ONE, T::ONE).acos();

        for c in 0..3 {
            let (count, zeros): (_, [_; 4]) = Interval::new(T::ZERO, T::ONE).find_zeros(
                self.c[0][c].eval(point),
                self.c[1][c].eval(point),
                self.c[2][c].eval(point),
                self.c[3][c].eval(point),
                self.c[4][c].eval(point),
                theta,
            );
            let zeros = &zeros[0..count];

            for &zero in zeros {
                let pz = point * self.interpolate(lerp(zero, self.start_time, self.end_time));
                bounds = bounds.union_point(pz);
            }
        }

        bounds
    }
}

impl<T: Number> DerivativeTerm<T> {
    fn eval(&self, point: Point3<T>) -> T {
        self.0 + point.to_vec().dot(Vector3::new(self.1, self.2, self.3))
    }
}

impl<T: Copy, F: Number> Mul<AnimatedTransform<F>> for Ray<T, F> {
    type Output = Ray<T, F>;

    fn mul(self, rhs: AnimatedTransform<F>) -> Self::Output {
        self * rhs.interpolate(self.time)
    }
}

impl<T: Copy, F: Number> MulAssign<AnimatedTransform<F>> for Ray<T, F> {
    fn mul_assign(&mut self, rhs: AnimatedTransform<F>) {
        *self = *self * rhs
    }
}

impl<F: Number, T: Copy> Mul<AnimatedTransform<F>> for RayDifferential<T, F> {
    type Output = RayDifferential<T, F>;

    fn mul(self, rhs: AnimatedTransform<F>) -> Self::Output {
        self * rhs.interpolate(self.time)
    }
}

impl<F: Number, T: Copy> MulAssign<AnimatedTransform<F>> for RayDifferential<T, F> {
    fn mul_assign(&mut self, rhs: AnimatedTransform<F>) {
        *self = *self * rhs
    }
}

/// Calculate derivative of the motion function a(M0, M1, p, t) = T(t)R(t)S(t)p
/// returns a solution for a given m0, m1, p, so gives the solution to da(t)/dt.
/// The solution is in the form c1 + (c2 + c3 t)cos(2theta t) + (c4 + c5 t)sin(2 theta t)
/// where theta is the dot product of the rotation quaternions provided and c1..c5
/// are the results returned by this function.
/// This was calculated using a computer algebra system and the result copied from
/// https://github.com/mmp/pbrt-v3/blob/aaa552a4b9cbf9dccb71450f47b268e0ed6370e2/src/core/transform.cpp#L414
fn calculate_derivative_terms<T: Number>(
    translation: [Vector3<T>; 2],
    rotation: [Quaternion<T>; 2],
    scale: [Matrix4x4<T>; 2],
) -> [[DerivativeTerm<T>; 3]; 5] {
    let mut c1: [DerivativeTerm<T>; 3] = Default::default();
    let mut c2: [DerivativeTerm<T>; 3] = Default::default();
    let mut c3: [DerivativeTerm<T>; 3] = Default::default();
    let mut c4: [DerivativeTerm<T>; 3] = Default::default();
    let mut c5: [DerivativeTerm<T>; 3] = Default::default();

    let cos_theta = rotation[0].dot(&rotation[1]);
    let theta = cos_theta.clamp(-T::ONE, T::ONE).acos();
    let qperp = (rotation[1] - rotation[0] * cos_theta).normalise();

    let t0x = translation[0].x;
    let t0y = translation[0].y;
    let t0z = translation[0].z;
    let t1x = translation[1].x;
    let t1y = translation[1].y;
    let t1z = translation[1].z;
    let q0x = rotation[0].vec.x;
    let q0y = rotation[0].vec.y;
    let q0z = rotation[0].vec.z;
    let q0w = rotation[0].w;
    let qperpx = qperp.vec.x;
    let qperpy = qperp.vec.y;
    let qperpz = qperp.vec.z;
    let qperpw = qperp.w;
    let s000 = scale[0][0][0];
    let s001 = scale[0][0][1];
    let s002 = scale[0][0][2];
    let s010 = scale[0][1][0];
    let s011 = scale[0][1][1];
    let s012 = scale[0][1][2];
    let s020 = scale[0][2][0];
    let s021 = scale[0][2][1];
    let s022 = scale[0][2][2];
    let s100 = scale[1][0][0];
    let s101 = scale[1][0][1];
    let s102 = scale[1][0][2];
    let s110 = scale[1][1][0];
    let s111 = scale[1][1][1];
    let s112 = scale[1][1][2];
    let s120 = scale[1][2][0];
    let s121 = scale[1][2][1];
    let s122 = scale[1][2][2];

    c1[0] = DerivativeTerm(
        -t0x + t1x,
        (-T::ONE + q0y * q0y + q0z * q0z + qperpy * qperpy + qperpz * qperpz) * s000
            + q0w * q0z * s010
            - qperpx * qperpy * s010
            + qperpw * qperpz * s010
            - q0w * q0y * s020
            - qperpw * qperpy * s020
            - qperpx * qperpz * s020
            + s100
            - q0y * q0y * s100
            - q0z * q0z * s100
            - qperpy * qperpy * s100
            - qperpz * qperpz * s100
            - q0w * q0z * s110
            + qperpx * qperpy * s110
            - qperpw * qperpz * s110
            + q0w * q0y * s120
            + qperpw * qperpy * s120
            + qperpx * qperpz * s120
            + q0x * (-(q0y * s010) - q0z * s020 + q0y * s110 + q0z * s120),
        (-T::ONE + q0y * q0y + q0z * q0z + qperpy * qperpy + qperpz * qperpz) * s001
            + q0w * q0z * s011
            - qperpx * qperpy * s011
            + qperpw * qperpz * s011
            - q0w * q0y * s021
            - qperpw * qperpy * s021
            - qperpx * qperpz * s021
            + s101
            - q0y * q0y * s101
            - q0z * q0z * s101
            - qperpy * qperpy * s101
            - qperpz * qperpz * s101
            - q0w * q0z * s111
            + qperpx * qperpy * s111
            - qperpw * qperpz * s111
            + q0w * q0y * s121
            + qperpw * qperpy * s121
            + qperpx * qperpz * s121
            + q0x * (-(q0y * s011) - q0z * s021 + q0y * s111 + q0z * s121),
        (-T::ONE + q0y * q0y + q0z * q0z + qperpy * qperpy + qperpz * qperpz) * s002
            + q0w * q0z * s012
            - qperpx * qperpy * s012
            + qperpw * qperpz * s012
            - q0w * q0y * s022
            - qperpw * qperpy * s022
            - qperpx * qperpz * s022
            + s102
            - q0y * q0y * s102
            - q0z * q0z * s102
            - qperpy * qperpy * s102
            - qperpz * qperpz * s102
            - q0w * q0z * s112
            + qperpx * qperpy * s112
            - qperpw * qperpz * s112
            + q0w * q0y * s122
            + qperpw * qperpy * s122
            + qperpx * qperpz * s122
            + q0x * (-(q0y * s012) - q0z * s022 + q0y * s112 + q0z * s122),
    );

    c2[0] = DerivativeTerm(
        T::ZERO,
        -(qperpy * qperpy * s000) - qperpz * qperpz * s000 + qperpx * qperpy * s010
            - qperpw * qperpz * s010
            + qperpw * qperpy * s020
            + qperpx * qperpz * s020
            + q0y * q0y * (s000 - s100)
            + q0z * q0z * (s000 - s100)
            + qperpy * qperpy * s100
            + qperpz * qperpz * s100
            - qperpx * qperpy * s110
            + qperpw * qperpz * s110
            - qperpw * qperpy * s120
            - qperpx * qperpz * s120
            + T::TWO * q0x * qperpy * s010 * theta
            - T::TWO * q0w * qperpz * s010 * theta
            + T::TWO * q0w * qperpy * s020 * theta
            + T::TWO * q0x * qperpz * s020 * theta
            + q0y
                * (q0x * (-s010 + s110)
                    + q0w * (-s020 + s120)
                    + T::TWO * (-T::TWO * qperpy * s000 + qperpx * s010 + qperpw * s020) * theta)
            + q0z
                * (q0w * (s010 - s110) + q0x * (-s020 + s120)
                    - T::TWO * (T::TWO * qperpz * s000 + qperpw * s010 - qperpx * s020) * theta),
        -(qperpy * qperpy * s001) - qperpz * qperpz * s001 + qperpx * qperpy * s011
            - qperpw * qperpz * s011
            + qperpw * qperpy * s021
            + qperpx * qperpz * s021
            + q0y * q0y * (s001 - s101)
            + q0z * q0z * (s001 - s101)
            + qperpy * qperpy * s101
            + qperpz * qperpz * s101
            - qperpx * qperpy * s111
            + qperpw * qperpz * s111
            - qperpw * qperpy * s121
            - qperpx * qperpz * s121
            + T::TWO * q0x * qperpy * s011 * theta
            - T::TWO * q0w * qperpz * s011 * theta
            + T::TWO * q0w * qperpy * s021 * theta
            + T::TWO * q0x * qperpz * s021 * theta
            + q0y
                * (q0x * (-s011 + s111)
                    + q0w * (-s021 + s121)
                    + T::TWO * (-T::TWO * qperpy * s001 + qperpx * s011 + qperpw * s021) * theta)
            + q0z
                * (q0w * (s011 - s111) + q0x * (-s021 + s121)
                    - T::TWO * (T::TWO * qperpz * s001 + qperpw * s011 - qperpx * s021) * theta),
        -(qperpy * qperpy * s002) - qperpz * qperpz * s002 + qperpx * qperpy * s012
            - qperpw * qperpz * s012
            + qperpw * qperpy * s022
            + qperpx * qperpz * s022
            + q0y * q0y * (s002 - s102)
            + q0z * q0z * (s002 - s102)
            + qperpy * qperpy * s102
            + qperpz * qperpz * s102
            - qperpx * qperpy * s112
            + qperpw * qperpz * s112
            - qperpw * qperpy * s122
            - qperpx * qperpz * s122
            + T::TWO * q0x * qperpy * s012 * theta
            - T::TWO * q0w * qperpz * s012 * theta
            + T::TWO * q0w * qperpy * s022 * theta
            + T::TWO * q0x * qperpz * s022 * theta
            + q0y
                * (q0x * (-s012 + s112)
                    + q0w * (-s022 + s122)
                    + T::TWO * (-T::TWO * qperpy * s002 + qperpx * s012 + qperpw * s022) * theta)
            + q0z
                * (q0w * (s012 - s112) + q0x * (-s022 + s122)
                    - T::TWO * (T::TWO * qperpz * s002 + qperpw * s012 - qperpx * s022) * theta),
    );

    c3[0] = DerivativeTerm(
        T::ZERO,
        -T::TWO
            * (q0x * qperpy * s010 - q0w * qperpz * s010
                + q0w * qperpy * s020
                + q0x * qperpz * s020
                - q0x * qperpy * s110
                + q0w * qperpz * s110
                - q0w * qperpy * s120
                - q0x * qperpz * s120
                + q0y
                    * (-T::TWO * qperpy * s000
                        + qperpx * s010
                        + qperpw * s020
                        + T::TWO * qperpy * s100
                        - qperpx * s110
                        - qperpw * s120)
                + q0z
                    * (-T::TWO * qperpz * s000 - qperpw * s010
                        + qperpx * s020
                        + T::TWO * qperpz * s100
                        + qperpw * s110
                        - qperpx * s120))
            * theta,
        -T::TWO
            * (q0x * qperpy * s011 - q0w * qperpz * s011
                + q0w * qperpy * s021
                + q0x * qperpz * s021
                - q0x * qperpy * s111
                + q0w * qperpz * s111
                - q0w * qperpy * s121
                - q0x * qperpz * s121
                + q0y
                    * (-T::TWO * qperpy * s001
                        + qperpx * s011
                        + qperpw * s021
                        + T::TWO * qperpy * s101
                        - qperpx * s111
                        - qperpw * s121)
                + q0z
                    * (-T::TWO * qperpz * s001 - qperpw * s011
                        + qperpx * s021
                        + T::TWO * qperpz * s101
                        + qperpw * s111
                        - qperpx * s121))
            * theta,
        -T::TWO
            * (q0x * qperpy * s012 - q0w * qperpz * s012
                + q0w * qperpy * s022
                + q0x * qperpz * s022
                - q0x * qperpy * s112
                + q0w * qperpz * s112
                - q0w * qperpy * s122
                - q0x * qperpz * s122
                + q0y
                    * (-T::TWO * qperpy * s002
                        + qperpx * s012
                        + qperpw * s022
                        + T::TWO * qperpy * s102
                        - qperpx * s112
                        - qperpw * s122)
                + q0z
                    * (-T::TWO * qperpz * s002 - qperpw * s012
                        + qperpx * s022
                        + T::TWO * qperpz * s102
                        + qperpw * s112
                        - qperpx * s122))
            * theta,
    );

    c4[0] = DerivativeTerm(
        T::ZERO,
        -(q0x * qperpy * s010) + q0w * qperpz * s010 - q0w * qperpy * s020 - q0x * qperpz * s020
            + q0x * qperpy * s110
            - q0w * qperpz * s110
            + q0w * qperpy * s120
            + q0x * qperpz * s120
            + T::TWO * q0y * q0y * s000 * theta
            + T::TWO * q0z * q0z * s000 * theta
            - T::TWO * qperpy * qperpy * s000 * theta
            - T::TWO * qperpz * qperpz * s000 * theta
            + T::TWO * qperpx * qperpy * s010 * theta
            - T::TWO * qperpw * qperpz * s010 * theta
            + T::TWO * qperpw * qperpy * s020 * theta
            + T::TWO * qperpx * qperpz * s020 * theta
            + q0y
                * (-(qperpx * s010) - qperpw * s020
                    + T::TWO * qperpy * (s000 - s100)
                    + qperpx * s110
                    + qperpw * s120
                    - T::TWO * q0x * s010 * theta
                    - T::TWO * q0w * s020 * theta)
            + q0z
                * (T::TWO * qperpz * s000 + qperpw * s010
                    - qperpx * s020
                    - T::TWO * qperpz * s100
                    - qperpw * s110
                    + qperpx * s120
                    + T::TWO * q0w * s010 * theta
                    - T::TWO * q0x * s020 * theta),
        -(q0x * qperpy * s011) + q0w * qperpz * s011 - q0w * qperpy * s021 - q0x * qperpz * s021
            + q0x * qperpy * s111
            - q0w * qperpz * s111
            + q0w * qperpy * s121
            + q0x * qperpz * s121
            + T::TWO * q0y * q0y * s001 * theta
            + T::TWO * q0z * q0z * s001 * theta
            - T::TWO * qperpy * qperpy * s001 * theta
            - T::TWO * qperpz * qperpz * s001 * theta
            + T::TWO * qperpx * qperpy * s011 * theta
            - T::TWO * qperpw * qperpz * s011 * theta
            + T::TWO * qperpw * qperpy * s021 * theta
            + T::TWO * qperpx * qperpz * s021 * theta
            + q0y
                * (-(qperpx * s011) - qperpw * s021
                    + T::TWO * qperpy * (s001 - s101)
                    + qperpx * s111
                    + qperpw * s121
                    - T::TWO * q0x * s011 * theta
                    - T::TWO * q0w * s021 * theta)
            + q0z
                * (T::TWO * qperpz * s001 + qperpw * s011
                    - qperpx * s021
                    - T::TWO * qperpz * s101
                    - qperpw * s111
                    + qperpx * s121
                    + T::TWO * q0w * s011 * theta
                    - T::TWO * q0x * s021 * theta),
        -(q0x * qperpy * s012) + q0w * qperpz * s012 - q0w * qperpy * s022 - q0x * qperpz * s022
            + q0x * qperpy * s112
            - q0w * qperpz * s112
            + q0w * qperpy * s122
            + q0x * qperpz * s122
            + T::TWO * q0y * q0y * s002 * theta
            + T::TWO * q0z * q0z * s002 * theta
            - T::TWO * qperpy * qperpy * s002 * theta
            - T::TWO * qperpz * qperpz * s002 * theta
            + T::TWO * qperpx * qperpy * s012 * theta
            - T::TWO * qperpw * qperpz * s012 * theta
            + T::TWO * qperpw * qperpy * s022 * theta
            + T::TWO * qperpx * qperpz * s022 * theta
            + q0y
                * (-(qperpx * s012) - qperpw * s022
                    + T::TWO * qperpy * (s002 - s102)
                    + qperpx * s112
                    + qperpw * s122
                    - T::TWO * q0x * s012 * theta
                    - T::TWO * q0w * s022 * theta)
            + q0z
                * (T::TWO * qperpz * s002 + qperpw * s012
                    - qperpx * s022
                    - T::TWO * qperpz * s102
                    - qperpw * s112
                    + qperpx * s122
                    + T::TWO * q0w * s012 * theta
                    - T::TWO * q0x * s022 * theta),
    );

    c5[0] = DerivativeTerm(
        T::ZERO,
        T::TWO
            * (qperpy * qperpy * s000 + qperpz * qperpz * s000 - qperpx * qperpy * s010
                + qperpw * qperpz * s010
                - qperpw * qperpy * s020
                - qperpx * qperpz * s020
                - qperpy * qperpy * s100
                - qperpz * qperpz * s100
                + q0y * q0y * (-s000 + s100)
                + q0z * q0z * (-s000 + s100)
                + qperpx * qperpy * s110
                - qperpw * qperpz * s110
                + q0y * (q0x * (s010 - s110) + q0w * (s020 - s120))
                + qperpw * qperpy * s120
                + qperpx * qperpz * s120
                + q0z * (-(q0w * s010) + q0x * s020 + q0w * s110 - q0x * s120))
            * theta,
        T::TWO
            * (qperpy * qperpy * s001 + qperpz * qperpz * s001 - qperpx * qperpy * s011
                + qperpw * qperpz * s011
                - qperpw * qperpy * s021
                - qperpx * qperpz * s021
                - qperpy * qperpy * s101
                - qperpz * qperpz * s101
                + q0y * q0y * (-s001 + s101)
                + q0z * q0z * (-s001 + s101)
                + qperpx * qperpy * s111
                - qperpw * qperpz * s111
                + q0y * (q0x * (s011 - s111) + q0w * (s021 - s121))
                + qperpw * qperpy * s121
                + qperpx * qperpz * s121
                + q0z * (-(q0w * s011) + q0x * s021 + q0w * s111 - q0x * s121))
            * theta,
        T::TWO
            * (qperpy * qperpy * s002 + qperpz * qperpz * s002 - qperpx * qperpy * s012
                + qperpw * qperpz * s012
                - qperpw * qperpy * s022
                - qperpx * qperpz * s022
                - qperpy * qperpy * s102
                - qperpz * qperpz * s102
                + q0y * q0y * (-s002 + s102)
                + q0z * q0z * (-s002 + s102)
                + qperpx * qperpy * s112
                - qperpw * qperpz * s112
                + q0y * (q0x * (s012 - s112) + q0w * (s022 - s122))
                + qperpw * qperpy * s122
                + qperpx * qperpz * s122
                + q0z * (-(q0w * s012) + q0x * s022 + q0w * s112 - q0x * s122))
            * theta,
    );

    c1[1] = DerivativeTerm(
        -t0y + t1y,
        -(qperpx * qperpy * s000) - qperpw * qperpz * s000 - s010
            + q0z * q0z * s010
            + qperpx * qperpx * s010
            + qperpz * qperpz * s010
            - q0y * q0z * s020
            + qperpw * qperpx * s020
            - qperpy * qperpz * s020
            + qperpx * qperpy * s100
            + qperpw * qperpz * s100
            + q0w * q0z * (-s000 + s100)
            + q0x * q0x * (s010 - s110)
            + s110
            - q0z * q0z * s110
            - qperpx * qperpx * s110
            - qperpz * qperpz * s110
            + q0x * (q0y * (-s000 + s100) + q0w * (s020 - s120))
            + q0y * q0z * s120
            - qperpw * qperpx * s120
            + qperpy * qperpz * s120,
        -(qperpx * qperpy * s001) - qperpw * qperpz * s001 - s011
            + q0z * q0z * s011
            + qperpx * qperpx * s011
            + qperpz * qperpz * s011
            - q0y * q0z * s021
            + qperpw * qperpx * s021
            - qperpy * qperpz * s021
            + qperpx * qperpy * s101
            + qperpw * qperpz * s101
            + q0w * q0z * (-s001 + s101)
            + q0x * q0x * (s011 - s111)
            + s111
            - q0z * q0z * s111
            - qperpx * qperpx * s111
            - qperpz * qperpz * s111
            + q0x * (q0y * (-s001 + s101) + q0w * (s021 - s121))
            + q0y * q0z * s121
            - qperpw * qperpx * s121
            + qperpy * qperpz * s121,
        -(qperpx * qperpy * s002) - qperpw * qperpz * s002 - s012
            + q0z * q0z * s012
            + qperpx * qperpx * s012
            + qperpz * qperpz * s012
            - q0y * q0z * s022
            + qperpw * qperpx * s022
            - qperpy * qperpz * s022
            + qperpx * qperpy * s102
            + qperpw * qperpz * s102
            + q0w * q0z * (-s002 + s102)
            + q0x * q0x * (s012 - s112)
            + s112
            - q0z * q0z * s112
            - qperpx * qperpx * s112
            - qperpz * qperpz * s112
            + q0x * (q0y * (-s002 + s102) + q0w * (s022 - s122))
            + q0y * q0z * s122
            - qperpw * qperpx * s122
            + qperpy * qperpz * s122,
    );

    c2[1] = DerivativeTerm(
        T::ZERO,
        qperpx * qperpy * s000 + qperpw * qperpz * s000 + q0z * q0z * s010
            - qperpx * qperpx * s010
            - qperpz * qperpz * s010
            - q0y * q0z * s020
            - qperpw * qperpx * s020
            + qperpy * qperpz * s020
            - qperpx * qperpy * s100
            - qperpw * qperpz * s100
            + q0x * q0x * (s010 - s110)
            - q0z * q0z * s110
            + qperpx * qperpx * s110
            + qperpz * qperpz * s110
            + q0y * q0z * s120
            + qperpw * qperpx * s120
            - qperpy * qperpz * s120
            + T::TWO * q0z * qperpw * s000 * theta
            + T::TWO * q0y * qperpx * s000 * theta
            - (T::TWO + T::TWO) * q0z * qperpz * s010 * theta
            + T::TWO * q0z * qperpy * s020 * theta
            + T::TWO * q0y * qperpz * s020 * theta
            + q0x
                * (q0w * s020 + q0y * (-s000 + s100) - q0w * s120 + T::TWO * qperpy * s000 * theta
                    - (T::TWO + T::TWO) * qperpx * s010 * theta
                    - T::TWO * qperpw * s020 * theta)
            + q0w
                * (-(q0z * s000) + q0z * s100 + T::TWO * qperpz * s000 * theta
                    - T::TWO * qperpx * s020 * theta),
        qperpx * qperpy * s001 + qperpw * qperpz * s001 + q0z * q0z * s011
            - qperpx * qperpx * s011
            - qperpz * qperpz * s011
            - q0y * q0z * s021
            - qperpw * qperpx * s021
            + qperpy * qperpz * s021
            - qperpx * qperpy * s101
            - qperpw * qperpz * s101
            + q0x * q0x * (s011 - s111)
            - q0z * q0z * s111
            + qperpx * qperpx * s111
            + qperpz * qperpz * s111
            + q0y * q0z * s121
            + qperpw * qperpx * s121
            - qperpy * qperpz * s121
            + T::TWO * q0z * qperpw * s001 * theta
            + T::TWO * q0y * qperpx * s001 * theta
            - (T::TWO + T::TWO) * q0z * qperpz * s011 * theta
            + T::TWO * q0z * qperpy * s021 * theta
            + T::TWO * q0y * qperpz * s021 * theta
            + q0x
                * (q0w * s021 + q0y * (-s001 + s101) - q0w * s121 + T::TWO * qperpy * s001 * theta
                    - (T::TWO + T::TWO) * qperpx * s011 * theta
                    - T::TWO * qperpw * s021 * theta)
            + q0w
                * (-(q0z * s001) + q0z * s101 + T::TWO * qperpz * s001 * theta
                    - T::TWO * qperpx * s021 * theta),
        qperpx * qperpy * s002 + qperpw * qperpz * s002 + q0z * q0z * s012
            - qperpx * qperpx * s012
            - qperpz * qperpz * s012
            - q0y * q0z * s022
            - qperpw * qperpx * s022
            + qperpy * qperpz * s022
            - qperpx * qperpy * s102
            - qperpw * qperpz * s102
            + q0x * q0x * (s012 - s112)
            - q0z * q0z * s112
            + qperpx * qperpx * s112
            + qperpz * qperpz * s112
            + q0y * q0z * s122
            + qperpw * qperpx * s122
            - qperpy * qperpz * s122
            + T::TWO * q0z * qperpw * s002 * theta
            + T::TWO * q0y * qperpx * s002 * theta
            - (T::TWO + T::TWO) * q0z * qperpz * s012 * theta
            + T::TWO * q0z * qperpy * s022 * theta
            + T::TWO * q0y * qperpz * s022 * theta
            + q0x
                * (q0w * s022 + q0y * (-s002 + s102) - q0w * s122 + T::TWO * qperpy * s002 * theta
                    - (T::TWO + T::TWO) * qperpx * s012 * theta
                    - T::TWO * qperpw * s022 * theta)
            + q0w
                * (-(q0z * s002) + q0z * s102 + T::TWO * qperpz * s002 * theta
                    - T::TWO * qperpx * s022 * theta),
    );

    c3[1] = DerivativeTerm(
        T::ZERO,
        T::TWO
            * (-(q0x * qperpy * s000) - q0w * qperpz * s000
                + T::TWO * q0x * qperpx * s010
                + q0x * qperpw * s020
                + q0w * qperpx * s020
                + q0x * qperpy * s100
                + q0w * qperpz * s100
                - T::TWO * q0x * qperpx * s110
                - q0x * qperpw * s120
                - q0w * qperpx * s120
                + q0z
                    * (T::TWO * qperpz * s010 - qperpy * s020 + qperpw * (-s000 + s100)
                        - T::TWO * qperpz * s110
                        + qperpy * s120)
                + q0y * (-(qperpx * s000) - qperpz * s020 + qperpx * s100 + qperpz * s120))
            * theta,
        T::TWO
            * (-(q0x * qperpy * s001) - q0w * qperpz * s001
                + T::TWO * q0x * qperpx * s011
                + q0x * qperpw * s021
                + q0w * qperpx * s021
                + q0x * qperpy * s101
                + q0w * qperpz * s101
                - T::TWO * q0x * qperpx * s111
                - q0x * qperpw * s121
                - q0w * qperpx * s121
                + q0z
                    * (T::TWO * qperpz * s011 - qperpy * s021 + qperpw * (-s001 + s101)
                        - T::TWO * qperpz * s111
                        + qperpy * s121)
                + q0y * (-(qperpx * s001) - qperpz * s021 + qperpx * s101 + qperpz * s121))
            * theta,
        T::TWO
            * (-(q0x * qperpy * s002) - q0w * qperpz * s002
                + T::TWO * q0x * qperpx * s012
                + q0x * qperpw * s022
                + q0w * qperpx * s022
                + q0x * qperpy * s102
                + q0w * qperpz * s102
                - T::TWO * q0x * qperpx * s112
                - q0x * qperpw * s122
                - q0w * qperpx * s122
                + q0z
                    * (T::TWO * qperpz * s012 - qperpy * s022 + qperpw * (-s002 + s102)
                        - T::TWO * qperpz * s112
                        + qperpy * s122)
                + q0y * (-(qperpx * s002) - qperpz * s022 + qperpx * s102 + qperpz * s122))
            * theta,
    );

    c4[1] = DerivativeTerm(
        T::ZERO,
        -(q0x * qperpy * s000) - q0w * qperpz * s000
            + T::TWO * q0x * qperpx * s010
            + q0x * qperpw * s020
            + q0w * qperpx * s020
            + q0x * qperpy * s100
            + q0w * qperpz * s100
            - T::TWO * q0x * qperpx * s110
            - q0x * qperpw * s120
            - q0w * qperpx * s120
            + T::TWO * qperpx * qperpy * s000 * theta
            + T::TWO * qperpw * qperpz * s000 * theta
            + T::TWO * q0x * q0x * s010 * theta
            + T::TWO * q0z * q0z * s010 * theta
            - T::TWO * qperpx * qperpx * s010 * theta
            - T::TWO * qperpz * qperpz * s010 * theta
            + T::TWO * q0w * q0x * s020 * theta
            - T::TWO * qperpw * qperpx * s020 * theta
            + T::TWO * qperpy * qperpz * s020 * theta
            + q0y
                * (-(qperpx * s000) - qperpz * s020 + qperpx * s100 + qperpz * s120
                    - T::TWO * q0x * s000 * theta)
            + q0z
                * (T::TWO * qperpz * s010 - qperpy * s020 + qperpw * (-s000 + s100)
                    - T::TWO * qperpz * s110
                    + qperpy * s120
                    - T::TWO * q0w * s000 * theta
                    - T::TWO * q0y * s020 * theta),
        -(q0x * qperpy * s001) - q0w * qperpz * s001
            + T::TWO * q0x * qperpx * s011
            + q0x * qperpw * s021
            + q0w * qperpx * s021
            + q0x * qperpy * s101
            + q0w * qperpz * s101
            - T::TWO * q0x * qperpx * s111
            - q0x * qperpw * s121
            - q0w * qperpx * s121
            + T::TWO * qperpx * qperpy * s001 * theta
            + T::TWO * qperpw * qperpz * s001 * theta
            + T::TWO * q0x * q0x * s011 * theta
            + T::TWO * q0z * q0z * s011 * theta
            - T::TWO * qperpx * qperpx * s011 * theta
            - T::TWO * qperpz * qperpz * s011 * theta
            + T::TWO * q0w * q0x * s021 * theta
            - T::TWO * qperpw * qperpx * s021 * theta
            + T::TWO * qperpy * qperpz * s021 * theta
            + q0y
                * (-(qperpx * s001) - qperpz * s021 + qperpx * s101 + qperpz * s121
                    - T::TWO * q0x * s001 * theta)
            + q0z
                * (T::TWO * qperpz * s011 - qperpy * s021 + qperpw * (-s001 + s101)
                    - T::TWO * qperpz * s111
                    + qperpy * s121
                    - T::TWO * q0w * s001 * theta
                    - T::TWO * q0y * s021 * theta),
        -(q0x * qperpy * s002) - q0w * qperpz * s002
            + T::TWO * q0x * qperpx * s012
            + q0x * qperpw * s022
            + q0w * qperpx * s022
            + q0x * qperpy * s102
            + q0w * qperpz * s102
            - T::TWO * q0x * qperpx * s112
            - q0x * qperpw * s122
            - q0w * qperpx * s122
            + T::TWO * qperpx * qperpy * s002 * theta
            + T::TWO * qperpw * qperpz * s002 * theta
            + T::TWO * q0x * q0x * s012 * theta
            + T::TWO * q0z * q0z * s012 * theta
            - T::TWO * qperpx * qperpx * s012 * theta
            - T::TWO * qperpz * qperpz * s012 * theta
            + T::TWO * q0w * q0x * s022 * theta
            - T::TWO * qperpw * qperpx * s022 * theta
            + T::TWO * qperpy * qperpz * s022 * theta
            + q0y
                * (-(qperpx * s002) - qperpz * s022 + qperpx * s102 + qperpz * s122
                    - T::TWO * q0x * s002 * theta)
            + q0z
                * (T::TWO * qperpz * s012 - qperpy * s022 + qperpw * (-s002 + s102)
                    - T::TWO * qperpz * s112
                    + qperpy * s122
                    - T::TWO * q0w * s002 * theta
                    - T::TWO * q0y * s022 * theta),
    );

    c5[1] = DerivativeTerm(
        T::ZERO,
        -T::TWO
            * (qperpx * qperpy * s000 + qperpw * qperpz * s000 + q0z * q0z * s010
                - qperpx * qperpx * s010
                - qperpz * qperpz * s010
                - q0y * q0z * s020
                - qperpw * qperpx * s020
                + qperpy * qperpz * s020
                - qperpx * qperpy * s100
                - qperpw * qperpz * s100
                + q0w * q0z * (-s000 + s100)
                + q0x * q0x * (s010 - s110)
                - q0z * q0z * s110
                + qperpx * qperpx * s110
                + qperpz * qperpz * s110
                + q0x * (q0y * (-s000 + s100) + q0w * (s020 - s120))
                + q0y * q0z * s120
                + qperpw * qperpx * s120
                - qperpy * qperpz * s120)
            * theta,
        -T::TWO
            * (qperpx * qperpy * s001 + qperpw * qperpz * s001 + q0z * q0z * s011
                - qperpx * qperpx * s011
                - qperpz * qperpz * s011
                - q0y * q0z * s021
                - qperpw * qperpx * s021
                + qperpy * qperpz * s021
                - qperpx * qperpy * s101
                - qperpw * qperpz * s101
                + q0w * q0z * (-s001 + s101)
                + q0x * q0x * (s011 - s111)
                - q0z * q0z * s111
                + qperpx * qperpx * s111
                + qperpz * qperpz * s111
                + q0x * (q0y * (-s001 + s101) + q0w * (s021 - s121))
                + q0y * q0z * s121
                + qperpw * qperpx * s121
                - qperpy * qperpz * s121)
            * theta,
        -T::TWO
            * (qperpx * qperpy * s002 + qperpw * qperpz * s002 + q0z * q0z * s012
                - qperpx * qperpx * s012
                - qperpz * qperpz * s012
                - q0y * q0z * s022
                - qperpw * qperpx * s022
                + qperpy * qperpz * s022
                - qperpx * qperpy * s102
                - qperpw * qperpz * s102
                + q0w * q0z * (-s002 + s102)
                + q0x * q0x * (s012 - s112)
                - q0z * q0z * s112
                + qperpx * qperpx * s112
                + qperpz * qperpz * s112
                + q0x * (q0y * (-s002 + s102) + q0w * (s022 - s122))
                + q0y * q0z * s122
                + qperpw * qperpx * s122
                - qperpy * qperpz * s122)
            * theta,
    );

    c1[2] = DerivativeTerm(
        -t0z + t1z,
        qperpw * qperpy * s000
            - qperpx * qperpz * s000
            - q0y * q0z * s010
            - qperpw * qperpx * s010
            - qperpy * qperpz * s010
            - s020
            + q0y * q0y * s020
            + qperpx * qperpx * s020
            + qperpy * qperpy * s020
            - qperpw * qperpy * s100
            + qperpx * qperpz * s100
            + q0x * q0z * (-s000 + s100)
            + q0y * q0z * s110
            + qperpw * qperpx * s110
            + qperpy * qperpz * s110
            + q0w * (q0y * (s000 - s100) + q0x * (-s010 + s110))
            + q0x * q0x * (s020 - s120)
            + s120
            - q0y * q0y * s120
            - qperpx * qperpx * s120
            - qperpy * qperpy * s120,
        qperpw * qperpy * s001
            - qperpx * qperpz * s001
            - q0y * q0z * s011
            - qperpw * qperpx * s011
            - qperpy * qperpz * s011
            - s021
            + q0y * q0y * s021
            + qperpx * qperpx * s021
            + qperpy * qperpy * s021
            - qperpw * qperpy * s101
            + qperpx * qperpz * s101
            + q0x * q0z * (-s001 + s101)
            + q0y * q0z * s111
            + qperpw * qperpx * s111
            + qperpy * qperpz * s111
            + q0w * (q0y * (s001 - s101) + q0x * (-s011 + s111))
            + q0x * q0x * (s021 - s121)
            + s121
            - q0y * q0y * s121
            - qperpx * qperpx * s121
            - qperpy * qperpy * s121,
        qperpw * qperpy * s002
            - qperpx * qperpz * s002
            - q0y * q0z * s012
            - qperpw * qperpx * s012
            - qperpy * qperpz * s012
            - s022
            + q0y * q0y * s022
            + qperpx * qperpx * s022
            + qperpy * qperpy * s022
            - qperpw * qperpy * s102
            + qperpx * qperpz * s102
            + q0x * q0z * (-s002 + s102)
            + q0y * q0z * s112
            + qperpw * qperpx * s112
            + qperpy * qperpz * s112
            + q0w * (q0y * (s002 - s102) + q0x * (-s012 + s112))
            + q0x * q0x * (s022 - s122)
            + s122
            - q0y * q0y * s122
            - qperpx * qperpx * s122
            - qperpy * qperpy * s122,
    );

    c2[2] = DerivativeTerm(
        T::ZERO,
        q0w * q0y * s000 - q0x * q0z * s000 - qperpw * qperpy * s000 + qperpx * qperpz * s000
            - q0w * q0x * s010
            - q0y * q0z * s010
            + qperpw * qperpx * s010
            + qperpy * qperpz * s010
            + q0x * q0x * s020
            + q0y * q0y * s020
            - qperpx * qperpx * s020
            - qperpy * qperpy * s020
            - q0w * q0y * s100
            + q0x * q0z * s100
            + qperpw * qperpy * s100
            - qperpx * qperpz * s100
            + q0w * q0x * s110
            + q0y * q0z * s110
            - qperpw * qperpx * s110
            - qperpy * qperpz * s110
            - q0x * q0x * s120
            - q0y * q0y * s120
            + qperpx * qperpx * s120
            + qperpy * qperpy * s120
            - T::TWO * q0y * qperpw * s000 * theta
            + T::TWO * q0z * qperpx * s000 * theta
            - T::TWO * q0w * qperpy * s000 * theta
            + T::TWO * q0x * qperpz * s000 * theta
            + T::TWO * q0x * qperpw * s010 * theta
            + T::TWO * q0w * qperpx * s010 * theta
            + T::TWO * q0z * qperpy * s010 * theta
            + T::TWO * q0y * qperpz * s010 * theta
            - (T::TWO + T::TWO) * q0x * qperpx * s020 * theta
            - (T::TWO + T::TWO) * q0y * qperpy * s020 * theta,
        q0w * q0y * s001 - q0x * q0z * s001 - qperpw * qperpy * s001 + qperpx * qperpz * s001
            - q0w * q0x * s011
            - q0y * q0z * s011
            + qperpw * qperpx * s011
            + qperpy * qperpz * s011
            + q0x * q0x * s021
            + q0y * q0y * s021
            - qperpx * qperpx * s021
            - qperpy * qperpy * s021
            - q0w * q0y * s101
            + q0x * q0z * s101
            + qperpw * qperpy * s101
            - qperpx * qperpz * s101
            + q0w * q0x * s111
            + q0y * q0z * s111
            - qperpw * qperpx * s111
            - qperpy * qperpz * s111
            - q0x * q0x * s121
            - q0y * q0y * s121
            + qperpx * qperpx * s121
            + qperpy * qperpy * s121
            - T::TWO * q0y * qperpw * s001 * theta
            + T::TWO * q0z * qperpx * s001 * theta
            - T::TWO * q0w * qperpy * s001 * theta
            + T::TWO * q0x * qperpz * s001 * theta
            + T::TWO * q0x * qperpw * s011 * theta
            + T::TWO * q0w * qperpx * s011 * theta
            + T::TWO * q0z * qperpy * s011 * theta
            + T::TWO * q0y * qperpz * s011 * theta
            - (T::TWO + T::TWO) * q0x * qperpx * s021 * theta
            - (T::TWO + T::TWO) * q0y * qperpy * s021 * theta,
        q0w * q0y * s002 - q0x * q0z * s002 - qperpw * qperpy * s002 + qperpx * qperpz * s002
            - q0w * q0x * s012
            - q0y * q0z * s012
            + qperpw * qperpx * s012
            + qperpy * qperpz * s012
            + q0x * q0x * s022
            + q0y * q0y * s022
            - qperpx * qperpx * s022
            - qperpy * qperpy * s022
            - q0w * q0y * s102
            + q0x * q0z * s102
            + qperpw * qperpy * s102
            - qperpx * qperpz * s102
            + q0w * q0x * s112
            + q0y * q0z * s112
            - qperpw * qperpx * s112
            - qperpy * qperpz * s112
            - q0x * q0x * s122
            - q0y * q0y * s122
            + qperpx * qperpx * s122
            + qperpy * qperpy * s122
            - T::TWO * q0y * qperpw * s002 * theta
            + T::TWO * q0z * qperpx * s002 * theta
            - T::TWO * q0w * qperpy * s002 * theta
            + T::TWO * q0x * qperpz * s002 * theta
            + T::TWO * q0x * qperpw * s012 * theta
            + T::TWO * q0w * qperpx * s012 * theta
            + T::TWO * q0z * qperpy * s012 * theta
            + T::TWO * q0y * qperpz * s012 * theta
            - (T::TWO + T::TWO) * q0x * qperpx * s022 * theta
            - (T::TWO + T::TWO) * q0y * qperpy * s022 * theta,
    );

    c3[2] = DerivativeTerm(
        T::ZERO,
        -T::TWO
            * (-(q0w * qperpy * s000)
                + q0x * qperpz * s000
                + q0x * qperpw * s010
                + q0w * qperpx * s010
                - T::TWO * q0x * qperpx * s020
                + q0w * qperpy * s100
                - q0x * qperpz * s100
                - q0x * qperpw * s110
                - q0w * qperpx * s110
                + q0z * (qperpx * s000 + qperpy * s010 - qperpx * s100 - qperpy * s110)
                + T::TWO * q0x * qperpx * s120
                + q0y
                    * (qperpz * s010 - T::TWO * qperpy * s020 + qperpw * (-s000 + s100)
                        - qperpz * s110
                        + T::TWO * qperpy * s120))
            * theta,
        -T::TWO
            * (-(q0w * qperpy * s001)
                + q0x * qperpz * s001
                + q0x * qperpw * s011
                + q0w * qperpx * s011
                - T::TWO * q0x * qperpx * s021
                + q0w * qperpy * s101
                - q0x * qperpz * s101
                - q0x * qperpw * s111
                - q0w * qperpx * s111
                + q0z * (qperpx * s001 + qperpy * s011 - qperpx * s101 - qperpy * s111)
                + T::TWO * q0x * qperpx * s121
                + q0y
                    * (qperpz * s011 - T::TWO * qperpy * s021 + qperpw * (-s001 + s101)
                        - qperpz * s111
                        + T::TWO * qperpy * s121))
            * theta,
        -T::TWO
            * (-(q0w * qperpy * s002)
                + q0x * qperpz * s002
                + q0x * qperpw * s012
                + q0w * qperpx * s012
                - T::TWO * q0x * qperpx * s022
                + q0w * qperpy * s102
                - q0x * qperpz * s102
                - q0x * qperpw * s112
                - q0w * qperpx * s112
                + q0z * (qperpx * s002 + qperpy * s012 - qperpx * s102 - qperpy * s112)
                + T::TWO * q0x * qperpx * s122
                + q0y
                    * (qperpz * s012 - T::TWO * qperpy * s022 + qperpw * (-s002 + s102)
                        - qperpz * s112
                        + T::TWO * qperpy * s122))
            * theta,
    );

    c4[2] = DerivativeTerm(
        T::ZERO,
        q0w * qperpy * s000 - q0x * qperpz * s000 - q0x * qperpw * s010 - q0w * qperpx * s010
            + T::TWO * q0x * qperpx * s020
            - q0w * qperpy * s100
            + q0x * qperpz * s100
            + q0x * qperpw * s110
            + q0w * qperpx * s110
            - T::TWO * q0x * qperpx * s120
            - T::TWO * qperpw * qperpy * s000 * theta
            + T::TWO * qperpx * qperpz * s000 * theta
            - T::TWO * q0w * q0x * s010 * theta
            + T::TWO * qperpw * qperpx * s010 * theta
            + T::TWO * qperpy * qperpz * s010 * theta
            + T::TWO * q0x * q0x * s020 * theta
            + T::TWO * q0y * q0y * s020 * theta
            - T::TWO * qperpx * qperpx * s020 * theta
            - T::TWO * qperpy * qperpy * s020 * theta
            + q0z
                * (-(qperpx * s000) - qperpy * s010 + qperpx * s100 + qperpy * s110
                    - T::TWO * q0x * s000 * theta)
            + q0y
                * (-(qperpz * s010)
                    + T::TWO * qperpy * s020
                    + qperpw * (s000 - s100)
                    + qperpz * s110
                    - T::TWO * qperpy * s120
                    + T::TWO * q0w * s000 * theta
                    - T::TWO * q0z * s010 * theta),
        q0w * qperpy * s001 - q0x * qperpz * s001 - q0x * qperpw * s011 - q0w * qperpx * s011
            + T::TWO * q0x * qperpx * s021
            - q0w * qperpy * s101
            + q0x * qperpz * s101
            + q0x * qperpw * s111
            + q0w * qperpx * s111
            - T::TWO * q0x * qperpx * s121
            - T::TWO * qperpw * qperpy * s001 * theta
            + T::TWO * qperpx * qperpz * s001 * theta
            - T::TWO * q0w * q0x * s011 * theta
            + T::TWO * qperpw * qperpx * s011 * theta
            + T::TWO * qperpy * qperpz * s011 * theta
            + T::TWO * q0x * q0x * s021 * theta
            + T::TWO * q0y * q0y * s021 * theta
            - T::TWO * qperpx * qperpx * s021 * theta
            - T::TWO * qperpy * qperpy * s021 * theta
            + q0z
                * (-(qperpx * s001) - qperpy * s011 + qperpx * s101 + qperpy * s111
                    - T::TWO * q0x * s001 * theta)
            + q0y
                * (-(qperpz * s011)
                    + T::TWO * qperpy * s021
                    + qperpw * (s001 - s101)
                    + qperpz * s111
                    - T::TWO * qperpy * s121
                    + T::TWO * q0w * s001 * theta
                    - T::TWO * q0z * s011 * theta),
        q0w * qperpy * s002 - q0x * qperpz * s002 - q0x * qperpw * s012 - q0w * qperpx * s012
            + T::TWO * q0x * qperpx * s022
            - q0w * qperpy * s102
            + q0x * qperpz * s102
            + q0x * qperpw * s112
            + q0w * qperpx * s112
            - T::TWO * q0x * qperpx * s122
            - T::TWO * qperpw * qperpy * s002 * theta
            + T::TWO * qperpx * qperpz * s002 * theta
            - T::TWO * q0w * q0x * s012 * theta
            + T::TWO * qperpw * qperpx * s012 * theta
            + T::TWO * qperpy * qperpz * s012 * theta
            + T::TWO * q0x * q0x * s022 * theta
            + T::TWO * q0y * q0y * s022 * theta
            - T::TWO * qperpx * qperpx * s022 * theta
            - T::TWO * qperpy * qperpy * s022 * theta
            + q0z
                * (-(qperpx * s002) - qperpy * s012 + qperpx * s102 + qperpy * s112
                    - T::TWO * q0x * s002 * theta)
            + q0y
                * (-(qperpz * s012)
                    + T::TWO * qperpy * s022
                    + qperpw * (s002 - s102)
                    + qperpz * s112
                    - T::TWO * qperpy * s122
                    + T::TWO * q0w * s002 * theta
                    - T::TWO * q0z * s012 * theta),
    );

    c5[2] = DerivativeTerm(
        T::ZERO,
        T::TWO
            * (qperpw * qperpy * s000 - qperpx * qperpz * s000 + q0y * q0z * s010
                - qperpw * qperpx * s010
                - qperpy * qperpz * s010
                - q0y * q0y * s020
                + qperpx * qperpx * s020
                + qperpy * qperpy * s020
                + q0x * q0z * (s000 - s100)
                - qperpw * qperpy * s100
                + qperpx * qperpz * s100
                + q0w * (q0y * (-s000 + s100) + q0x * (s010 - s110))
                - q0y * q0z * s110
                + qperpw * qperpx * s110
                + qperpy * qperpz * s110
                + q0y * q0y * s120
                - qperpx * qperpx * s120
                - qperpy * qperpy * s120
                + q0x * q0x * (-s020 + s120))
            * theta,
        T::TWO
            * (qperpw * qperpy * s001 - qperpx * qperpz * s001 + q0y * q0z * s011
                - qperpw * qperpx * s011
                - qperpy * qperpz * s011
                - q0y * q0y * s021
                + qperpx * qperpx * s021
                + qperpy * qperpy * s021
                + q0x * q0z * (s001 - s101)
                - qperpw * qperpy * s101
                + qperpx * qperpz * s101
                + q0w * (q0y * (-s001 + s101) + q0x * (s011 - s111))
                - q0y * q0z * s111
                + qperpw * qperpx * s111
                + qperpy * qperpz * s111
                + q0y * q0y * s121
                - qperpx * qperpx * s121
                - qperpy * qperpy * s121
                + q0x * q0x * (-s021 + s121))
            * theta,
        T::TWO
            * (qperpw * qperpy * s002 - qperpx * qperpz * s002 + q0y * q0z * s012
                - qperpw * qperpx * s012
                - qperpy * qperpz * s012
                - q0y * q0y * s022
                + qperpx * qperpx * s022
                + qperpy * qperpy * s022
                + q0x * q0z * (s002 - s102)
                - qperpw * qperpy * s102
                + qperpx * qperpz * s102
                + q0w * (q0y * (-s002 + s102) + q0x * (s012 - s112))
                - q0y * q0z * s112
                + qperpw * qperpx * s112
                + qperpy * qperpz * s112
                + q0y * q0y * s122
                - qperpx * qperpx * s122
                - qperpy * qperpy * s122
                + q0x * q0x * (-s022 + s122))
            * theta,
    );

    [c1, c2, c3, c4, c5]
}
