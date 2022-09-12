use std::ops::{Mul, MulAssign};

use crate::{
    float::gamma, number::Number, ray::RayDifferentials, Bounds3, Matrix4x4, Normal3, Point3, Ray,
    RayDifferential, Vector3,
};

/// A 3d transformation matrix representing an affine transformation
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Transform<T: Number> {
    mat: Matrix4x4<T>,
    inv: Matrix4x4<T>,
}

/// Apply a transform to any type
pub trait Applicable<T> {
    fn apply(&self, other: T) -> T;
}

/// Apply a transform to any type and get error bounds
pub trait ApplicableError<T, E, R> {
    fn apply(&self, other: T) -> (R, E);
}

impl<T: Number> Default for Transform<T> {
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl<T: Number> Transform<T> {
    /// Transform that does not move an object at all
    pub const IDENTITY: Self = Self {
        mat: Matrix4x4::IDENTITY,
        inv: Matrix4x4::IDENTITY,
    };

    /// Construct a transform with the given inverse matrix.  Assumes the inverse
    /// is correct.
    pub fn new(mat: &Matrix4x4<T>, inv: &Matrix4x4<T>) -> Self {
        Self {
            mat: *mat,
            inv: *inv,
        }
    }

    /// Create a transformation from its matrix data
    /// returns None if the matrix is singular
    pub fn from_array(data: &[[T; 4]; 4]) -> Option<Self> {
        let mat = Matrix4x4::new(data);
        Some(Self {
            mat,
            inv: mat.inverse()?,
        })
    }

    /// Create a transformation from its matrix data
    /// returns None if the matrix is singular
    pub fn from_mat(data: &Matrix4x4<T>) -> Option<Self> {
        Some(Self {
            mat: *data,
            inv: data.inverse()?,
        })
    }

    // note: using getters here so there are no setters to break inv = mat.inv()
    // if could expose read-only fields that would be better

    /// Get the transformation matrix for this transform
    pub fn mat(&self) -> &Matrix4x4<T> {
        &self.mat
    }

    /// Get the inverse of the transformation matrix of this transform
    pub fn inv(&self) -> &Matrix4x4<T> {
        &self.inv
    }

    /// Create the inverse of this transform
    pub fn inverse(&self) -> Self {
        Self {
            mat: self.inv,
            inv: self.mat,
        }
    }

    /// Transpose of the transform
    pub fn transpose(&self) -> Self {
        Self {
            mat: self.mat.transpose(),
            inv: self.inv.transpose(),
        }
    }

    /// Is this transform the identity transform
    pub fn is_identity(&self) -> bool {
        self == &Self::IDENTITY
    }

    /// Construct a translation matrix
    pub fn translation(delta: Vector3<T>) -> Transform<T> {
        let mat = Matrix4x4::new(&[
            [T::ONE, T::ZERO, T::ZERO, delta.x],
            [T::ZERO, T::ONE, T::ZERO, delta.y],
            [T::ZERO, T::ZERO, T::ONE, delta.z],
            [T::ZERO, T::ZERO, T::ZERO, T::ONE],
        ]);

        let inv = Matrix4x4::new(&[
            [T::ONE, T::ZERO, T::ZERO, -delta.x],
            [T::ZERO, T::ONE, T::ZERO, -delta.y],
            [T::ZERO, T::ZERO, T::ONE, -delta.z],
            [T::ZERO, T::ZERO, T::ZERO, T::ONE],
        ]);

        Transform::new(&mat, &inv)
    }

    /// Construct a scale matrix
    pub fn scale(delta: Vector3<T>) -> Self {
        let mat = Matrix4x4::new(&[
            [delta.x, T::ZERO, T::ZERO, T::ZERO],
            [T::ZERO, delta.y, T::ZERO, T::ZERO],
            [T::ZERO, T::ZERO, delta.z, T::ZERO],
            [T::ZERO, T::ZERO, T::ZERO, T::ONE],
        ]);

        let inv = Matrix4x4::new(&[
            [T::ONE / delta.x, T::ZERO, T::ZERO, T::ZERO],
            [T::ZERO, T::ONE / delta.y, T::ZERO, T::ZERO],
            [T::ZERO, T::ZERO, T::ONE / delta.z, T::ZERO],
            [T::ZERO, T::ZERO, T::ZERO, T::ONE],
        ]);

        Transform::new(&mat, &inv)
    }

    /// Does this transformation matrix have a scale component?
    pub fn has_scale(&self) -> bool {
        let a = (Vector3::X * *self).length_square();
        let b = (Vector3::Y * *self).length_square();
        let c = (Vector3::Z * *self).length_square();

        let not_one = |x| !((T::ONE - T::LARGE_EPSILON)..=(T::ONE + T::LARGE_EPSILON)).contains(&x);

        not_one(a) || not_one(b) || not_one(c)
    }

    /// Rotation matrix around the X axis, angle in degrees
    pub fn rotate_x(angle: T) -> Self {
        let sin = angle.to_radians().sin();
        let cos = angle.to_radians().cos();

        let mat = Matrix4x4::new(&[
            [T::ONE, T::ZERO, T::ZERO, T::ZERO],
            [T::ZERO, cos, -sin, T::ZERO],
            [T::ZERO, sin, cos, T::ZERO],
            [T::ZERO, T::ZERO, T::ZERO, T::ONE],
        ]);

        Self {
            mat,
            inv: mat.transpose(),
        }
    }

    /// Rotation matrix around the Y axis, angle in degrees
    pub fn rotate_y(angle: T) -> Self {
        let sin = angle.to_radians().sin();
        let cos = angle.to_radians().cos();

        let mat = Matrix4x4::new(&[
            [cos, T::ZERO, sin, T::ZERO],
            [T::ZERO, T::ONE, T::ZERO, T::ZERO],
            [-sin, T::ZERO, cos, T::ZERO],
            [T::ZERO, T::ZERO, T::ZERO, T::ONE],
        ]);

        Self {
            mat,
            inv: mat.transpose(),
        }
    }

    /// Rotation matrix around the Z axis, angle in degrees
    pub fn rotate_z(angle: T) -> Self {
        let sin = angle.to_radians().sin();
        let cos = angle.to_radians().cos();

        let mat = Matrix4x4::new(&[
            [cos, -sin, T::ZERO, T::ZERO],
            [sin, cos, T::ZERO, T::ZERO],
            [T::ZERO, T::ZERO, T::ONE, T::ZERO],
            [T::ZERO, T::ZERO, T::ZERO, T::ONE],
        ]);

        Self {
            mat,
            inv: mat.transpose(),
        }
    }

    /// Calculate a rotation around an arbitrary axis, angle in degrees
    pub fn rotate(axis: Vector3<T>, angle: T) -> Self {
        let a = axis.normalise();
        let sin = angle.to_radians().sin();
        let cos = angle.to_radians().cos();

        let mat = Matrix4x4::new(&[
            [
                a.x * a.x + (T::ONE - a.x * a.x) * cos,
                a.x * a.y * (T::ONE - cos) - a.z * sin,
                a.x * a.z * (T::ONE - cos) + a.y * sin,
                T::ZERO,
            ],
            [
                a.x * a.y * (T::ONE - cos) + a.z * sin,
                a.y * a.y + (T::ONE - a.y * a.y) * cos,
                a.y * a.z * (T::ONE - cos) - a.x * sin,
                T::ZERO,
            ],
            [
                a.x * a.z * (T::ONE - cos) - a.y * sin,
                a.y * a.z * (T::ONE - cos) + a.x * sin,
                a.z * a.z + (T::ONE - a.z * a.z) * cos,
                T::ZERO,
            ],
            [T::ZERO, T::ZERO, T::ZERO, T::ONE],
        ]);

        Self {
            mat,
            inv: mat.transpose(),
        }
    }

    /// A transformation matrix representing looking at the given point from pos
    /// where the upwards direction is the provided vector
    /// Returns none if the transform creates a singular matrix
    pub fn look_at(pos: Point3<T>, look: Point3<T>, up: Vector3<T>) -> Option<Self> {
        let dir = (look - pos).normalise();
        let right = up.normalise().cross(dir).normalise();
        let new_up = dir.cross(right);

        let mat = Matrix4x4::new(&[
            [right.x, right.y, right.z, T::ZERO],
            [new_up.x, new_up.y, new_up.z, T::ZERO],
            [dir.x, dir.y, dir.z, T::ZERO],
            [pos.x, pos.y, pos.z, T::ONE],
        ]);

        Some(Transform {
            mat: mat.inverse()?,
            inv: mat,
        })
    }

    /// Does this transformation convert left handed to right handed coordinates?
    pub fn swaps_handedness(&self) -> bool {
        let mat = &self.mat;
        let det = mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1])
            - mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0])
            + mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);

        det < T::ZERO
    }

    /// Apply this transformation to an object.  Equivalent to other * self
    pub fn apply<U>(&self, other: U) -> U
    where
        Self: Applicable<U>,
    {
        <Self as Applicable<U>>::apply(self, other)
    }

    /// Apply this transformation to an object.  Returns the result along with its error bounds
    pub fn apply_err<U, E, R>(&self, other: U) -> (R, E)
    where
        Self: ApplicableError<U, E, R>,
    {
        <Self as ApplicableError<U, E, R>>::apply(self, other)
    }
}

impl<T: Number> Mul<Transform<T>> for Transform<T> {
    type Output = Transform<T>;

    fn mul(self, rhs: Transform<T>) -> Self::Output {
        Transform {
            mat: self.mat * rhs.mat,
            inv: rhs.inv * self.inv,
        }
    }
}

impl<T: Number> MulAssign<Transform<T>> for Transform<T> {
    fn mul_assign(&mut self, rhs: Transform<T>) {
        *self = *self * rhs
    }
}

impl<T: Number> Applicable<Transform<T>> for Transform<T> {
    fn apply(&self, other: Transform<T>) -> Transform<T> {
        *self * other
    }
}

impl<T: Number> Mul<Transform<T>> for Point3<T> {
    type Output = Point3<T>;

    fn mul(self, rhs: Transform<T>) -> Self::Output {
        let m = &rhs.mat;
        let [x, y, z] = [self.x, self.y, self.z];
        let xp = m[0][0] * x + m[0][1] * y + m[0][2] * z + m[0][3];
        let yp = m[1][0] * x + m[1][1] * y + m[1][2] * z + m[1][3];
        let zp = m[2][0] * x + m[2][1] * y + m[2][2] * z + m[2][3];
        let wp = m[3][0] * x + m[3][1] * y + m[3][2] * z + m[3][3];

        if wp == T::ONE {
            Point3::new(xp, yp, zp)
        } else {
            Point3::new(xp, yp, zp) / wp
        }
    }
}

impl<T: Number> MulAssign<Transform<T>> for Point3<T> {
    fn mul_assign(&mut self, rhs: Transform<T>) {
        *self = *self * rhs
    }
}

impl<T: Number> Applicable<Point3<T>> for Transform<T> {
    fn apply(&self, other: Point3<T>) -> Point3<T> {
        other * *self
    }
}

impl<T: Number> ApplicableError<Point3<T>, Vector3<T>, Point3<T>> for Transform<T> {
    fn apply(&self, other: Point3<T>) -> (Point3<T>, Vector3<T>) {
        let m = &self.mat;
        let [x, y, z] = [other.x, other.y, other.z];
        let xp = m[0][0] * x + m[0][1] * y + m[0][2] * z + m[0][3];
        let yp = m[1][0] * x + m[1][1] * y + m[1][2] * z + m[1][3];
        let zp = m[2][0] * x + m[2][1] * y + m[2][2] * z + m[2][3];
        let wp = m[3][0] * x + m[3][1] * y + m[3][2] * z + m[3][3];

        let abs_sum = [
            (m[0][0] * x).abs() + (m[0][1] * y).abs() + (m[0][2] * z).abs() + m[0][3].abs(),
            (m[1][0] * x).abs() + (m[1][1] * y).abs() + (m[1][2] * z).abs() + m[1][3].abs(),
            (m[2][0] * x).abs() + (m[2][1] * y).abs() + (m[2][2] * z).abs() + m[2][3].abs(),
        ];

        let err = Vector3::from_array(abs_sum) * gamma::<T>(3);

        debug_assert_ne!(wp, T::ZERO);

        if wp == T::ONE {
            (Point3::new(xp, yp, zp), err)
        } else {
            (Point3::new(xp, yp, zp) / wp, err)
        }
    }
}

impl<T: Number> ApplicableError<(Point3<T>, Vector3<T>), Vector3<T>, Point3<T>> for Transform<T> {
    fn apply(&self, (point, error): (Point3<T>, Vector3<T>)) -> (Point3<T>, Vector3<T>) {
        let m = &self.mat;
        let [x, y, z] = [point.x, point.y, point.z];
        let xp = m[0][0] * x + m[0][1] * y + m[0][2] * z + m[0][3];
        let yp = m[1][0] * x + m[1][1] * y + m[1][2] * z + m[1][3];
        let zp = m[2][0] * x + m[2][1] * y + m[2][2] * z + m[2][3];
        let wp = m[3][0] * x + m[3][1] * y + m[3][2] * z + m[3][3];

        let err_sum = [
            (m[0][0] * error.x).abs()
                + (m[0][1] * error.y).abs()
                + (m[0][2] * error.z).abs()
                + m[0][3].abs(),
            (m[1][0] * error.x).abs()
                + (m[1][1] * error.y).abs()
                + (m[1][2] * error.z).abs()
                + m[1][3].abs(),
            (m[2][0] * error.x).abs()
                + (m[2][1] * error.y).abs()
                + (m[2][2] * error.z).abs()
                + m[2][3].abs(),
        ];
        let abs_sum = [
            (m[0][0] * x).abs() + (m[0][1] * y).abs() + (m[0][2] * z).abs() + m[0][3].abs(),
            (m[1][0] * x).abs() + (m[1][1] * y).abs() + (m[1][2] * z).abs() + m[1][3].abs(),
            (m[2][0] * x).abs() + (m[2][1] * y).abs() + (m[2][2] * z).abs() + m[2][3].abs(),
        ];

        let err = Vector3::from_array(err_sum) * (gamma::<T>(3) + T::ONE)
            + Vector3::from_array(abs_sum) * gamma::<T>(3);

        debug_assert_ne!(wp, T::ZERO);

        if wp == T::ONE {
            (Point3::new(xp, yp, zp), err)
        } else {
            (Point3::new(xp, yp, zp) / wp, err)
        }
    }
}

impl<T: Number> Mul<Transform<T>> for Vector3<T> {
    type Output = Vector3<T>;

    fn mul(self, rhs: Transform<T>) -> Self::Output {
        let m = &rhs.mat;
        let [x, y, z] = [self.x, self.y, self.z];
        let xp = m[0][0] * x + m[0][1] * y + m[0][2] * z;
        let yp = m[1][0] * x + m[1][1] * y + m[1][2] * z;
        let zp = m[2][0] * x + m[2][1] * y + m[2][2] * z;

        Vector3::new(xp, yp, zp)
    }
}

impl<T: Number> MulAssign<Transform<T>> for Vector3<T> {
    fn mul_assign(&mut self, rhs: Transform<T>) {
        *self = *self * rhs
    }
}

impl<T: Number> Applicable<Vector3<T>> for Transform<T> {
    fn apply(&self, other: Vector3<T>) -> Vector3<T> {
        other * *self
    }
}

impl<T: Number> ApplicableError<Vector3<T>, Vector3<T>, Vector3<T>> for Transform<T> {
    fn apply(&self, other: Vector3<T>) -> (Vector3<T>, Vector3<T>) {
        let m = &self.mat;
        let [x, y, z] = [other.x, other.y, other.z];
        let err = [
            (m[0][0] * x).abs() + (m[0][1] * y).abs() + (m[0][2] * z).abs() + m[0][3].abs(),
            (m[1][0] * x).abs() + (m[1][1] * y).abs() + (m[1][2] * z).abs() + m[1][3].abs(),
            (m[2][0] * x).abs() + (m[2][1] * y).abs() + (m[2][2] * z).abs() + m[2][3].abs(),
        ];

        (self.apply(other), Vector3::from_array(err) * gamma::<T>(3))
    }
}

impl<T: Number> ApplicableError<(Vector3<T>, Vector3<T>), Vector3<T>, Vector3<T>> for Transform<T> {
    fn apply(&self, (other, error): (Vector3<T>, Vector3<T>)) -> (Vector3<T>, Vector3<T>) {
        let m = &self.mat;
        let [x, y, z] = [other.x, other.y, other.z];

        let err_sum = [
            (m[0][0] * error.x).abs()
                + (m[0][1] * error.y).abs()
                + (m[0][2] * error.z).abs()
                + m[0][3].abs(),
            (m[1][0] * error.x).abs()
                + (m[1][1] * error.y).abs()
                + (m[1][2] * error.z).abs()
                + m[1][3].abs(),
            (m[2][0] * error.x).abs()
                + (m[2][1] * error.y).abs()
                + (m[2][2] * error.z).abs()
                + m[2][3].abs(),
        ];
        let abs_sum = [
            (m[0][0] * x).abs() + (m[0][1] * y).abs() + (m[0][2] * z).abs() + m[0][3].abs(),
            (m[1][0] * x).abs() + (m[1][1] * y).abs() + (m[1][2] * z).abs() + m[1][3].abs(),
            (m[2][0] * x).abs() + (m[2][1] * y).abs() + (m[2][2] * z).abs() + m[2][3].abs(),
        ];

        let err = Vector3::from_array(err_sum) * (gamma::<T>(3) + T::ONE)
            + Vector3::from_array(abs_sum) * gamma::<T>(3);

        (self.apply(other), err)
    }
}

impl<T: Number> Mul<Transform<T>> for Normal3<T> {
    type Output = Normal3<T>;

    fn mul(self, rhs: Transform<T>) -> Self::Output {
        let m = &rhs.inv;
        let [x, y, z] = [self.x, self.y, self.z];
        let xp = m[0][0] * x + m[1][0] * y + m[2][0] * z;
        let yp = m[0][1] * x + m[1][1] * y + m[2][1] * z;
        let zp = m[0][2] * x + m[1][2] * y + m[2][2] * z;

        Normal3::new(xp, yp, zp)
    }
}

impl<T: Number> MulAssign<Transform<T>> for Normal3<T> {
    fn mul_assign(&mut self, rhs: Transform<T>) {
        *self = *self * rhs
    }
}

impl<T: Number> Applicable<Normal3<T>> for Transform<T> {
    fn apply(&self, other: Normal3<T>) -> Normal3<T> {
        other * *self
    }
}

impl<T, F: Number> Mul<Transform<F>> for Ray<T, F> {
    type Output = Self;

    fn mul(self, rhs: Transform<F>) -> Self::Output {
        let (mut origin, o_err) = rhs.apply_err(self.origin);
        let direction = self.direction * rhs;

        let length_square = direction.length_square();
        let mut t_max = self.t_max;
        if length_square > F::ZERO {
            let dt = direction.abs().dot(o_err) / length_square;
            origin += direction * dt;
            t_max -= dt;
        }

        Ray {
            origin,
            direction,
            t_max,
            ..self
        }
    }
}

impl<T, F: Number> MulAssign<Transform<F>> for Ray<T, F> {
    fn mul_assign(&mut self, rhs: Transform<F>) {
        self.origin *= rhs;
        self.direction *= rhs;
    }
}

impl<T, F: Number> Applicable<Ray<T, F>> for Transform<F> {
    fn apply(&self, other: Ray<T, F>) -> Ray<T, F> {
        other * *self
    }
}

impl<T, F: Number> ApplicableError<Ray<T, F>, (Vector3<F>, Vector3<F>), Ray<T, F>>
    for Transform<F>
{
    fn apply(&self, other: Ray<T, F>) -> (Ray<T, F>, (Vector3<F>, Vector3<F>)) {
        let (mut origin, origin_err) = self.apply_err(other.origin);
        let (dir, dir_err) = self.apply_err(other.direction);
        let len_sq = dir.length_square();
        if len_sq > F::ZERO {
            let dt = dir.abs().dot(origin_err) / len_sq;
            origin += dir * dt;
        }

        (
            Ray {
                origin,
                direction: dir,
                t_max: other.t_max,
                time: other.time,
                material: other.material,
            },
            (origin_err, dir_err),
        )
    }
}

impl<T, F: Number>
    ApplicableError<(Ray<T, F>, Vector3<F>, Vector3<F>), (Vector3<F>, Vector3<F>), Ray<T, F>>
    for Transform<F>
{
    fn apply(
        &self,
        (ray, origin_err, direction_err): (Ray<T, F>, Vector3<F>, Vector3<F>),
    ) -> (Ray<T, F>, (Vector3<F>, Vector3<F>)) {
        let (mut origin, origin_err) = self.apply_err((ray.origin, origin_err));
        let (dir, dir_err) = self.apply_err((ray.direction, direction_err));
        let len_sq = dir.length_square();
        if len_sq > F::ZERO {
            let dt = dir.abs().dot(origin_err) / len_sq;
            origin += dir * dt;
        }

        (
            Ray {
                origin,
                direction: dir,
                t_max: ray.t_max,
                time: ray.time,
                material: ray.material,
            },
            (origin_err, dir_err),
        )
    }
}

impl<F: Number, T> Mul<Transform<F>> for RayDifferential<T, F> {
    type Output = Self;

    fn mul(self, rhs: Transform<F>) -> Self::Output {
        let tr = self.main * rhs;
        RayDifferential {
            main: tr,
            differentials: self.differentials.map(|old| RayDifferentials {
                rx_origin: old.rx_origin * rhs,
                ry_origin: old.ry_origin * rhs,
                rx_direction: old.rx_direction * rhs,
                ry_direction: old.ry_direction * rhs,
            }),
        }
    }
}

impl<F: Number, T> MulAssign<Transform<F>> for RayDifferential<T, F> {
    fn mul_assign(&mut self, rhs: Transform<F>) {
        self.main *= rhs;
        self.differentials = self.differentials.map(|old| RayDifferentials {
            rx_origin: old.rx_origin * rhs,
            ry_origin: old.ry_origin * rhs,
            rx_direction: old.rx_direction * rhs,
            ry_direction: old.ry_direction * rhs,
        })
    }
}

impl<F: Number, T> Applicable<RayDifferential<T, F>> for Transform<F> {
    fn apply(&self, other: RayDifferential<T, F>) -> RayDifferential<T, F> {
        other * *self
    }
}

impl<T: Number> Mul<Transform<T>> for Bounds3<T> {
    type Output = Bounds3<T>;

    fn mul(self, rhs: Transform<T>) -> Self::Output {
        let m = |x| x * rhs;

        let ret = Bounds3::at_point(m(Point3::new(self.min.x, self.min.y, self.min.z)));
        let ret = ret.union_point(m(Point3::new(self.max.x, self.min.y, self.min.z)));
        let ret = ret.union_point(m(Point3::new(self.min.x, self.max.y, self.min.z)));
        let ret = ret.union_point(m(Point3::new(self.max.x, self.min.y, self.max.z)));
        let ret = ret.union_point(m(Point3::new(self.min.x, self.max.y, self.max.z)));
        let ret = ret.union_point(m(Point3::new(self.max.x, self.max.y, self.min.z)));
        let ret = ret.union_point(m(Point3::new(self.max.x, self.min.y, self.max.z)));

        ret.union_point(m(Point3::new(self.max.x, self.max.y, self.max.z)))
    }
}

impl<T: Number> MulAssign<Transform<T>> for Bounds3<T> {
    fn mul_assign(&mut self, rhs: Transform<T>) {
        *self = *self * rhs
    }
}

impl<T: Number> Applicable<Bounds3<T>> for Transform<T> {
    fn apply(&self, other: Bounds3<T>) -> Bounds3<T> {
        other * *self
    }
}
