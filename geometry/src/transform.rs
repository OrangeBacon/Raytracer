use std::ops::{Mul, MulAssign};

use crate::{
    number::Number, ray::RayDifferentials, Bounds3, Float, Matrix4x4, Normal3, Point3, Point3f,
    Ray, RayDifferential, Vector3, Vector3f,
};

/// A 3d transformation matrix representing an affine transformation
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Transform {
    mat: Matrix4x4,
    inv: Matrix4x4,
}

impl Default for Transform {
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl Transform {
    /// Transform that does not move an object at all
    pub const IDENTITY: Self = Self {
        mat: Matrix4x4::IDENTITY,
        inv: Matrix4x4::IDENTITY,
    };

    /// Construct a transform with the given inverse matrix.  Assumes the inverse
    /// is correct.
    pub fn new(mat: &Matrix4x4, inv: &Matrix4x4) -> Self {
        Self {
            mat: *mat,
            inv: *inv,
        }
    }

    /// Create a transformation from its matrix data
    /// returns None if the matrix is singular
    pub fn from_array(data: &[[Float; 4]; 4]) -> Option<Self> {
        let mat = Matrix4x4::new(data);
        Some(Self {
            mat,
            inv: mat.inverse()?,
        })
    }

    /// Create a transformation from its matrix data
    /// returns None if the matrix is singular
    pub fn from_mat(data: &Matrix4x4) -> Option<Self> {
        Some(Self {
            mat: *data,
            inv: data.inverse()?,
        })
    }

    // note: using getters here so there are no setters to break inv = mat.inv()
    // if could expose read-only fields that would be better

    /// Get the transformation matrix for this transform
    pub fn mat(&self) -> &Matrix4x4 {
        &self.mat
    }

    /// Get the inverse of the transformation matrix of this transform
    pub fn inv(&self) -> &Matrix4x4 {
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
    pub fn translation(delta: Vector3f) -> Transform {
        let mat = Matrix4x4::new(&[
            [1.0, 0.0, 0.0, delta.x],
            [0.0, 1.0, 0.0, delta.y],
            [0.0, 0.0, 1.0, delta.z],
            [0.0, 0.0, 0.0, 1.0],
        ]);

        let inv = Matrix4x4::new(&[
            [1.0, 0.0, 0.0, -delta.x],
            [0.0, 1.0, 0.0, -delta.y],
            [0.0, 0.0, 1.0, -delta.z],
            [0.0, 0.0, 0.0, 1.0],
        ]);

        Transform::new(&mat, &inv)
    }

    /// Construct a scale matrix
    pub fn scale(delta: Vector3f) -> Transform {
        let mat = Matrix4x4::new(&[
            [delta.x, 0.0, 0.0, 0.0],
            [0.0, delta.y, 0.0, 0.0],
            [0.0, 0.0, delta.z, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]);

        let inv = Matrix4x4::new(&[
            [1.0 / delta.x, 0.0, 0.0, 0.0],
            [0.0, 1.0 / delta.y, 0.0, 0.0],
            [0.0, 0.0, 1.0 / delta.z, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]);

        Transform::new(&mat, &inv)
    }

    /// Does this transformation matrix have a scale component?
    pub fn has_scale(&self) -> bool {
        let a = (Vector3f::X * *self).length_square();
        let b = (Vector3f::Y * *self).length_square();
        let c = (Vector3f::Z * *self).length_square();

        let not_one = |x| x < 0.999 || x > 1.001;

        not_one(a) || not_one(b) || not_one(c)
    }

    /// Rotation matrix around the X axis, angle in degrees
    pub fn rotate_x(angle: Float) -> Transform {
        let sin = angle.to_radians().sin();
        let cos = angle.to_radians().cos();

        let mat = Matrix4x4::new(&[
            [1.0, 0.0, 0.0, 0.0],
            [0.0, cos, -sin, 0.0],
            [0.0, sin, cos, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]);

        Self {
            mat,
            inv: mat.transpose(),
        }
    }

    /// Rotation matrix around the Y axis, angle in degrees
    pub fn rotate_y(angle: Float) -> Transform {
        let sin = angle.to_radians().sin();
        let cos = angle.to_radians().cos();

        let mat = Matrix4x4::new(&[
            [cos, 0.0, sin, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [-sin, 0.0, cos, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]);

        Self {
            mat,
            inv: mat.transpose(),
        }
    }

    /// Rotation matrix around the Z axis, angle in degrees
    pub fn rotate_z(angle: Float) -> Transform {
        let sin = angle.to_radians().sin();
        let cos = angle.to_radians().cos();

        let mat = Matrix4x4::new(&[
            [cos, -sin, 0.0, 0.0],
            [sin, cos, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]);

        Self {
            mat,
            inv: mat.transpose(),
        }
    }

    /// Calculate a rotation around an arbitrary axis, angle in degrees
    pub fn rotate(axis: Vector3f, angle: Float) -> Transform {
        let a = axis.normalise();
        let sin = angle.to_radians().sin();
        let cos = angle.to_radians().cos();

        let mat = Matrix4x4::new(&[
            [
                a.x * a.x + (1.0 - a.x * a.x) * cos,
                a.x * a.y * (1.0 - cos) - a.z * sin,
                a.x * a.z * (1.0 - cos) + a.y * sin,
                0.0,
            ],
            [
                a.x * a.y * (1.0 - cos) + a.z * sin,
                a.y * a.y + (1.0 - a.y * a.y) * cos,
                a.y * a.z * (1.0 - cos) - a.x * sin,
                0.0,
            ],
            [
                a.x * a.z * (1.0 - cos) - a.y * sin,
                a.y * a.z * (1.0 - cos) + a.x * sin,
                a.z * a.z + (1.0 - a.z * a.z) * cos,
                0.0,
            ],
            [0.0, 0.0, 0.0, 1.0],
        ]);

        Self {
            mat,
            inv: mat.transpose(),
        }
    }

    /// A transformation matrix representing looking at the given point from pos
    /// where the upwards direction is the provided vector
    /// Returns none if the transform creates a singular matrix
    pub fn look_at(pos: Point3f, look: Point3f, up: Vector3f) -> Option<Transform> {
        let dir = (look - pos).normalise();
        let right = up.normalise().cross(dir).normalise();
        let new_up = dir.cross(right);

        let mat = Matrix4x4::new(&[
            [right.x, right.y, right.z, 0.0],
            [new_up.x, new_up.y, new_up.z, 0.0],
            [dir.x, dir.y, dir.z, 0.0],
            [pos.x, pos.y, pos.z, 1.0],
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

        det < 0.0
    }
}

impl Mul<Transform> for Transform {
    type Output = Transform;

    fn mul(self, rhs: Transform) -> Self::Output {
        Transform {
            mat: self.mat * rhs.mat,
            inv: rhs.inv * self.inv,
        }
    }
}

impl MulAssign<Transform> for Transform {
    fn mul_assign(&mut self, rhs: Transform) {
        *self = *self * rhs
    }
}

impl<T: Number> Mul<Transform> for Point3<T> {
    type Output = Point3<T>;

    fn mul(self, rhs: Transform) -> Self::Output {
        let m = &rhs.mat;
        let [x, y, z] = [
            Float::cast(self.x),
            Float::cast(self.y),
            Float::cast(self.z),
        ];
        let xp = T::cast(m[0][0] * x + m[0][1] * y + m[0][2] * z + m[0][3]);
        let yp = T::cast(m[1][0] * x + m[1][1] * y + m[1][2] * z + m[1][3]);
        let zp = T::cast(m[2][0] * x + m[2][1] * y + m[2][2] * z + m[2][3]);
        let wp = T::cast(m[3][0] * x + m[3][1] * y + m[3][2] * z + m[3][3]);

        if wp == T::ONE {
            Point3::new(xp, yp, zp)
        } else {
            Point3::new(xp, yp, zp) / wp
        }
    }
}

impl<T: Number> MulAssign<Transform> for Point3<T> {
    fn mul_assign(&mut self, rhs: Transform) {
        *self = *self * rhs
    }
}

impl<T: Number> Mul<Transform> for Vector3<T> {
    type Output = Vector3<T>;

    fn mul(self, rhs: Transform) -> Self::Output {
        let m = &rhs.mat;
        let [x, y, z] = [
            Float::cast(self.x),
            Float::cast(self.y),
            Float::cast(self.z),
        ];
        let xp = T::cast(m[0][0] * x + m[0][1] * y + m[0][2] * z);
        let yp = T::cast(m[1][0] * x + m[1][1] * y + m[1][2] * z);
        let zp = T::cast(m[2][0] * x + m[2][1] * y + m[2][2] * z);

        Vector3::new(xp, yp, zp)
    }
}

impl<T: Number> MulAssign<Transform> for Vector3<T> {
    fn mul_assign(&mut self, rhs: Transform) {
        *self = *self * rhs
    }
}

impl<T: Number> Mul<Transform> for Normal3<T> {
    type Output = Normal3<T>;

    fn mul(self, rhs: Transform) -> Self::Output {
        let m = &rhs.inv;
        let [x, y, z] = [
            Float::cast(self.x),
            Float::cast(self.y),
            Float::cast(self.z),
        ];
        let xp = T::cast(m[0][0] * x + m[1][0] * y + m[2][0] * z);
        let yp = T::cast(m[0][1] * x + m[1][1] * y + m[2][1] * z);
        let zp = T::cast(m[0][2] * x + m[1][2] * y + m[2][2] * z);

        Normal3::new(xp, yp, zp)
    }
}

impl<T: Number> MulAssign<Transform> for Normal3<T> {
    fn mul_assign(&mut self, rhs: Transform) {
        *self = *self * rhs
    }
}

impl Mul<Transform> for Ray {
    type Output = Ray;

    fn mul(self, rhs: Transform) -> Self::Output {
        let origin = self.origin * rhs;
        let direction = self.direction * rhs;

        Ray {
            origin,
            direction,
            ..self
        }
    }
}

impl MulAssign<Transform> for Ray {
    fn mul_assign(&mut self, rhs: Transform) {
        *self = *self * rhs
    }
}

impl Mul<Transform> for RayDifferential {
    type Output = RayDifferential;

    fn mul(self, rhs: Transform) -> Self::Output {
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

impl<T: Number> Mul<Transform> for Bounds3<T> {
    type Output = Bounds3<T>;

    fn mul(self, rhs: Transform) -> Self::Output {
        let m = |x| x * rhs;

        let ret = Bounds3::at_point(m(Point3::new(self.min.x, self.min.y, self.min.z)));
        let ret = ret.union_point(m(Point3::new(self.max.x, self.min.y, self.min.z)));
        let ret = ret.union_point(m(Point3::new(self.min.x, self.max.y, self.min.z)));
        let ret = ret.union_point(m(Point3::new(self.max.x, self.min.y, self.max.z)));
        let ret = ret.union_point(m(Point3::new(self.min.x, self.max.y, self.max.z)));
        let ret = ret.union_point(m(Point3::new(self.max.x, self.max.y, self.min.z)));
        let ret = ret.union_point(m(Point3::new(self.max.x, self.min.y, self.max.z)));
        let ret = ret.union_point(m(Point3::new(self.max.x, self.max.y, self.max.z)));

        ret
    }
}

impl<T: Number> MulAssign<Transform> for Bounds3<T> {
    fn mul_assign(&mut self, rhs: Transform) {
        *self = *self * rhs
    }
}
