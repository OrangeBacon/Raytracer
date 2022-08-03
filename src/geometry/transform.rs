use crate::geometry::{
    number::Number, Bounds3, Float, Matrix4x4, Normal3, Point3, Point3f, Ray, Vector3, Vector3f,
};

/// A 3d transformation matrix
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Transform {
    mat: Matrix4x4,
    inv: Matrix4x4,
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
            [1.0, 0.0, 1.0, delta.z],
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
            [1.0, 0.0, delta.z, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]);

        let inv = Matrix4x4::new(&[
            [1.0 / delta.x, 0.0, 0.0, 0.0],
            [0.0, 1.0 / delta.y, 0.0, 0.0],
            [1.0, 0.0, 1.0 / delta.z, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]);

        Transform::new(&mat, &inv)
    }

    /// Does this transformation matrix have a scale component?
    pub fn has_scale(&self) -> bool {
        let a = self.apply(Vector3f::X).length_square();
        let b = self.apply(Vector3f::Y).length_square();
        let c = self.apply(Vector3f::Z).length_square();

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
            [-sin, 1.0, cos, 0.0],
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
                a.x * a.y * (1.0 - cos) - a.z * sin,
                a.y * a.y + (1.0 - a.y * a.y) * cos,
                a.y * a.z * (1.0 - cos) + a.y * sin,
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

    /// Apply a transformation matrix to any applicable type
    pub fn apply<T: Applicable>(&self, val: T) -> <T as Applicable>::Ret {
        val.apply(self)
    }
}

/// Allows a type to have a transformation matrix be applied to it
pub trait Applicable {
    /// The output of the transformation
    type Ret;

    /// Apply the transformation
    fn apply(&self, transform: &Transform) -> Self::Ret;
}

impl<T: Number> Applicable for Point3<T> {
    type Ret = Point3<T>;

    fn apply(&self, transform: &Transform) -> Self::Ret {
        let m = &transform.mat;
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

impl<T: Number> Applicable for Vector3<T> {
    type Ret = Vector3<T>;

    fn apply(&self, transform: &Transform) -> Self::Ret {
        let m = &transform.mat;
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

impl<T: Number> Applicable for Normal3<T> {
    type Ret = Normal3<T>;

    fn apply(&self, transform: &Transform) -> Self::Ret {
        let m = &transform.inv;
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

impl Applicable for Ray {
    type Ret = Ray;

    fn apply(&self, transform: &Transform) -> Self::Ret {
        let origin = transform.apply(self.origin);
        let direction = transform.apply(self.direction);

        Ray {
            origin,
            direction,
            ..*self
        }
    }
}

impl<T: Number> Applicable for Bounds3<T> {
    type Ret = Bounds3<T>;

    fn apply(&self, transform: &Transform) -> Self::Ret {
        let m = |x| transform.apply(x);

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
