use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::{ConstZero, Matrix4x4, Number, Transform, Vector3};

/// Quaternion (4 component rotation)
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Quaternion<F: Number> {
    pub vec: Vector3<F>,
    pub w: F,
}

impl<F: Number> Default for Quaternion<F> {
    fn default() -> Self {
        Self::ZERO
    }
}

impl<F: Number> ConstZero for Quaternion<F> {
    const ZERO: Self = Self {
        vec: Vector3::ZERO,
        w: F::ZERO,
    };
}

impl<F: Number> Quaternion<F> {
    /// Create a quaternion from its four components (x,y,z,w)
    pub fn new(data: &[F; 4]) -> Self {
        Self {
            vec: Vector3::new(data[0], data[1], data[2]),
            w: data[3],
        }
    }

    /// Inner product (dot product) of two quaternions
    pub fn dot(&self, rhs: &Self) -> F {
        self.vec.dot(rhs.vec) + self.w * rhs.w
    }

    /// Unit length quaternion with the same direction
    pub fn normalise(&self) -> Self {
        *self / self.dot(self).sqrt()
    }

    /// Convert a quaternion rotation into a similar matrix rotation
    pub fn to_transform(&self) -> Transform<F> {
        let xx = self.vec.x * self.vec.x;
        let yy = self.vec.y * self.vec.y;
        let zz = self.vec.z * self.vec.z;
        let xy = self.vec.x * self.vec.y;
        let xz = self.vec.x * self.vec.z;
        let yz = self.vec.y * self.vec.z;
        let ws = self.vec * self.w;

        let mat = Matrix4x4::from_array(&[
            F::ONE - F::TWO * (yy + zz),
            F::TWO * (xy + ws.z),
            F::TWO * (xz - ws.y),
            F::ZERO,
            F::TWO * (xy - ws.z),
            F::ONE - F::TWO * (xx + zz),
            F::TWO * (yz + ws.x),
            F::ZERO,
            F::TWO * (xz + ws.y),
            F::TWO * (yz - ws.x),
            F::ONE - F::TWO * (xx + yy),
            F::ZERO,
            F::ZERO,
            F::ZERO,
            F::ZERO,
            F::ONE,
        ]);

        Transform::new(&mat.transpose(), &mat)
    }

    /// Spherical linear interpolation between two quaternions
    /// Uses minimal rotation amount great circles
    pub fn slerp(&self, t: F, other: &Self) -> Self {
        let cos = self.dot(other);

        if cos > (F::ONE - F::LARGE_EPSILON) {
            (*self * (F::ONE - t) + *other * t).normalise()
        } else {
            let theta = cos.clamp(-F::ONE, F::ONE).acos();
            let theta_t = theta * t;
            let q_perpendicular = (*other - *self * cos).normalise();

            *self * theta_t.cos() + q_perpendicular * theta_t.sin()
        }
    }
}

impl<F: Number> Transform<F> {
    /// Convert a transformation matrix into a quaternion representing its
    /// rotation component.  See [`here`] for where the implementation was taken from
    ///
    /// [`here`]: https://github.com/mmp/pbrt-v3/blob/aaa552a4b9cbf9dccb71450f47b268e0ed6370e2/src/core/quaternion.cpp#L61
    pub fn to_quaternion(&self) -> Quaternion<F> {
        let m = self.mat();
        let trace = m[0][0] + m[1][1] + m[2][2];

        if trace > F::ZERO {
            // Compute w from matrix trace, then xyz
            // 4w^2 = m[0][0] + m[1][1] + m[2][2] + m[3][3] (but m[3][3] == 1)
            let s = (trace + F::ONE).sqrt();
            let w = s / F::TWO;
            let s = F::HALF / s;
            let vec = Vector3::new(
                (m[2][1] - m[1][2]) * s,
                (m[0][2] - m[2][0]) * s,
                (m[1][0] - m[0][1]) * s,
            );

            Quaternion { vec, w }
        } else {
            // Compute largest of x, y, or z, then remaining components
            let nxt = [1, 2, 0];
            let mut q = [F::ZERO; 3];
            let mut i = 0;
            if m[1][1] > m[0][0] {
                i = 1;
            }
            if m[2][2] > m[i][i] {
                i = 2;
            }
            let j = nxt[i];
            let k = nxt[j];
            let mut s = ((m[i][i] - (m[j][j] + m[k][k])) + F::ONE).sqrt();
            q[i] = s * F::HALF;
            if s != F::ZERO {
                s = F::HALF / s;
            }

            let w = (m[k][j] - m[j][k]) * s;
            q[j] = (m[j][i] + m[i][j]) * s;
            q[k] = (m[k][i] + m[i][k]) * s;
            let vec = Vector3::from_array(q);

            Quaternion { vec, w }
        }
    }
}

impl<F: Number> Add<Quaternion<F>> for Quaternion<F> {
    type Output = Quaternion<F>;

    fn add(self, rhs: Quaternion<F>) -> Self::Output {
        Self {
            vec: self.vec + rhs.vec,
            w: self.w + rhs.w,
        }
    }
}

impl<F: Number> AddAssign<Quaternion<F>> for Quaternion<F> {
    fn add_assign(&mut self, rhs: Quaternion<F>) {
        *self = *self + rhs
    }
}

impl<F: Number> Sub<Quaternion<F>> for Quaternion<F> {
    type Output = Quaternion<F>;

    fn sub(self, rhs: Quaternion<F>) -> Self::Output {
        Self {
            vec: self.vec - rhs.vec,
            w: self.w - rhs.w,
        }
    }
}

impl<F: Number> SubAssign<Quaternion<F>> for Quaternion<F> {
    fn sub_assign(&mut self, rhs: Quaternion<F>) {
        *self = *self - rhs
    }
}

impl<F: Number> Neg for Quaternion<F> {
    type Output = Quaternion<F>;

    fn neg(self) -> Self::Output {
        Self {
            vec: -self.vec,
            w: -self.w,
        }
    }
}

impl<F: Number> Mul<F> for Quaternion<F> {
    type Output = Quaternion<F>;

    fn mul(self, rhs: F) -> Self::Output {
        Self {
            vec: self.vec * rhs,
            w: self.w * rhs,
        }
    }
}

impl<F: Number> MulAssign<F> for Quaternion<F> {
    fn mul_assign(&mut self, rhs: F) {
        *self = *self * rhs
    }
}

impl<F: Number> Div<F> for Quaternion<F> {
    type Output = Quaternion<F>;

    fn div(self, rhs: F) -> Self::Output {
        Self {
            vec: self.vec / rhs,
            w: self.w / rhs,
        }
    }
}

impl<F: Number> DivAssign<F> for Quaternion<F> {
    fn div_assign(&mut self, rhs: F) {
        *self = *self / rhs
    }
}
