use std::{
    marker::PhantomData,
    ops::{Add, AddAssign, Div, DivAssign, Index, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::{number::Number, ConstZero, Float, Normal3};

/// Three component numeric vector
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Vector3<T: Number> {
    pub x: T,
    pub y: T,
    pub z: T,
    _remove_constructors: PhantomData<()>,
}

pub type Vector3f = Vector3<Float>;
pub type Vector3i = Vector3<i32>;

impl<T: Number> Default for Vector3<T> {
    fn default() -> Self {
        Self::ZERO
    }
}

impl<T: Number> ConstZero for Vector3<T> {
    const ZERO: Self = Self {
        x: T::ZERO,
        y: T::ZERO,
        z: T::ZERO,
        _remove_constructors: PhantomData,
    };
}

impl<T: Number> Vector3<T> {
    /// Vector pointing along the X axis
    pub const X: Self = Self {
        x: T::ONE,
        y: T::ZERO,
        z: T::ZERO,
        _remove_constructors: PhantomData,
    };

    /// Vector pointing along the X axis
    pub const Y: Self = Self {
        x: T::ZERO,
        y: T::ONE,
        z: T::ZERO,
        _remove_constructors: PhantomData,
    };

    /// Vector pointing along the X axis
    pub const Z: Self = Self {
        x: T::ZERO,
        y: T::ZERO,
        z: T::ONE,
        _remove_constructors: PhantomData,
    };

    /// Create a new vector with the given components
    #[inline]
    pub fn new(x: T, y: T, z: T) -> Self {
        debug_assert!(!x.is_nan());
        debug_assert!(!y.is_nan());
        debug_assert!(!z.is_nan());

        Self {
            x,
            y,
            z,
            _remove_constructors: PhantomData,
        }
    }

    /// Create a vector with all components being equal
    pub fn splat(val: T) -> Self {
        Self::new(val, val, val)
    }

    /// Cast the location of a vector to another numeric type
    pub fn cast<U: Number>(&self) -> Vector3<U> {
        Vector3::new(U::cast(self.x), U::cast(self.y), U::cast(self.z))
    }

    /// Convert a vector3 into the equivalent normal3
    pub fn to_normal(&self) -> Normal3<T> {
        Normal3::new(self.x, self.y, self.z)
    }

    /// Array of all components of the vector
    pub fn to_array(&self) -> [T; 3] {
        [self.x, self.y, self.z]
    }

    /// Construct a vector from its (x,y,z) components
    pub fn from_array(data: [T; 3]) -> Self {
        Self::new(data[0], data[1], data[2])
    }

    /// Are any of the components of this vector NaN
    pub fn has_nan(&self) -> bool {
        self.x.is_nan() || self.y.is_nan() || self.z.is_nan()
    }

    /// Component-wise absolute value of the vector
    pub fn abs(&self) -> Self {
        Vector3::new(self.x.abs(), self.y.abs(), self.z.abs())
    }

    /// Calculate the dot product of two vectors
    pub fn dot(&self, rhs: Self) -> T {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z
    }

    /// Calculate the absolute value of the dot product of two vectors
    pub fn absdot(&self, rhs: Self) -> T {
        self.dot(rhs).abs()
    }

    /// Calculate the cross product of two vectors
    pub fn cross(&self, rhs: Self) -> Self {
        Vector3::new(
            (self.y * rhs.z) - (self.z * rhs.y),
            (self.z * rhs.x) - (self.x * rhs.z),
            (self.x * rhs.y) - (self.y * rhs.x),
        )
    }

    /// Square of the length of the vector
    pub fn length_square(&self) -> T {
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    /// Length of the vector
    pub fn length(&self) -> T {
        self.length_square().sqrt()
    }

    /// Unit vector/Normalise the vector
    pub fn normalise(&self) -> Self {
        *self / self.length()
    }

    /// Value of the smallest component of the vector
    pub fn min_component(&self) -> T {
        self.to_array().into_iter().min_by(T::order).unwrap()
    }

    /// Value of the largest component of the vector
    pub fn max_component(&self) -> T {
        self.to_array().into_iter().max_by(T::order).unwrap()
    }

    /// Index of the largest component of the vector
    pub fn max_dimension(&self) -> usize {
        self.to_array()
            .into_iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| T::order(a, b))
            .unwrap()
            .0
    }

    /// Component-wise minimum of two vectors
    pub fn min(&self, rhs: Self) -> Self {
        Vector3::new(
            std::cmp::min_by(self.x, rhs.x, T::order),
            std::cmp::min_by(self.y, rhs.y, T::order),
            std::cmp::min_by(self.z, rhs.z, T::order),
        )
    }

    /// Component-wise maximum of two vectors
    pub fn max(&self, rhs: Self) -> Self {
        Vector3::new(
            std::cmp::max_by(self.x, rhs.x, T::order),
            std::cmp::max_by(self.y, rhs.y, T::order),
            std::cmp::max_by(self.z, rhs.z, T::order),
        )
    }

    /// Re-order a vector given indices of the where to move each component to
    pub fn permute(&self, idx: [usize; 3]) -> Self {
        Vector3::new(self[idx[0]], self[idx[1]], self[idx[2]])
    }

    /// Construct a coordinate system from a single axis
    pub fn coordinate_system(&self) -> (Self, Self) {
        let v2 = if self.x.abs() > self.y.abs() {
            Vector3::new(-self.z, T::ZERO, self.x) / (self.x * self.x + self.z * self.z).sqrt()
        } else {
            Vector3::new(T::ZERO, self.z, -self.y) / (self.y * self.y + self.z * self.z).sqrt()
        };

        (v2, self.cross(v2))
    }
}

impl<T: Number> Index<usize> for Vector3<T> {
    type Output = T;

    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        [&self.x, &self.y, &self.z][index]
    }
}

impl<T: Number> Add for Vector3<T> {
    type Output = Vector3<T>;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Vector3::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl<T: Number> AddAssign for Vector3<T> {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl<T: Number> Sub for Vector3<T> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        Vector3::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl<T: Number> SubAssign for Vector3<T> {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl<T: Number> Mul<T> for Vector3<T> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: T) -> Self::Output {
        Vector3::new(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}

impl<T: Number> MulAssign<T> for Vector3<T> {
    #[inline]
    fn mul_assign(&mut self, rhs: T) {
        *self = *self * rhs;
    }
}

impl<T: Number> Div<T> for Vector3<T> {
    type Output = Self;

    #[inline]
    fn div(self, rhs: T) -> Self::Output {
        debug_assert_ne!(rhs, T::ZERO);
        let inv = T::ONE / rhs;
        Vector3::new(self.x * inv, self.y * inv, self.z * inv)
    }
}

impl<T: Number> DivAssign<T> for Vector3<T> {
    #[inline]
    fn div_assign(&mut self, rhs: T) {
        *self = *self / rhs;
    }
}

impl<T: Number> Neg for Vector3<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Self::new(-self.x, -self.y, -self.z)
    }
}

macro_rules! ForNumbers {
    ($ty:ty) => {
        impl Mul<Vector3<$ty>> for $ty {
            type Output = Vector3<$ty>;

            #[inline]
            fn mul(self, rhs: Vector3<$ty>) -> Self::Output {
                rhs * self
            }
        }
    };
    ($type:ty, $($other:ty),+ $(,)?) => {
        ForNumbers!($type); $(ForNumbers!($other);)+
    }
}

ForNumbers!(f32, f64, i32);
