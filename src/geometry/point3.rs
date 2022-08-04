use std::{
    marker::PhantomData,
    ops::{Add, AddAssign, Div, DivAssign, Index, Mul, MulAssign, Sub, SubAssign},
};

use crate::geometry::{number::Number, Float, Point2, Vector3};

/// Three dimensional cartesian coordinate
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Point3<T: Number> {
    pub x: T,
    pub y: T,
    pub z: T,
    _remove_constructors: PhantomData<()>,
}

pub type Point3f = Point3<Float>;
pub type Point3i = Point3<i32>;

impl<T: Number> Default for Point3<T> {
    fn default() -> Self {
        Self::ZERO
    }
}

impl<T: Number> Point3<T> {
    /// A point at (0, 0, 0)
    pub const ZERO: Self = Self {
        x: T::ZERO,
        y: T::ZERO,
        z: T::ZERO,
        _remove_constructors: PhantomData,
    };

    /// Smallest possible point
    pub const MIN: Self = Self {
        x: T::MIN,
        y: T::MIN,
        z: T::MIN,
        _remove_constructors: PhantomData,
    };

    /// Largest possible point
    pub const MAX: Self = Self {
        x: T::MAX,
        y: T::MAX,
        z: T::MAX,
        _remove_constructors: PhantomData,
    };

    /// Create a new point at the given location
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

    /// Array of all coordinates of the point
    pub fn to_array(&self) -> [T; 3] {
        [self.x, self.y, self.z]
    }

    /// Are any of the coordinates of this point NaN
    pub fn has_nan(&self) -> bool {
        self.x.is_nan() || self.y.is_nan() || self.z.is_nan()
    }

    /// Drop the z coordinate of a point
    pub fn as_point2(&self) -> Point2<T> {
        Point2::new(self.x, self.y)
    }

    /// Cast the coordinates of a point to another numeric type
    pub fn cast<U: Number>(&self) -> Point3<U> {
        Point3::new(U::cast(self.x), U::cast(self.y), U::cast(self.z))
    }

    /// Convert the point to a vector with the same values
    pub fn to_vec(&self) -> Vector3<T> {
        Vector3::new(self.x, self.y, self.z)
    }

    /// Square of the distance between two points
    pub fn distance_square(&self, other: Point3<T>) -> T {
        (*self - other).length_square()
    }

    /// Distance between two points
    pub fn distance(&self, other: Point3<T>) -> T {
        (*self - other).length()
    }

    /// Linear interpolation between two points
    pub fn lerp(&self, other: Point3<T>, t: T) -> Self {
        other * t + *self * (T::ONE - t)
    }

    /// Component-wise minimum of two points
    pub fn min(&self, rhs: Self) -> Self {
        Self::new(
            std::cmp::min_by(self.x, rhs.x, T::order),
            std::cmp::min_by(self.y, rhs.y, T::order),
            std::cmp::min_by(self.z, rhs.z, T::order),
        )
    }

    /// Component-wise maximum of two points
    pub fn max(&self, rhs: Self) -> Self {
        Self::new(
            std::cmp::max_by(self.x, rhs.x, T::order),
            std::cmp::max_by(self.y, rhs.y, T::order),
            std::cmp::max_by(self.z, rhs.z, T::order),
        )
    }

    /// Component-wise floor of a point
    pub fn floor(&self) -> Self {
        Self::new(self.x.floor(), self.y.floor(), self.z.floor())
    }

    /// Component-wise ceiling of a point
    pub fn ceil(&self) -> Self {
        Self::new(self.x.ceil(), self.y.ceil(), self.z.ceil())
    }

    /// Component-wise absolute value of a point
    pub fn abs(&self) -> Self {
        Self::new(self.x.abs(), self.y.abs(), self.z.abs())
    }

    /// Re-order a point given indices of the where to move each component to
    pub fn permute(&self, idx: [usize; 3]) -> Self {
        Self::new(self[idx[0]], self[idx[1]], self[idx[2]])
    }
}

impl<T: Number> Index<usize> for Point3<T> {
    type Output = T;

    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        [&self.x, &self.y, &self.z][index]
    }
}

impl<T: Number> Add for Point3<T> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Point3::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl<T: Number> AddAssign for Point3<T> {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl<T: Number> Add<Vector3<T>> for Point3<T> {
    type Output = Self;

    fn add(self, rhs: Vector3<T>) -> Self::Output {
        Point3::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z)
    }
}

impl<T: Number> AddAssign<Vector3<T>> for Point3<T> {
    fn add_assign(&mut self, rhs: Vector3<T>) {
        *self = *self + rhs
    }
}

impl<T: Number> Sub for Point3<T> {
    type Output = Vector3<T>;

    fn sub(self, rhs: Self) -> Self::Output {
        Vector3::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl<T: Number> Sub<Vector3<T>> for Point3<T> {
    type Output = Self;

    fn sub(self, rhs: Vector3<T>) -> Self::Output {
        Point3::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z)
    }
}

impl<T: Number> SubAssign<Vector3<T>> for Point3<T> {
    fn sub_assign(&mut self, rhs: Vector3<T>) {
        *self = *self - rhs
    }
}

impl<T: Number> Mul<T> for Point3<T> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: T) -> Self::Output {
        Point3::new(self.x * rhs, self.y * rhs, self.z * rhs)
    }
}

impl<T: Number> MulAssign<T> for Point3<T> {
    #[inline]
    fn mul_assign(&mut self, rhs: T) {
        *self = *self * rhs;
    }
}

impl<T: Number> Div<T> for Point3<T> {
    type Output = Self;

    #[inline]
    fn div(self, rhs: T) -> Self::Output {
        debug_assert_ne!(rhs, T::ZERO);
        let inv = T::ONE / rhs;
        Point3::new(self.x * inv, self.y * inv, self.z * inv)
    }
}

impl<T: Number> DivAssign<T> for Point3<T> {
    #[inline]
    fn div_assign(&mut self, rhs: T) {
        *self = *self / rhs;
    }
}

macro_rules! ForNumbers {
    ($ty:ty) => {
        impl Mul<Point3<$ty>> for $ty {
            type Output = Point3<$ty>;

            #[inline]
            fn mul(self, rhs: Point3<$ty>) -> Self::Output {
                rhs * self
            }
        }
    };
    ($type:ty, $($other:ty),+ $(,)?) => {
        ForNumbers!($type); $(ForNumbers!($other);)+
    }
}

ForNumbers!(f32, f64, i32);
