use std::{
    marker::PhantomData,
    ops::{Add, AddAssign, Index, Mul, MulAssign, Sub, SubAssign},
};

use crate::geometry::{number::Number, Float, Vector2, Vector3};

/// Two dimensional cartesian coordinate
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Point2<T: Number> {
    pub x: T,
    pub y: T,
    _remove_constructors: PhantomData<()>,
}

pub type Point2f = Point2<Float>;
pub type Point2i = Point2<i32>;

impl<T: Number> Point2<T> {
    /// A point at (0, 0)
    pub const ZERO: Self = Self {
        x: T::ZERO,
        y: T::ZERO,
        _remove_constructors: PhantomData,
    };

    /// Create a new point at the given location
    #[inline]
    pub fn new(x: T, y: T) -> Self {
        debug_assert!(!x.is_nan());
        debug_assert!(!y.is_nan());

        Self {
            x,
            y,
            _remove_constructors: PhantomData,
        }
    }

    /// Array of all coordinates of the point
    pub fn to_array(&self) -> [T; 2] {
        [self.x, self.y]
    }

    /// Are any of the coordinates of this point NaN
    pub fn has_nan(&self) -> bool {
        self.x.is_nan() || self.y.is_nan()
    }

    /// Cast the coordinates of a point to another numeric type
    pub fn cast<U: Number>(&self) -> Point2<U> {
        Point2::new(U::cast(self.x), U::cast(self.y))
    }

    /// Convert the point to a vector with the same values
    pub fn to_vec(&self) -> Vector2<T> {
        Vector2::new(self.x, self.y)
    }

    /// Square of the distance between two points
    pub fn distance_square(&self, other: Point2<T>) -> T {
        (*self - other).length_square()
    }

    /// Distance between two points
    pub fn distance(&self, other: Point2<T>) -> T {
        (*self - other).length()
    }

    /// Linear interpolation between two points
    pub fn lerp(&self, other: Point2<T>, t: T) -> Self {
        other * t + *self * (T::ONE - t)
    }

    /// Component-wise minimum of two points
    pub fn min(&self, rhs: Self) -> Self {
        Self::new(
            std::cmp::min_by(self.x, rhs.x, T::order),
            std::cmp::min_by(self.y, rhs.y, T::order),
        )
    }

    /// Component-wise maximum of two points
    pub fn max(&self, rhs: Self) -> Self {
        Self::new(
            std::cmp::max_by(self.x, rhs.x, T::order),
            std::cmp::max_by(self.y, rhs.y, T::order),
        )
    }

    /// Component-wise floor of a point
    pub fn floor(&self) -> Self {
        Self::new(self.x.floor(), self.y.floor())
    }

    /// Component-wise ceiling of a point
    pub fn ceil(&self) -> Self {
        Self::new(self.x.ceil(), self.y.ceil())
    }

    /// Component-wise absolute value of a point
    pub fn abs(&self) -> Self {
        Self::new(self.x.abs(), self.y.abs())
    }

    /// Re-order a point given indices of the where to move each component to
    pub fn permute(&self, idx: [usize; 3]) -> Self {
        Self::new(self[idx[0]], self[idx[1]])
    }
}

impl<T: Number> Index<usize> for Point2<T> {
    type Output = T;

    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        [&self.x, &self.y][index]
    }
}

impl<T: Number> Add for Point2<T> {
    type Output = Self;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Point2::new(self.x + rhs.x, self.y + rhs.y)
    }
}

impl<T: Number> AddAssign for Point2<T> {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl<T: Number> Add<Vector3<T>> for Point2<T> {
    type Output = Self;

    fn add(self, rhs: Vector3<T>) -> Self::Output {
        Point2::new(self.x + rhs.x, self.y + rhs.y)
    }
}

impl<T: Number> AddAssign<Vector3<T>> for Point2<T> {
    fn add_assign(&mut self, rhs: Vector3<T>) {
        *self = *self + rhs
    }
}

impl<T: Number> Sub for Point2<T> {
    type Output = Vector2<T>;

    fn sub(self, rhs: Self) -> Self::Output {
        Vector2::new(self.x - rhs.x, self.y - rhs.y)
    }
}

impl<T: Number> Sub<Vector3<T>> for Point2<T> {
    type Output = Self;

    fn sub(self, rhs: Vector3<T>) -> Self::Output {
        Point2::new(self.x - rhs.x, self.y - rhs.y)
    }
}

impl<T: Number> SubAssign<Vector3<T>> for Point2<T> {
    fn sub_assign(&mut self, rhs: Vector3<T>) {
        *self = *self - rhs
    }
}

impl<T: Number> Mul<T> for Point2<T> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: T) -> Self::Output {
        Point2::new(self.x * rhs, self.y * rhs)
    }
}

impl<T: Number> MulAssign<T> for Point2<T> {
    #[inline]
    fn mul_assign(&mut self, rhs: T) {
        *self = *self * rhs;
    }
}

macro_rules! ForNumbers {
    ($ty:ty) => {
        impl Mul<Point2<$ty>> for $ty {
            type Output = Point2<$ty>;

            #[inline]
            fn mul(self, rhs: Point2<$ty>) -> Self::Output {
                rhs * self
            }
        }
    };
    ($type:ty, $($other:ty),+ $(,)?) => {
        ForNumbers!($type); $(ForNumbers!($other);)+
    }
}

ForNumbers!(f32, f64, i32);
