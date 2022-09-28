use std::{
    marker::PhantomData,
    ops::{Add, AddAssign, Div, DivAssign, Index, Mul, MulAssign, Neg, Sub, SubAssign},
};

use crate::{number::Number, ConstZero, Float, Point2};

/// Two component numeric vector
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Vector2<T: Number> {
    pub x: T,
    pub y: T,
    _remove_constructors: PhantomData<()>,
}

pub type Vector2f = Vector2<Float>;
pub type Vector2i = Vector2<i32>;

impl<T: Number> Default for Vector2<T> {
    fn default() -> Self {
        Self::ZERO
    }
}

impl<T: Number> ConstZero for Vector2<T> {
    const ZERO: Self = Self {
        x: T::ZERO,
        y: T::ZERO,
        _remove_constructors: PhantomData,
    };
}

impl<T: Number> Vector2<T> {
    /// Create a new vector with the given components
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

    /// Create a vector with all components being equal
    pub fn splat(val: T) -> Self {
        Self::new(val, val)
    }

    /// Cast the location of a vector to another numeric type
    pub fn cast<U: Number>(&self) -> Vector2<U> {
        Vector2::new(U::cast(self.x), U::cast(self.y))
    }

    /// Array of all components of the vector
    pub fn to_array(&self) -> [T; 2] {
        [self.x, self.y]
    }

    /// Convert this vector into a single point
    pub fn to_point(&self) -> Point2<T> {
        Point2::new(self.x, self.y)
    }

    /// Are any of the components of this vector NaN
    pub fn has_nan(&self) -> bool {
        self.x.is_nan() || self.y.is_nan()
    }

    /// Component-wise absolute value of the vector
    pub fn abs(&self) -> Self {
        Vector2::new(self.x.abs(), self.y.abs())
    }

    /// Calculate the dot product of two vectors
    pub fn dot(&self, rhs: Self) -> T {
        self.x * rhs.x + self.y * rhs.y
    }

    /// Calculate the absolute value of the dot product of two vectors
    pub fn absdot(&self, rhs: Self) -> T {
        self.dot(rhs).abs()
    }

    /// Square of the length of the vector
    pub fn length_square(&self) -> T {
        self.x * self.x + self.y * self.y
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
        Vector2::new(
            std::cmp::min_by(self.x, rhs.x, T::order),
            std::cmp::min_by(self.y, rhs.y, T::order),
        )
    }

    /// Component-wise maximum of two vectors
    pub fn max(&self, rhs: Self) -> Self {
        Vector2::new(
            std::cmp::max_by(self.x, rhs.x, T::order),
            std::cmp::max_by(self.y, rhs.y, T::order),
        )
    }

    /// Re-order a vector given indices of the where to move each component to
    pub fn permute(&self, idx: [usize; 2]) -> Self {
        Vector2::new(self[idx[0]], self[idx[1]])
    }
}

impl<T: Number> Index<usize> for Vector2<T> {
    type Output = T;

    #[inline]
    fn index(&self, index: usize) -> &Self::Output {
        [&self.x, &self.y][index]
    }
}

impl<T: Number> Add for Vector2<T> {
    type Output = Vector2<T>;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Vector2::new(self.x + rhs.x, self.y + rhs.y)
    }
}

impl<T: Number> AddAssign for Vector2<T> {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl<T: Number> Sub for Vector2<T> {
    type Output = Self;

    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        Vector2::new(self.x - rhs.x, self.y - rhs.y)
    }
}

impl<T: Number> SubAssign for Vector2<T> {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl<T: Number> Mul<T> for Vector2<T> {
    type Output = Self;

    #[inline]
    fn mul(self, rhs: T) -> Self::Output {
        Vector2::new(self.x * rhs, self.y * rhs)
    }
}

impl<T: Number> MulAssign<T> for Vector2<T> {
    #[inline]
    fn mul_assign(&mut self, rhs: T) {
        *self = *self * rhs;
    }
}

impl<T: Number> Div<T> for Vector2<T> {
    type Output = Self;

    #[inline]
    fn div(self, rhs: T) -> Self::Output {
        debug_assert_ne!(rhs, T::ZERO);
        let inv = T::ONE / rhs;
        Vector2::new(self.x * inv, self.y * inv)
    }
}

impl<T: Number> DivAssign<T> for Vector2<T> {
    #[inline]
    fn div_assign(&mut self, rhs: T) {
        *self = *self / rhs;
    }
}

impl<T: Number> Neg for Vector2<T> {
    type Output = Self;

    #[inline]
    fn neg(self) -> Self::Output {
        Self::new(-self.x, -self.y)
    }
}

macro_rules! ForNumbers {
    ($ty:ty) => {
        impl Mul<Vector2<$ty>> for $ty {
            type Output = Vector2<$ty>;

            #[inline]
            fn mul(self, rhs: Vector2<$ty>) -> Self::Output {
                rhs * self
            }
        }
    };
    ($type:ty, $($other:ty),+ $(,)?) => {
        ForNumbers!($type); $(ForNumbers!($other);)+
    }
}

ForNumbers!(f32, f64, i32);
