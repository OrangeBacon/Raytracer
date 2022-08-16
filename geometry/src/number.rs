use std::{
    cmp::Ordering,
    fmt::Debug,
    ops::{Add, Div, Mul, Neg, Sub},
};

/// any type that has a const default value
pub trait ConstZero {
    const ZERO: Self;
}

/// Any type that can be used as a number
pub trait Number:
    Debug
    + Copy
    + PartialEq
    + PartialOrd
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
{
    /// The Zero value for the type
    const ZERO: Self;

    /// The One value for the type
    const ONE: Self;

    /// The two value for the type
    const TWO: Self;

    /// The smallest value for this type
    const MIN: Self;

    /// The largest value for this type
    const MAX: Self;

    /// The closest value to pi for this type
    const PI: Self;

    /// Infinity or Self::MAX if there is no representable infinity
    const INFINITY: Self;

    /// 0.0005
    const LARGE_EPSILON: Self;

    /// Is this type a NaN value, always false if type is not floating
    fn is_nan(&self) -> bool;

    /// Absolute value of the input
    fn abs(&self) -> Self;

    /// Square root of the number
    fn sqrt(&self) -> Self;

    /// floor of the number
    fn floor(&self) -> Self;

    /// ceiling of the number
    fn ceil(&self) -> Self;

    /// Ordering of two values
    fn order(&self, rhs: &Self) -> Ordering;

    /// Cast to self
    fn cast<T: Number>(num: T) -> Self;

    /// Cast to double
    fn f64(&self) -> f64;

    /// Cast to float
    fn f32(&self) -> f32;

    /// Cast to integer
    fn i32(&self) -> i32;

    /// Clamp the value between two other values
    fn clamp(&self, min: Self, max: Self) -> Self;

    /// Sine of the number (radians)
    fn sin(&self) -> Self;

    /// Cos of the number (radians)
    fn cos(&self) -> Self;

    /// Inverse cos of the number (radians)
    fn acos(&self) -> Self;
}

/// Marker trait for integers
pub trait Integer: Number {}

impl<T: Number> ConstZero for T {
    const ZERO: Self = <T as Number>::ZERO;
}

macro_rules! NumberFloat {
    (($type:ty, $name:ident)) => {
        impl Number for $type {
            const ZERO: Self = 0.0;
            const ONE: Self = 1.0;
            const TWO: Self = 2.0;
            const MIN: Self = <$type>::MIN;
            const MAX: Self = <$type>::MAX;
            const PI: Self = std::$name::consts::PI;
            const INFINITY: Self = <$type>::INFINITY;
            const LARGE_EPSILON: Self = 0.0005;

            #[inline]
            fn is_nan(&self) -> bool {
                <$type>::is_nan(*self)
            }

            #[inline]
            fn order(&self, rhs: &Self) -> Ordering {
                self.total_cmp(rhs)
            }

            #[inline]
            fn cast<T: Number>(num: T) -> Self {
                num.$name()
            }

            #[inline]
            fn f64(&self) -> f64 {
                *self as _
            }

            #[inline]
            fn f32(&self) -> f32 {
                *self as _
            }

            #[inline]
            fn i32(&self) -> i32 {
                *self as _
            }

            #[inline]
            fn clamp(&self, min: Self, max: Self) -> Self {
                <$type>::clamp(*self, min, max)
            }

            NumberFloat! { @fns($type) sin, cos, acos, ceil, floor, sqrt, abs }
        }
    };
    (($type:ty, $name:ident), $(($other_type:ty, $other_name:ident)),+ $(,)?) => {
        NumberFloat!(($type, $name)); $(NumberFloat!(($other_type, $other_name));)+
    };
    (@fns($type:ty) $($name:ident),+ $(,)?) => {
        $(
            #[inline]
            fn $name(&self) -> Self {
                <$type>::$name(*self)
            }
        )+
    }
}

NumberFloat!((f32, f32), (f64, f64));

macro_rules! NumberInteger {
    (($type:ty, $name:ident)) => {
        impl Integer for $type {}
        impl Number for $type {
            const ZERO: Self = 0;
            const ONE: Self = 1;
            const TWO: Self = 2;
            const MIN: Self = <$type>::MIN;
            const MAX: Self = <$type>::MAX;
            const PI: Self = std::f64::consts::PI as _;
            const INFINITY: Self = Self::MAX;
            const LARGE_EPSILON: Self = 0;

            #[inline]
            fn is_nan(&self) -> bool {
                false
            }

            #[inline]
            fn abs(&self) -> Self {
                <$type>::abs(*self)
            }

            #[inline]
            fn floor(&self) -> Self {
                *self
            }

            #[inline]
            fn ceil(&self) -> Self {
                *self
            }

            #[inline]
            fn order(&self, rhs: &Self) -> Ordering {
                self.cmp(rhs)
            }

            #[inline]
            fn cast<T: Number>(num: T) -> Self {
                num.$name()
            }

            #[inline]
            fn f64(&self) -> f64 {
                *self as _
            }

            #[inline]
            fn f32(&self) -> f32 {
                *self as _
            }

            #[inline]
            fn i32(&self) -> i32 {
                *self as _
            }

            #[inline]
            fn clamp(&self, min: Self, max: Self) -> Self {
                <$type as Ord>::clamp(*self, min, max)
            }

            NumberInteger! { @fns($type) sqrt, sin, cos, acos }
        }
    };
    (($type:ty, $name:ident), $(($other_type:ty, $other_name:ident)),+ $(,)?) => {
        NumberFloat!(($type, $name)); $(NumberFloat!(($other_type, $other_name));)+
    };
    (@fns($type:ty) $($name:ident),+ $(,)?) => {
        $(
            #[inline]
            fn $name(&self) -> Self {
                debug_assert!(self.abs() < (Self::ONE + Self::ONE).pow(f64::MANTISSA_DIGITS as _));

                (*self as f64).$name() as Self
            }
        )+
    }
}

NumberInteger!((i32, i32));
