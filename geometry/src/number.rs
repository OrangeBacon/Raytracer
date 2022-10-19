use std::{
    cmp::Ordering,
    fmt::Debug,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
    sync::atomic::{AtomicI32, AtomicI64},
};

use crate::parallel::{AtomicF32, AtomicF64, AtomicNumber};

/// any type that has a const default value
pub trait ConstZero {
    const ZERO: Self;
}

/// Any type that can be used as a number
pub trait Number:
    Send
    + Sync
    + Debug
    + Copy
    + Default
    + PartialEq
    + PartialOrd
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
    + AddAssign<Self>
    + SubAssign<Self>
    + MulAssign<Self>
    + DivAssign<Self>
    + 'static
{
    /// The Bit representation of the type
    type Bits: Number;

    /// The atomic version of this type
    type Atomic: AtomicNumber<Self>;

    /// The Zero value for the type
    const ZERO: Self;

    /// The One value for the type
    const ONE: Self;

    /// The two value for the type
    const TWO: Self;

    /// Self::ONE / Self::TWO
    const HALF: Self;

    /// The smallest value for this type
    const MIN: Self;

    /// The largest value for this type
    const MAX: Self;

    /// The closest value to pi for this type
    const PI: Self;

    /// Infinity or Self::MAX if there is no representable infinity
    const INFINITY: Self;

    /// Machine epsilon value
    const EPSILON: Self;

    /// 0.0005
    const LARGE_EPSILON: Self;

    /// square root of 2, ~= 1.414
    const SQRT_2: Self;

    /// The largest floating point value that is less than one
    const ONE_MINUS_EPSILON: Self;

    /// Is this type a NaN value, always false if type is not floating
    fn is_nan(&self) -> bool;

    /// Is this number not infinity
    fn is_finite(&self) -> bool;

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
    fn cast<T: Number>(x: T) -> Self;

    /// Cast to double
    fn f64(&self) -> f64;

    /// Cast to float
    fn f32(&self) -> f32;

    /// Cast to integer
    fn i32(&self) -> i32;

    /// Cast to 64 bit integer
    fn i64(&self) -> i64;

    /// Clamp the value between two other values
    fn clamp(&self, min: Self, max: Self) -> Self;

    /// Largest of the input and self
    #[inline]
    fn max(&self, rhs: Self) -> Self {
        if *self > rhs {
            *self
        } else {
            rhs
        }
    }

    /// Smallest of the input and self
    #[inline]
    fn min(&self, rhs: Self) -> Self {
        if *self < rhs {
            *self
        } else {
            rhs
        }
    }

    /// Sine of the number (radians)
    fn sin(&self) -> Self;

    /// Cos of the number (radians)
    fn cos(&self) -> Self;

    /// Inverse cos of the number (radians)
    fn acos(&self) -> Self;

    /// Base 2 logarithm of the number
    fn log2(&self) -> Self;

    /// Base 10 logarithm of the number
    fn log10(&self) -> Self;

    /// Natural logarithm of the number
    fn ln(&self) -> Self;

    /// e ^ self
    fn exp(&self) -> Self;

    /// Round to the nearest integer
    fn round(&self) -> Self;

    /// Convert a value to its [`Self::Bits`] type
    fn to_bits(&self) -> Self::Bits;

    /// Convert a value from its bits to self
    fn from_bits(bits: Self::Bits) -> Self;

    /// Convert a number in degrees to radians
    fn to_radians(&self) -> Self;

    /// Four quadrant inverse tangent of rhs / self
    fn atan2(&self, rhs: Self) -> Self;

    /// This number to the power of the argument
    fn pow(&self, exp: Self) -> Self;

    /// Returns the logarithm of the number in a given base
    fn log(&self, base: Self) -> Self;
}

/// Marker trait for integers
pub trait Integer: Number {}

impl<T: Number> ConstZero for T {
    const ZERO: Self = <T as Number>::ZERO;
}

macro_rules! NumberFloat {
    (($type:ty, $name:ident, $bits:ty, $atomic:ty)) => {
        impl Number for $type {
            type Bits = $bits;
            type Atomic = $atomic;
            const ZERO: Self = 0.0;
            const ONE: Self = 1.0;
            const TWO: Self = 2.0;
            const MIN: Self = <$type>::MIN;
            const MAX: Self = <$type>::MAX;
            const PI: Self = std::$name::consts::PI;
            const INFINITY: Self = <$type>::INFINITY;
            const EPSILON: Self = <$type>::EPSILON;
            const LARGE_EPSILON: Self = 0.0005;
            const HALF: Self = 0.5;
            const SQRT_2: Self = std::$name::consts::SQRT_2;
            const ONE_MINUS_EPSILON: Self = 1.0 - <$type>::EPSILON;

            #[inline]
            fn is_nan(&self) -> bool {
                <$type>::is_nan(*self)
            }

            #[inline]
            fn is_finite(&self) -> bool {
                <$type>::is_finite(*self)
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
            fn i64(&self) -> i64 {
                *self as _
            }

            #[inline]
            fn clamp(&self, min: Self, max: Self) -> Self {
                <$type>::clamp(*self, min, max)
            }

            #[inline]
            fn to_bits(&self) -> Self::Bits {
                <$type>::to_bits(*self) as _
            }

            #[inline]
            fn from_bits(bits: Self::Bits) -> Self {
                <$type>::from_bits(bits as _)
            }

            #[inline]
            fn atan2(&self, rhs: Self) -> Self {
                <$type>::atan2(*self, rhs)
            }

            #[inline]
            fn pow(&self, exp: Self) -> Self {
                <$type>::powf(*self, exp)
            }

            #[inline]
            fn log(&self, base: Self) -> Self {
                <$type>::log(*self, base)
            }

            NumberFloat! { @fns($type) sin, cos, acos, ceil, floor, sqrt, abs, to_radians, log2, log10, ln, round, exp }
        }
    };
    (($type:ty, $name:ident, $bits:ty, $atomic:ty), $(($other_type:ty, $other_name:ident, $other_bits:ty, $other_atomic:ty)),+ $(,)?) => {
        NumberFloat!(($type, $name, $bits, $atomic)); $(NumberFloat!(($other_type, $other_name, $other_bits, $other_atomic));)+
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

NumberFloat!((f32, f32, i32, AtomicF32), (f64, f64, i64, AtomicF64));

macro_rules! NumberInteger {
    (($type:ty, $name:ident, $atomic:ty)) => {
        impl Integer for $type {}
        impl Number for $type {
            type Bits = Self;
            type Atomic = $atomic;
            const ZERO: Self = 0;
            const ONE: Self = 1;
            const TWO: Self = 2;
            const MIN: Self = <$type>::MIN;
            const MAX: Self = <$type>::MAX;
            const PI: Self = std::f64::consts::PI as _;
            const INFINITY: Self = Self::MAX;
            const EPSILON: Self = 0;
            const LARGE_EPSILON: Self = 0;
            const HALF: Self = 0.5 as _;
            const SQRT_2: Self = std::f64::consts::SQRT_2 as _;
            const ONE_MINUS_EPSILON: Self = (1.0 - f64::EPSILON) as _;

            #[inline]
            fn is_nan(&self) -> bool {
                false
            }

            #[inline]
            fn is_finite(&self) -> bool {
                true
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
            fn round(&self) -> Self {
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
            fn i64(&self) -> i64 {
                *self as _
            }

            #[inline]
            fn clamp(&self, min: Self, max: Self) -> Self {
                <$type as Ord>::clamp(*self, min, max)
            }

            #[inline]
            fn to_bits(&self) -> Self::Bits {
                // SAFETY: by definition, all numbers are valid bit patterns for
                // other numbers
                unsafe { std::mem::transmute(*self) }
            }

            #[inline]
            fn from_bits(bits: Self::Bits) -> Self {
                // SAFETY: by definition, all numbers are valid bit patterns for
                // other numbers
                unsafe { std::mem::transmute(bits) }
            }

            #[inline]
            fn atan2(&self, rhs: Self) -> Self {
                debug_assert!(self.abs() < Self::TWO.pow(f64::MANTISSA_DIGITS as _));
                debug_assert!(rhs.abs() < Self::TWO.pow(f64::MANTISSA_DIGITS as _));

                (*self as f64).atan2(rhs as f64) as Self
            }

            #[inline]
            fn pow(&self, exp: Self) -> Self {
                debug_assert!(self.abs() < Self::TWO.pow(f64::MANTISSA_DIGITS as _));
                debug_assert!(exp.abs() < Self::TWO.pow(f64::MANTISSA_DIGITS as _));

                (*self as f64).pow(exp as f64) as Self
            }

            #[inline]
            fn log(&self, base: Self) -> Self {
                debug_assert!(self.abs() < Self::TWO.pow(f64::MANTISSA_DIGITS as _));
                debug_assert!(base.abs() < Self::TWO.pow(f64::MANTISSA_DIGITS as _));

                (*self as f64).log(base as f64) as Self
            }

            NumberInteger! { @fns($type) sqrt, sin, cos, acos, to_radians, log2, log10, ln, exp }
        }
    };
    (($type:ty, $name:ident, $atomic:ty), $(($other_type:ty, $other_name:ident, $other_atomic:ty)),+ $(,)?) => {
        NumberInteger!(($type, $name, $atomic)); $(NumberInteger!(($other_type, $other_name, $other_atomic));)+
    };
    (@fns($type:ty) $($name:ident),+ $(,)?) => {
        $(
            #[inline]
            fn $name(&self) -> Self {
                debug_assert!(self.abs() < Self::TWO.pow(f64::MANTISSA_DIGITS as _));

                (*self as f64).$name() as Self
            }
        )+
    }
}

NumberInteger!((i32, i32, AtomicI32), (i64, i64, AtomicI64));
