use std::{
    cmp::Ordering,
    fmt::Debug,
    ops::{Add, Div, Mul, Neg, Sub},
};

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
}

macro_rules! NumberFloat {
    (($type:ty, $name:ident)) => {
        impl Number for $type {
            const ZERO: Self = 0.0;
            const ONE: Self = 1.0;

            #[inline]
            fn is_nan(&self) -> bool {
                <$type>::is_nan(*self)
            }

            #[inline]
            fn abs(&self) -> Self {
                <$type>::abs(*self)
            }

            #[inline]
            fn sqrt(&self) -> Self {
                <$type>::sqrt(*self)
            }

            #[inline]
            fn floor(&self) -> Self {
                <$type>::floor(*self)
            }

            #[inline]
            fn ceil(&self) -> Self {
                <$type>::ceil(*self)
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
        }
    };
    (($type:ty, $name:ident), $(($other_type:ty, $other_name:ident)),+ $(,)?) => {
        NumberFloat!(($type, $name)); $(NumberFloat!(($other_type, $other_name));)+
    }
}

NumberFloat!((f32, f32), (f64, f64));

macro_rules! NumberInteger {
    (($type:ty, $name:ident)) => {
        impl Number for $type {
            const ZERO: Self = 0;
            const ONE: Self = 1;

            #[inline]
            fn is_nan(&self) -> bool {
                false
            }

            #[inline]
            fn abs(&self) -> Self {
                <$type>::abs(*self)
            }

            #[inline]
            fn sqrt(&self) -> Self {
                debug_assert!(self.abs() < (Self::ONE + Self::ONE).pow(f64::MANTISSA_DIGITS as _));

                (*self as f64).sqrt() as Self
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
        }
    };
    (($type:ty, $name:ident), $(($other_type:ty, $other_name:ident)),+ $(,)?) => {
        NumberFloat!(($type, $name)); $(NumberFloat!(($other_type, $other_name));)+
    }
}

NumberInteger!((i32, i32));
