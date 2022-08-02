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

    /// Ordering of two values
    fn order(&self, rhs: &Self) -> Ordering;
}

macro_rules! NumberFloat {
    ($type:ty) => {
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
            fn order(&self, rhs: &Self) -> Ordering {
                self.total_cmp(rhs)
            }
        }
    };
    ($type:ty, $($other:ty),+ $(,)?) => {
        NumberFloat!($type); $(NumberFloat!($other);)+
    }
}

NumberFloat!(f32, f64);

macro_rules! NumberInteger {
    ($type:ty) => {
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
            fn order(&self, rhs: &Self) -> Ordering {
                self.cmp(rhs)
            }
        }
    };
    ($type:ty, $($other:ty),+ $(,)?) => {
        NumberInteger!($type); $(NumberInteger!($other);)+
    }
}

NumberInteger!(i32);
