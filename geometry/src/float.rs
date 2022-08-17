//! Functions for dealing with floating point error

use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::{Float, Number};

const EPSILON: Float = Float::EPSILON * 0.5;

/// Convert a floating point value into its bit representation.
/// Safety: Assumes that the bit representation of the input is a valid value
/// of the output type.
#[inline]
pub unsafe fn float_to_bits(f: Float) -> <Float as Number>::Bits {
    std::mem::transmute_copy(&f)
}

/// Convert the bits of a float into a float.
/// Safety: Assumes that the bit representation of the input is a valid value
/// of the output type.
#[inline]
pub unsafe fn bits_to_float(f: <Float as Number>::Bits) -> Float {
    std::mem::transmute_copy(&f)
}

/// Calculate the next float after the input
pub fn next_float_up(v: Float) -> Float {
    if v.is_infinite() && v > 0.0 {
        return v;
    }
    if v == -0.0 {
        return 0.0;
    }

    let mut bits = unsafe { float_to_bits(v) };
    if v >= 0.0 {
        bits += 1;
    } else {
        bits -= 1;
    }
    unsafe { bits_to_float(bits) }
}

/// Calculate the next float before the input
pub fn next_float_down(v: Float) -> Float {
    if v.is_infinite() && v < 0.0 {
        return v;
    }
    if v == 0.0 {
        return -0.0;
    }

    let mut bits = unsafe { float_to_bits(v) };
    if v > 0.0 {
        bits -= 1;
    } else {
        bits += 1;
    }
    unsafe { bits_to_float(bits) }
}

/// Gamma floating point error bound
pub fn gamma(n: i32) -> Float {
    (n as Float * EPSILON) / (1.0 - n as Float * EPSILON)
}

/// Float that tracks its error due to floating point inaccuracies
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct EFloat {
    value: Float,
    low: Float,
    high: Float,

    #[cfg(debug_assertions)]
    ld: f64,
}

impl EFloat {
    /// Create a new Error float with 0 error
    pub fn new(v: Float) -> Self {
        Self {
            value: v,
            low: v,
            high: v,
            #[cfg(debug_assertions)]
            ld: v as _,
        }
        .assert_correct()
    }

    /// Create a new Error float with a given error bound
    pub fn new_with_err(v: Float, err: Float) -> Self {
        if err == 0.0 {
            return Self::new(v);
        }

        Self {
            value: v,
            low: next_float_down(v - err),
            high: next_float_up(v + err),

            #[cfg(debug_assertions)]
            ld: v as _,
        }
        .assert_correct()
    }

    /// The lower bound of the error stored in the float
    pub fn lower_bound(&self) -> Float {
        self.low
    }

    /// The upper bound of the error stored in the float
    pub fn upper_bound(&self) -> Float {
        self.high
    }

    /// Assert that this error float is valid
    fn assert_correct(&self) -> Self {
        if self.low.is_finite()
            && !self.low.is_nan()
            && self.high.is_finite()
            && !self.high.is_nan()
        {
            debug_assert!(self.low <= self.high);
        }

        #[cfg(debug_assertions)]
        if self.value.is_finite() && !self.value.is_nan() {
            debug_assert!(self.low as f64 <= self.ld);
            debug_assert!(self.ld <= self.high as f64);
        }

        *self
    }

    /// The absolute value of the error
    pub fn absolute_error(&self) -> Float {
        next_float_up(
            (self.high - self.value)
                .abs()
                .max((self.value - self.low).abs()),
        )
    }

    #[cfg(debug_assertions)]
    pub fn relative_error(&self) -> Float {
        ((self.ld - self.value as f64) / self.ld).abs() as _
    }

    /// Get the contained floating point value, discarding the error bounds
    pub fn value(&self) -> Float {
        self.value
    }

    /// Square root of the float and calculate error bounds
    pub fn sqrt(&self) -> Self {
        Self {
            value: self.value.sqrt(),
            low: next_float_down(self.low.sqrt()),
            high: next_float_up(self.high.sqrt()),

            #[cfg(debug_assertions)]
            ld: self.ld.sqrt(),
        }
        .assert_correct()
    }

    /// Absolute value of the float, including dealing with the bounds
    pub fn abs(&self) -> Self {
        if self.low >= 0.0 {
            *self // All above zero so nothing changes
        } else if self.high <= 0.0 {
            Self {
                // All bellow zero so negate
                value: -self.value,
                low: -self.high,
                high: -self.low,

                #[cfg(debug_assertions)]
                ld: -self.ld,
            }
            .assert_correct()
        } else {
            // Interval approx equal zero
            Self {
                value: self.value.abs(),
                low: 0.0,
                high: self.high.max(-self.low),

                #[cfg(debug_assertions)]
                ld: self.ld.abs(),
            }
            .assert_correct()
        }
    }

    /// Solve a quadratic equation ax^2 + bx + c = 0
    /// If no real solutions are found, returns None.  If both solutions are the
    /// same, that solution will be both return values.
    pub fn quadratic(a: EFloat, b: EFloat, c: EFloat) -> Option<(EFloat, EFloat)> {
        let discriminant = b.value * b.value - 4.0 * a.value * c.value;
        if discriminant < 0.0 {
            return None;
        }
        let root = discriminant.sqrt();
        let root = EFloat::new_with_err(root, EPSILON * root);

        let q = if b.value < 0.0 {
            -0.5 * (b - root)
        } else {
            -0.5 * (b + root)
        };

        let t0: EFloat = q / a;
        let t1: EFloat = c / q;

        if t0.value() > t1.value() {
            Some((t1, t0))
        } else {
            Some((t0, t1))
        }
    }
}

impl Add for EFloat {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            value: self.value + rhs.value,
            low: next_float_down(self.low + rhs.low),
            high: next_float_up(self.high + rhs.high),

            #[cfg(debug_assertions)]
            ld: self.ld + rhs.ld,
        }
        .assert_correct()
    }
}

impl AddAssign for EFloat {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs
    }
}

impl Sub for EFloat {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            value: self.value - rhs.value,
            low: next_float_down(self.low - rhs.low),
            high: next_float_up(self.high - rhs.high),

            #[cfg(debug_assertions)]
            ld: self.ld - rhs.ld,
        }
        .assert_correct()
    }
}

impl SubAssign for EFloat {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs
    }
}

impl Mul for EFloat {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let prod = [
            self.low * rhs.low,
            self.high * rhs.low,
            self.low * rhs.high,
            self.high * rhs.high,
        ];
        let min = prod[0].min(prod[1]).min(prod[2]).min(prod[3]);
        let max = prod[0].max(prod[1]).max(prod[2]).max(prod[3]);

        Self {
            value: self.value * rhs.value,
            low: next_float_down(min),
            high: next_float_up(max),

            #[cfg(debug_assertions)]
            ld: self.ld * rhs.ld,
        }
        .assert_correct()
    }
}

impl MulAssign for EFloat {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs
    }
}

impl Div for EFloat {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        let (low, high) = if rhs.low < 0.0 && rhs.high > 0.0 {
            // denominator is approx 0 so say infinite error
            (-Float::INFINITY, Float::INFINITY)
        } else {
            let div = [
                self.low / rhs.low,
                self.high / rhs.low,
                self.low / rhs.high,
                self.high / rhs.high,
            ];
            let min = div[0].min(div[1]).min(div[2]).min(div[3]);
            let max = div[0].max(div[1]).max(div[2]).max(div[3]);

            (next_float_down(min), next_float_up(max))
        };

        Self {
            value: self.value / rhs.value,
            low,
            high,

            #[cfg(debug_assertions)]
            ld: self.ld / rhs.ld,
        }
        .assert_correct()
    }
}

impl DivAssign for EFloat {
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs
    }
}

impl Neg for EFloat {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            value: -self.value,
            low: -self.high,
            high: -self.low,

            #[cfg(debug_assertions)]
            ld: -self.ld,
        }
        .assert_correct()
    }
}

impl<F: Number> Add<F> for EFloat {
    type Output = EFloat;

    fn add(self, rhs: F) -> Self::Output {
        self + EFloat::new(rhs.f32())
    }
}

impl<F: Number> AddAssign<F> for EFloat {
    fn add_assign(&mut self, rhs: F) {
        *self = *self + rhs
    }
}

impl<F: Number> Sub<F> for EFloat {
    type Output = EFloat;

    fn sub(self, rhs: F) -> Self::Output {
        self - EFloat::new(rhs.f32())
    }
}

impl<F: Number> SubAssign<F> for EFloat {
    fn sub_assign(&mut self, rhs: F) {
        *self = *self - rhs
    }
}

impl<F: Number> Mul<F> for EFloat {
    type Output = EFloat;

    fn mul(self, rhs: F) -> Self::Output {
        self * EFloat::new(rhs.f32())
    }
}

impl<F: Number> MulAssign<F> for EFloat {
    fn mul_assign(&mut self, rhs: F) {
        *self = *self * rhs
    }
}

impl<F: Number> Div<F> for EFloat {
    type Output = EFloat;

    fn div(self, rhs: F) -> Self::Output {
        self / EFloat::new(rhs.f32())
    }
}

impl<F: Number> DivAssign<F> for EFloat {
    fn div_assign(&mut self, rhs: F) {
        *self = *self / rhs
    }
}

macro_rules! number_impls {
    ($lhs:ty) => {
        impl Add<EFloat> for $lhs {
            type Output = EFloat;

            fn add(self, rhs: EFloat) -> Self::Output {
                EFloat::new(self as Float) + rhs
            }
        }

        impl Sub<EFloat> for $lhs {
            type Output = EFloat;

            fn sub(self, rhs: EFloat) -> Self::Output {
                EFloat::new(self as Float) - rhs
            }
        }

        impl Mul<EFloat> for $lhs {
            type Output = EFloat;

            fn mul(self, rhs: EFloat) -> Self::Output {
                EFloat::new(self as Float) * rhs
            }
        }

        impl Div<EFloat> for $lhs {
            type Output = EFloat;

            fn div(self, rhs: EFloat) -> Self::Output {
                EFloat::new(self as Float) / rhs
            }
        }
    };

    ($($lhs:ty),+ $(,)?) => {
        $(number_impls!($lhs);)+
    }
}

number_impls!(f32, f64, i32);
