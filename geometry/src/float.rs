//! Functions for dealing with floating point error

use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::Number;

/// Calculate the next float after the input
pub fn next_float_up<T: Number>(v: T) -> T {
    if !v.is_finite() && v > T::ZERO {
        return v;
    }
    if v == -T::ZERO {
        return T::ZERO;
    }

    let mut bits = v.to_bits();
    if v >= T::ZERO {
        bits += <T as Number>::Bits::ONE;
    } else {
        bits -= <T as Number>::Bits::ONE;
    }
    T::from_bits(bits)
}

/// Calculate the next float before the input
pub fn next_float_down<T: Number>(v: T) -> T {
    if !v.is_finite() && v < T::ZERO {
        return v;
    }
    if v == T::ZERO {
        return -T::ZERO;
    }

    let mut bits = v.to_bits();
    if v > T::ZERO {
        bits += <T as Number>::Bits::ONE;
    } else {
        bits -= <T as Number>::Bits::ONE;
    }
    T::from_bits(bits)
}

/// Gamma floating point error bound
pub fn gamma<T: Number>(n: i32) -> T {
    (T::cast(n) * T::EPSILON) / (T::ONE - T::cast(n) * T::EPSILON)
}

/// Float that tracks its error due to floating point inaccuracies
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct EFloat<T: Number> {
    value: T,
    low: T,
    high: T,

    #[cfg(debug_assertions)]
    ld: f64,
}

impl<T: Number> EFloat<T> {
    /// Create a new Error float with 0 error
    pub fn new(v: T) -> Self {
        Self {
            value: v,
            low: v,
            high: v,
            #[cfg(debug_assertions)]
            ld: v.f64(),
        }
        .assert_correct()
    }

    /// Create a new Error float with a given error bound
    pub fn new_with_err(v: T, err: T) -> Self {
        if err == T::ZERO {
            return Self::new(v);
        }

        Self {
            value: v,
            low: next_float_down(v - err),
            high: next_float_up(v + err),

            #[cfg(debug_assertions)]
            ld: v.f64(),
        }
        .assert_correct()
    }

    /// The lower bound of the error stored in the float
    pub fn lower_bound(&self) -> T {
        self.low
    }

    /// The upper bound of the error stored in the float
    pub fn upper_bound(&self) -> T {
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
            debug_assert!(self.low.f64() <= self.ld);
            debug_assert!(self.ld <= self.high.f64());
        }

        *self
    }

    /// The absolute value of the error
    pub fn absolute_error(&self) -> T {
        next_float_up(
            (self.high - self.value)
                .abs()
                .max((self.value - self.low).abs()),
        )
    }

    #[cfg(debug_assertions)]
    pub fn relative_error(&self) -> f64 {
        ((self.ld - self.value.f64()) / self.ld).abs()
    }

    /// Get the contained floating point value, discarding the error bounds
    pub fn value(&self) -> T {
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
        if self.low >= T::ZERO {
            *self // All above zero so nothing changes
        } else if self.high <= T::ZERO {
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
                low: T::ZERO,
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
    pub fn quadratic(a: Self, b: Self, c: Self) -> Option<(Self, Self)> {
        let discriminant = b.value * b.value - T::cast(4) * a.value * c.value;
        if discriminant < T::ZERO {
            return None;
        }
        let root = discriminant.sqrt();
        let root = EFloat::new_with_err(root, T::EPSILON * (T::ONE / T::ZERO) * root);

        let q = if b.value < T::ZERO {
            (b - root) * -(T::HALF)
        } else {
            (b + root) * -(T::HALF)
        };

        let t0 = q / a;
        let t1 = c / q;

        if t0.value() > t1.value() {
            Some((t1, t0))
        } else {
            Some((t0, t1))
        }
    }
}

impl<T: Number> Add for EFloat<T> {
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

impl<T: Number> AddAssign for EFloat<T> {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs
    }
}

impl<T: Number> Sub for EFloat<T> {
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

impl<T: Number> SubAssign for EFloat<T> {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs
    }
}

impl<T: Number> Mul for EFloat<T> {
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

impl<T: Number> MulAssign for EFloat<T> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs
    }
}

impl<T: Number> Div for EFloat<T> {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        let (low, high) = if rhs.low < T::ZERO && rhs.high > T::ZERO {
            // denominator is approx 0 so say infinite error
            (-T::INFINITY, T::INFINITY)
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

impl<T: Number> DivAssign for EFloat<T> {
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs
    }
}

impl<T: Number> Neg for EFloat<T> {
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

impl<F: Number> Add<F> for EFloat<F> {
    type Output = Self;

    fn add(self, rhs: F) -> Self::Output {
        self + EFloat::new(rhs)
    }
}

impl<F: Number> AddAssign<F> for EFloat<F> {
    fn add_assign(&mut self, rhs: F) {
        *self = *self + rhs
    }
}

impl<F: Number> Sub<F> for EFloat<F> {
    type Output = Self;

    fn sub(self, rhs: F) -> Self::Output {
        self - EFloat::new(rhs)
    }
}

impl<F: Number> SubAssign<F> for EFloat<F> {
    fn sub_assign(&mut self, rhs: F) {
        *self = *self - rhs
    }
}

impl<F: Number> Mul<F> for EFloat<F> {
    type Output = Self;

    fn mul(self, rhs: F) -> Self::Output {
        self * EFloat::new(rhs)
    }
}

impl<F: Number> MulAssign<F> for EFloat<F> {
    fn mul_assign(&mut self, rhs: F) {
        *self = *self * rhs
    }
}

impl<F: Number> Div<F> for EFloat<F> {
    type Output = Self;

    fn div(self, rhs: F) -> Self::Output {
        self / EFloat::new(rhs)
    }
}

impl<F: Number> DivAssign<F> for EFloat<F> {
    fn div_assign(&mut self, rhs: F) {
        *self = *self / rhs
    }
}

macro_rules! number_impls {
    ($lhs:ty) => {
        impl Add<EFloat<$lhs>> for $lhs {
            type Output = EFloat<$lhs>;

            fn add(self, rhs: EFloat<$lhs>) -> Self::Output {
                EFloat::new(self) + rhs
            }
        }

        impl Sub<EFloat<$lhs>> for $lhs {
            type Output = EFloat<$lhs>;

            fn sub(self, rhs: EFloat<$lhs>) -> Self::Output {
                EFloat::new(self) - rhs
            }
        }

        impl Mul<EFloat<$lhs>> for $lhs {
            type Output = EFloat<$lhs>;

            fn mul(self, rhs: EFloat<$lhs>) -> Self::Output {
                EFloat::new(self) * rhs
            }
        }

        impl Div<EFloat<$lhs>> for $lhs {
            type Output = EFloat<$lhs>;

            fn div(self, rhs: EFloat<$lhs>) -> Self::Output {
                EFloat::new(self) / rhs
            }
        }
    };

    ($($lhs:ty),+ $(,)?) => {
        $(number_impls!($lhs);)+
    }
}

number_impls!(f32, f64, i32);
