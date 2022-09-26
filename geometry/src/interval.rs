use std::{
    marker::PhantomData,
    ops::{Add, Mul, Sub},
};

use crate::number::Number;

/// class to simplify calculations with intervals of real numbers
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Interval<T: Number> {
    low: T,
    high: T,
    _remove_constructors: PhantomData<()>,
}

impl<T: Number> Default for Interval<T> {
    fn default() -> Self {
        Self {
            low: T::ZERO,
            high: T::ZERO,
            _remove_constructors: PhantomData,
        }
    }
}

impl<T: Number> Interval<T> {
    /// Create a new interval between two points
    pub fn new(v0: T, v1: T) -> Self {
        debug_assert!(!v0.is_nan());
        debug_assert!(!v1.is_nan());

        Self {
            low: v0.min(v1),
            high: v0.max(v1),
            _remove_constructors: PhantomData,
        }
    }

    /// Create a new interval where the bounds are equal
    pub fn new_eq(val: T) -> Self {
        Self {
            low: val,
            high: val,
            _remove_constructors: PhantomData,
        }
    }

    /// Calculate the sin of the interval, output is undefined when the interval is
    /// outside of the bounds [0, 2pi]
    pub fn sin(&self) -> Self {
        debug_assert!(self.low >= T::ZERO);
        debug_assert!(self.high <= T::cast(2.0001) * T::PI);
        let mut sin_low = self.low.sin();
        let mut sin_high = self.high.sin();

        if sin_low > sin_high {
            std::mem::swap(&mut sin_low, &mut sin_high);
        }
        if self.low < T::PI / T::TWO && self.high > T::PI / T::TWO {
            sin_high = T::ONE;
        }
        if self.low < T::cast(3.0 / 2.0) * T::PI && self.high > T::cast(3.0 / 2.0) * T::PI {
            sin_low = -T::ONE;
        }

        Interval::new(sin_low, sin_high)
    }

    /// Calculate the cos of the interval, output is undefined when the interval is
    /// outside of the bounds [0, 2pi]
    pub fn cos(&self) -> Self {
        debug_assert!(self.low >= T::ZERO);
        debug_assert!(self.high <= T::cast(2.0001) * T::PI);

        let mut cos_low = self.low.cos();
        let mut cos_high = self.high.cos();

        if cos_low > cos_high {
            std::mem::swap(&mut cos_low, &mut cos_high);
        }
        if self.low < T::PI && self.high > T::PI {
            cos_low = -T::ONE;
        }

        Interval::new(cos_low, cos_high)
    }

    /// Find up to N motion derivative zeros within the this interval
    pub fn find_zeros(
        &self,
        c: [T; 5],
        theta: T,
        zeros: &mut [T],
        zero_count: &mut usize,
        depth: usize,
    ) {
        let [c1, c2, c3, c4, c5] = c;
        let range = Interval::new_eq(c1)
            + (Interval::new_eq(c2) + Interval::new_eq(c3) * *self)
                * (Interval::new_eq(T::TWO * theta) * *self).cos()
            + (Interval::new_eq(c4) + Interval::new_eq(c5) * *self)
                * (Interval::new_eq(T::TWO * theta) * *self).sin();
        if range.low > T::ZERO || range.high < T::ZERO || range.low == range.high {
            return;
        }

        if depth > 0 {
            let mid = (self.low + self.high) * T::HALF;
            Interval::new(self.low, mid).find_zeros(c, theta, zeros, zero_count, depth - 1);
            Interval::new(mid, self.high).find_zeros(c, theta, zeros, zero_count, depth - 1);
        } else {
            // use Newton's method to refine zero
            let mut t_newton = (self.low + self.high) * T::HALF;
            for _ in 0..4 {
                let f_newton = c1
                    + (c2 + c3 * t_newton) * (T::TWO * theta * t_newton).cos()
                    + (c4 + c5 * t_newton) * (T::TWO * theta * t_newton).sin();
                let f_prime_newton = (c3 + T::TWO * (c4 + c5 * t_newton) * theta)
                    * (T::TWO * t_newton * theta).cos()
                    + (c5 - T::TWO * (c2 + c3 * t_newton) * theta)
                        * (T::TWO * t_newton * theta).sin();
                if f_newton == T::ZERO || f_prime_newton == T::ZERO {
                    break;
                }
                t_newton = t_newton - f_newton / f_prime_newton;
            }
            if t_newton >= self.low - T::cast(1e-3) && t_newton < self.high + T::cast(1e-3) {
                zeros[*zero_count] = t_newton;
                *zero_count += 1;
            }
        }
    }
}

// NOTE: CPU rounding method should be round down for low bound, round up for high bound
// this is not done currently/at all
impl<T: Number> Add<Interval<T>> for Interval<T> {
    type Output = Self;

    fn add(self, rhs: Interval<T>) -> Self {
        Self::new(self.low + rhs.low, self.high + rhs.high)
    }
}

impl<T: Number> Sub<Interval<T>> for Interval<T> {
    type Output = Self;

    fn sub(self, rhs: Interval<T>) -> Self {
        Self::new(self.low - rhs.high, self.high - rhs.low)
    }
}

impl<T: Number> Mul<Interval<T>> for Interval<T> {
    type Output = Self;

    fn mul(self, rhs: Interval<T>) -> Self {
        let min = (self.low * rhs.low)
            .min(self.high * rhs.low)
            .min(self.low * rhs.high)
            .min(self.high * rhs.high);
        let max = (self.low * rhs.low)
            .max(self.high * rhs.low)
            .max(self.low * rhs.high)
            .max(self.high * rhs.high);
        Interval::new(min, max)
    }
}
