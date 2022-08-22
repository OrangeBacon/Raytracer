use std::{
    marker::PhantomData,
    ops::{Add, Index, Mul, Sub},
};

use crate::number::Number;

/// class to simplify calculations with intervals of real numbers
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
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
    pub fn find_zeros<const N: usize>(
        &self,
        c1: T,
        c2: T,
        c3: T,
        c4: T,
        c5: T,
        theta: T,
    ) -> (usize, [T; N]) {
        #[derive(Clone, Copy)]
        struct Params<T: Number>([T; 5]);
        impl<T: Number> Index<usize> for Params<T> {
            type Output = T;

            fn index(&self, index: usize) -> &Self::Output {
                &self.0[index - 1]
            }
        }

        fn inner<const N: usize, T: Number>(
            interval: Interval<T>,
            c: Params<T>,
            theta: T,
            depth: i32,
            res: &mut [T; N],
            result_count: &mut usize,
        ) {
            if depth < 0 || *result_count == N {
                return;
            }

            let n = Interval::new_eq;
            let range = n(c[1])
                + (n(c[2]) + n(c[3]) * interval) * (n(T::TWO * theta) * interval).cos()
                + (n(c[4]) + n(c[5]) * interval) * (n(T::TWO * theta) * interval).sin();

            if range.low > T::ZERO || range.high < T::ZERO || range.low == range.high {
                return;
            }

            if depth > 0 {
                // split interval and check both parts
                let mid = (interval.low + interval.high) * (T::ONE / T::TWO);
                inner(
                    Interval::new(interval.low, mid),
                    c,
                    theta,
                    depth - 1,
                    res,
                    result_count,
                );
                inner(
                    Interval::new(mid, interval.high),
                    c,
                    theta,
                    depth - 1,
                    res,
                    result_count,
                );
            } else {
                // refine zero using newton's method
                let mut t_newton = (interval.low + interval.high) * (T::ONE / T::TWO);
                for _ in 0..4 {
                    let f_newton = c[1]
                        + (c[2] + c[3] * t_newton) * (T::TWO * theta * t_newton).cos()
                        + (c[4] + c[5] * t_newton) * (T::TWO * theta * t_newton).sin();
                    let f_prime_newton = (c[3] + T::TWO * (c[4] + c[5] * t_newton) * theta)
                        * (T::TWO * t_newton * theta).cos()
                        + (c[5] - T::TWO * (c[2] + c[3] * t_newton) * theta)
                            * (T::TWO * t_newton * theta).sin();
                    if f_newton == T::ZERO || f_prime_newton == T::ZERO {
                        return;
                    }
                    t_newton -= f_newton / f_prime_newton;
                }
                if t_newton >= interval.low - T::cast(1e-3)
                    && t_newton < interval.high + T::cast(1e-3)
                {
                    res[*result_count] = t_newton;
                    *result_count += 1;
                }
            }
        }

        let mut result = [T::ZERO; N];
        let mut result_count = 0;
        inner(
            *self,
            Params([c1, c2, c3, c4, c5]),
            theta,
            8,
            &mut result,
            &mut result_count,
        );
        (result_count, result)
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
