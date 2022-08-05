use std::{
    f32::consts::PI,
    ops::{Add, Index, Mul, Sub},
};

use crate::Float;

/// class to simplify calculations with intervals of real numbers
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Interval {
    low: Float,
    high: Float,
}

impl Default for Interval {
    fn default() -> Self {
        Self {
            low: 0.0,
            high: 0.0,
        }
    }
}

impl Interval {
    /// Create a new interval between two points
    pub fn new(v0: Float, v1: Float) -> Self {
        Self {
            low: v0.min(v1),
            high: v0.max(v1),
        }
    }

    /// Create a new interval where the bounds are equal
    pub fn new_eq(val: Float) -> Self {
        Self {
            low: val,
            high: val,
        }
    }

    /// Calculate the sin of the interval, output is undefined when the interval is
    /// outside of the bounds [0, 2pi]
    pub fn sin(&self) -> Self {
        let mut sin_low = self.low.sin();
        let mut sin_high = self.high.sin();

        if sin_low > sin_high {
            std::mem::swap(&mut sin_low, &mut sin_high);
        }
        if self.low < PI / 2.0 && self.high > PI / 2.0 {
            sin_high = 1.0;
        }
        if self.low < (3.0 / 2.0) * PI && self.high > (3.0 / 2.0) * PI {
            sin_low = -1.0;
        }

        Interval::new(sin_low, sin_high)
    }

    /// Calculate the cos of the interval, output is undefined when the interval is
    /// outside of the bounds [0, 2pi]
    pub fn cos(&self) -> Self {
        let mut cos_low = self.low.cos();
        let mut cos_high = self.high.cos();

        if cos_low > cos_high {
            std::mem::swap(&mut cos_low, &mut cos_high);
        }
        if self.low < PI && self.high > PI {
            cos_low = -1.0;
        }

        Interval::new(cos_low, cos_high)
    }

    /// Find up to N motion derivative zeros within the this interval
    pub fn find_zeros<const N: usize>(
        &self,
        c1: Float,
        c2: Float,
        c3: Float,
        c4: Float,
        c5: Float,
        theta: Float,
    ) -> (usize, [Float; N]) {
        #[derive(Clone, Copy)]
        struct Params([Float; 5]);
        impl Index<usize> for Params {
            type Output = Float;

            fn index(&self, index: usize) -> &Self::Output {
                &self.0[index - 1]
            }
        }

        fn inner<const N: usize>(
            interval: Interval,
            c: Params,
            theta: Float,
            depth: i32,
            res: &mut [Float; N],
            result_count: &mut usize,
        ) {
            if depth < 0 || *result_count == N {
                return;
            }

            let n = Interval::new_eq;
            let range = n(c[1])
                + (n(c[2]) + n(c[3]) * interval) * (n(2.0 * theta) * interval).cos()
                + (n(c[4]) + n(c[5]) * interval) * (n(2.0 * theta) * interval).sin();

            if range.low > 0.0 || range.high < 0.0 || range.low == range.high {
                return;
            }

            if depth > 0 {
                // split interval and check both parts
                let mid = (interval.low + interval.high) / 2.0;
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
                let mut t_newton = (interval.low + interval.high) / 2.0;
                for _ in 0..4 {
                    let f_newton = c[1]
                        + (c[2] + c[3] * t_newton) * (2.0 * theta * t_newton).cos()
                        + (c[4] + c[5] * t_newton) * (2.0 * theta * t_newton).sin();
                    let f_prime_newton = (c[3] + 2.0 * (c[4] + c[5] * t_newton) * theta)
                        * (2.0 * t_newton * theta).cos()
                        + (c[5] - 2.0 * (c[2] + c[3] * t_newton) * theta)
                            * (2.0 * t_newton * theta).sin();
                    if f_newton == 0.0 || f_prime_newton == 0.0 {
                        return;
                    }
                    t_newton = t_newton - f_newton / f_prime_newton;
                }
                res[*result_count] = t_newton;
                *result_count += 1;
            }
        }

        let mut result = [0.0; N];
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
impl Add<Interval> for Interval {
    type Output = Self;

    fn add(self, rhs: Interval) -> Self {
        Self::new(self.low + rhs.low, self.high + rhs.high)
    }
}

impl Sub<Interval> for Interval {
    type Output = Self;

    fn sub(self, rhs: Interval) -> Self {
        Self::new(self.low - rhs.high, self.high - rhs.low)
    }
}

impl Mul<Interval> for Interval {
    type Output = Self;

    fn mul(self, rhs: Interval) -> Self {
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
