//! Functions for dealing with floating point error

use std::ops::{Add, Sub, Mul, Div};

use crate::{Float, FloatBits};

const EPSILON: Float = Float::EPSILON * 0.5;

/// Convert a floating point value into its bit representation
pub fn float_to_bits(f: Float) -> FloatBits {
    unsafe { std::mem::transmute_copy(&f) }
}

/// Convert the bits of a float into a float
pub fn bits_to_float(f: FloatBits) -> Float {
    unsafe { std::mem::transmute_copy(&f) }
}

/// Calculate the next float after the input
pub fn next_float_up(v: Float) -> Float {
    if v.is_infinite() && v > 0.0 {
        return v;
    }
    if v == -0.0 {
        return 0.0;
    }

    let mut bits = float_to_bits(v);
    if v >= 0.0 {
        bits += 1;
    } else {
        bits -= 1;
    }
    bits_to_float(bits)
}

/// Calculate the next float before the input
pub fn next_float_down(v: Float) -> Float {
    if v.is_infinite() && v < 0.0 {
        return v;
    }
    if v == 0.0 {
        return -0.0;
    }

    let mut bits = float_to_bits(v);
    if v > 0.0 {
        bits -= 1;
    } else {
        bits += 1;
    }
    bits_to_float(bits)
}

/// Gamma floating point error bound
pub fn gamma(n: i32) -> Float {
    (n as Float * EPSILON) / (1.0 - n as Float * EPSILON)
}

/// Float that tracks its error due to floating point inaccuracies
pub struct EFloat {
    value: Float,
    err: Float,

    #[cfg(debug_assertions)]
    ld: f64,
}

impl EFloat {
    /// Create a new Error float with 0 error
    pub fn new(v: Float) -> Self {
        Self {
            value: v,
            err: 0.0,
            #[cfg(debug_assertions)]
            ld: v as _,
        }
    }
}

impl Add for EFloat {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            value: self.value + rhs.value,
            err: self.err
                + rhs.err
                + gamma(1) * ((self.value + rhs.value).abs() + self.err + rhs.err),

            #[cfg(debug_assertions)]
            ld: self.ld + rhs.ld,
        }
    }
}

impl Sub for EFloat {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            value: self.value - rhs.value,
            err: self.err
                + rhs.err
                + gamma(1) * ((self.value - rhs.value).abs() + self.err + rhs.err),

            #[cfg(debug_assertions)]
            ld: self.ld - rhs.ld,
        }
    }
}

impl Mul for EFloat {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self {
            value: self.value * rhs.value,
            err: self.err
                + rhs.err
                + gamma(1) * ((self.value * rhs.value).abs() + self.err + rhs.err),

            #[cfg(debug_assertions)]
            ld: self.ld * rhs.ld,
        }
    }
}

impl Div for EFloat {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        Self {
            value: self.value / rhs.value,
            err: self.err
                + rhs.err
                + gamma(1) * ((self.value / rhs.value).abs() + self.err + rhs.err),

            #[cfg(debug_assertions)]
            ld: self.ld / rhs.ld,
        }
    }
}
