use std::ops::{
    Add, AddAssign, Div, DivAssign, Index, IndexMut, Mul, MulAssign, Neg, Sub, SubAssign,
};

use geometry::Number;

/// Base implementation of a colour spectrum
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CoefficientSpectrum<const N: usize, T: Number> {
    pub samples: [T; N],
}

impl<const N: usize, T: Number> CoefficientSpectrum<N, T> {
    pub const SAMPLE_COUNT: usize = N;
    /// Create a new spectrum with a given value at all points
    pub fn new(base: T) -> Self {
        assert!(!base.is_nan());
        Self { samples: [base; N] }
    }

    /// Is this spectrum zero everywhere?
    pub fn is_black(&self) -> bool {
        self.samples == [T::ZERO; N]
    }

    /// Take the square root of the given spectrum
    pub fn sqrt(mut self) -> Self {
        for lhs in &mut self.samples {
            *lhs = lhs.sqrt();
        }
        assert!(!self.has_nan());
        self
    }

    /// Raise the current spectrum to the power of a given value
    pub fn pow(mut self, num: T) -> Self {
        for lhs in &mut self.samples {
            // only bothering to qualify this due to rust-analyser being weird
            *lhs = Number::pow(lhs, num);
        }
        assert!(!self.has_nan());
        self
    }

    /// Calculate e ^ the current spectrum
    pub fn exp(mut self) -> Self {
        for lhs in &mut self.samples {
            *lhs = lhs.exp();
        }
        assert!(!self.has_nan());
        self
    }

    /// Linear interpolate between two spectrums
    pub fn lerp(self, other: Self, t: T) -> Self {
        self * (T::ONE - t) + other * t
    }

    /// Clamp the values in the spectrum between two values
    pub fn clamp(mut self, low: T, high: T) -> Self {
        for lhs in &mut self.samples {
            // only bothering to qualify this due to rust-analyser being weird
            *lhs = Number::clamp(lhs, low, high);
        }
        self
    }

    /// Are any of the values stored NaN?
    pub fn has_nan(&self) -> bool {
        self.samples.iter().any(|a| a.is_nan())
    }

    /// Cast all the numbers in the spectrum to another number type
    pub fn cast<U: Number>(&self) -> CoefficientSpectrum<N, U> {
        CoefficientSpectrum {
            samples: self.samples.map(|x|U::cast(x))
        }
    }
}

/// How is the light from this spectrum generated
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum SpectrumType {
    Reflectance,
    Illuminant,
}

/// Convert an XYZ colour to an RGB colour
pub fn xyz_to_rgb<T: Number>([x, y, z]: [T; 3]) -> [T; 3] {
    [
        T::cast(3.240479_f64) * x - T::cast(1.537150_f64) * y - T::cast(0.498535_f64) * z,
        T::cast(-0.969256_f64) * x + T::cast(1.875991_f64) * y + T::cast(0.041556_f64) * z,
        T::cast(0.055648_f64) * x - T::cast(0.204043_f64) * y + T::cast(1.057311_f64) * z,
    ]
}

/// Convert an RGB colour to an XYZ colour
pub fn rgb_to_xyz<T: Number>([r, g, b]: [T; 3]) -> [T; 3] {
    [
        T::cast(0.412453_f64) * r + T::cast(0.357580_f64) * g + T::cast(0.180423_f64) * b,
        T::cast(0.212671_f64) * r + T::cast(0.715160_f64) * g + T::cast(0.072169_f64) * b,
        T::cast(0.019334_f64) * r + T::cast(0.119193_f64) * g + T::cast(0.950227_f64) * b,
    ]
}

impl<const N: usize, T: Number> Default for CoefficientSpectrum<N, T> {
    fn default() -> Self {
        Self {
            samples: [T::ZERO; N],
        }
    }
}

impl<const N: usize, T: Number> Add<CoefficientSpectrum<N, T>> for CoefficientSpectrum<N, T> {
    type Output = Self;

    fn add(mut self, rhs: CoefficientSpectrum<N, T>) -> Self::Output {
        for (lhs, rhs) in self.samples.iter_mut().zip(rhs.samples) {
            *lhs += rhs;
        }
        self
    }
}

impl<const N: usize, T: Number> AddAssign<CoefficientSpectrum<N, T>> for CoefficientSpectrum<N, T> {
    fn add_assign(&mut self, rhs: CoefficientSpectrum<N, T>) {
        *self = *self + rhs;
    }
}

impl<const N: usize, T: Number> Add<T> for CoefficientSpectrum<N, T> {
    type Output = Self;

    fn add(mut self, rhs: T) -> Self::Output {
        for lhs in &mut self.samples {
            *lhs += rhs;
        }
        self
    }
}

impl<const N: usize, T: Number> AddAssign<T> for CoefficientSpectrum<N, T> {
    fn add_assign(&mut self, rhs: T) {
        *self = *self + rhs;
    }
}

impl<const N: usize, T: Number> Sub<CoefficientSpectrum<N, T>> for CoefficientSpectrum<N, T> {
    type Output = Self;

    fn sub(mut self, rhs: CoefficientSpectrum<N, T>) -> Self::Output {
        for (lhs, rhs) in self.samples.iter_mut().zip(rhs.samples) {
            *lhs -= rhs;
        }
        self
    }
}

impl<const N: usize, T: Number> SubAssign<CoefficientSpectrum<N, T>> for CoefficientSpectrum<N, T> {
    fn sub_assign(&mut self, rhs: CoefficientSpectrum<N, T>) {
        *self = *self - rhs;
    }
}

impl<const N: usize, T: Number> Sub<T> for CoefficientSpectrum<N, T> {
    type Output = Self;

    fn sub(mut self, rhs: T) -> Self::Output {
        for lhs in &mut self.samples {
            *lhs -= rhs;
        }
        self
    }
}

impl<const N: usize, T: Number> SubAssign<T> for CoefficientSpectrum<N, T> {
    fn sub_assign(&mut self, rhs: T) {
        *self = *self - rhs;
    }
}

impl<const N: usize, T: Number> Mul<CoefficientSpectrum<N, T>> for CoefficientSpectrum<N, T> {
    type Output = Self;

    fn mul(mut self, rhs: CoefficientSpectrum<N, T>) -> Self::Output {
        for (lhs, rhs) in self.samples.iter_mut().zip(rhs.samples) {
            *lhs *= rhs;
        }
        self
    }
}

impl<const N: usize, T: Number> MulAssign<CoefficientSpectrum<N, T>> for CoefficientSpectrum<N, T> {
    fn mul_assign(&mut self, rhs: CoefficientSpectrum<N, T>) {
        *self = *self * rhs;
    }
}

impl<const N: usize, T: Number> Mul<T> for CoefficientSpectrum<N, T> {
    type Output = Self;

    fn mul(mut self, rhs: T) -> Self::Output {
        for lhs in &mut self.samples {
            *lhs *= rhs;
        }
        self
    }
}

impl<const N: usize, T: Number> MulAssign<T> for CoefficientSpectrum<N, T> {
    fn mul_assign(&mut self, rhs: T) {
        *self = *self * rhs;
    }
}

impl<const N: usize, T: Number> Div<CoefficientSpectrum<N, T>> for CoefficientSpectrum<N, T> {
    type Output = Self;

    fn div(mut self, rhs: CoefficientSpectrum<N, T>) -> Self::Output {
        for (lhs, rhs) in self.samples.iter_mut().zip(rhs.samples) {
            *lhs /= rhs;
        }
        assert!(!self.has_nan());
        self
    }
}

impl<const N: usize, T: Number> DivAssign<CoefficientSpectrum<N, T>> for CoefficientSpectrum<N, T> {
    fn div_assign(&mut self, rhs: CoefficientSpectrum<N, T>) {
        *self = *self / rhs;
        assert!(!self.has_nan());
    }
}

impl<const N: usize, T: Number> Div<T> for CoefficientSpectrum<N, T> {
    type Output = Self;

    fn div(mut self, rhs: T) -> Self::Output {
        for lhs in &mut self.samples {
            *lhs /= rhs;
        }
        assert!(!self.has_nan());
        self
    }
}

impl<const N: usize, T: Number> DivAssign<T> for CoefficientSpectrum<N, T> {
    fn div_assign(&mut self, rhs: T) {
        *self = *self / rhs;
        assert!(!self.has_nan());
    }
}

impl<const N: usize, T: Number> Neg for CoefficientSpectrum<N, T> {
    type Output = Self;

    fn neg(mut self) -> Self::Output {
        for lhs in &mut self.samples {
            *lhs = -*lhs;
        }
        self
    }
}

impl<const N: usize, T: Number> Index<usize> for CoefficientSpectrum<N, T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.samples[index]
    }
}

impl<const N: usize, T: Number> IndexMut<usize> for CoefficientSpectrum<N, T> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.samples[index]
    }
}
