use data::cie::{CIE_LAMBDA, CIE_X, CIE_Y, CIE_Y_INTEGRAL, CIE_Z};
use geometry::{find_interval, lerp, Number};

use super::{
    coefficient_spectrum::{rgb_to_xyz, xyz_to_rgb},
    CoefficientSpectrum, SpectrumType,
};

/// Colour spectrum that only stores 3 RGB values
pub type RGBSpectrum<T> = CoefficientSpectrum<3, T>;

impl<T: Number> RGBSpectrum<T> {
    pub fn from_sampled(lambda: &[T], value: &[T]) -> Self {
        let mut data: Vec<_> = lambda.iter().copied().zip(value.iter().copied()).collect();
        data.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
        let (lambda, value): (Vec<_>, Vec<_>) = data.into_iter().unzip();

        let mut xyz = [T::ZERO; 3];
        for i in 0..CIE_LAMBDA.len() {
            let val = Self::interpolate_spectrum_samples(&lambda, &value, T::cast(CIE_LAMBDA[i]));
            xyz[0] += val * T::cast(CIE_X[i]);
            xyz[1] += val * T::cast(CIE_Y[i]);
            xyz[2] += val * T::cast(CIE_Z[i]);
        }

        let scale = T::cast(
            (CIE_LAMBDA[CIE_LAMBDA.len() - 1] - CIE_LAMBDA[0])
                / (CIE_Y_INTEGRAL * CIE_LAMBDA.len() as f64),
        );

        xyz[0] *= scale;
        xyz[1] *= scale;
        xyz[2] *= scale;

        Self::from_xyz(xyz)
    }

    fn interpolate_spectrum_samples(lambda: &[T], value: &[T], l: T) -> T {
        if l <= lambda[0] {
            return value[0];
        }
        if l >= lambda[lambda.len() - 1] {
            return lambda[lambda.len() - 1];
        }

        let offset = find_interval(lambda.len(), |idx| lambda[idx] <= l);

        let t = (l - lambda[offset]) / (lambda[offset + 1] - lambda[offset]);

        lerp(t, value[offset], value[offset + 1])
    }

    pub fn from_rgb(rgb: [T; 3], _kind: SpectrumType) -> Self {
        let mut new = Self::new(T::ZERO);
        new.samples = rgb;
        new
    }

    pub fn to_rgb(&self) -> [T; 3] {
        self.samples
    }

    pub fn from_xyz(xyz: [T; 3]) -> Self {
        let mut new = Self::new(T::ZERO);
        new.samples = xyz_to_rgb(xyz);
        new
    }

    pub fn to_xyz(&self) -> [T; 3] {
        rgb_to_xyz(self.samples)
    }
}
