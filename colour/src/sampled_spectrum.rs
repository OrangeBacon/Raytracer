use data::cie::{CIE_LAMBDA, CIE_X, CIE_Y, CIE_Z};
use geometry::{lerp, Number};
use once_cell::sync::OnceCell;

use crate::CoefficientSpectrum;

/// uniformly sampled spectrum
pub type SampledSpectrum<T> = CoefficientSpectrum<60, T>;

static XYZ_CURVES64: OnceCell<[SampledSpectrum<f64>; 3]> = OnceCell::new();

impl<T: Number> SampledSpectrum<T> {
    /// Starting wavelength
    pub const START: usize = 400;

    /// Ending wavelength
    pub const END: usize = 700;

    /// Convert a list of (wavelength, value) samples into a uniform spectrum
    pub fn from_sampled(lambda: &[T], value: &[T]) -> Self {
        let mut data: Vec<_> = lambda.iter().copied().zip(value.iter().copied()).collect();
        data.sort_unstable_by(|a, b| a.partial_cmp(&b).unwrap());
        let (lambda, value): (Vec<_>, Vec<_>) = data.into_iter().unzip();

        let mut spectrum = CoefficientSpectrum::new(T::ZERO);
        for idx in 0..Self::SAMPLE_COUNT {
            let lambda0 = lerp(
                T::cast(idx as i32) / T::cast(Self::SAMPLE_COUNT as i32),
                T::cast(Self::START as i32),
                T::cast(Self::END as i32),
            );
            let lambda1 = lerp(
                T::cast(idx as i32 + 1) / T::cast(Self::SAMPLE_COUNT as i32),
                T::cast(Self::START as i32),
                T::cast(Self::END as i32),
            );
            spectrum[idx] = Self::average_spectrum_samples(&lambda, &value, lambda0, lambda1);
        }

        spectrum
    }

    /// Sample a spectrum between given wavelengths to calculate the average value.
    /// Data is provided in (lambda, value) tuples
    fn average_spectrum_samples(lambda: &[T], value: &[T], lambda_start: T, lambda_end: T) -> T {
        assert_eq!(lambda.len(), value.len());

        let n = lambda.len();

        if lambda_end <= lambda[0] {
            return value[0];
        }
        if lambda_start >= lambda[n - 1] {
            return value[n - 1];
        }
        if n == 1 {
            return value[0];
        }

        let mut sum = T::ZERO;
        if lambda_start < lambda[0] {
            sum += value[0] * (lambda[0] - lambda_start);
        }
        if lambda_end > lambda[n - 1] {
            sum += lambda[n - 1] * (lambda_end - lambda[n - 1]);
        }

        let mut i = 0;
        while lambda_start > lambda[i + 1] {
            i += 1;
        }

        let interp = |w: T, i: usize| {
            lerp(
                (w - lambda[i]) / (lambda[i + 1] - lambda[i]),
                value[i],
                value[i + 1],
            )
        };

        loop {
            let start = lambda_start.max(lambda[i]);
            let end = lambda_end.min(lambda[i + 1]);

            sum += T::HALF * (interp(start, i) + interp(end, i)) * (end - start);

            if i + 1 >= n || lambda_end < lambda[i] {
                break;
            }
            i += 1;
        }

        sum / (lambda_end - lambda_start)
    }

    /// Calculate the XYZ color coefficients [x, y, z] for this spectrum
    pub fn to_xyz(&self) -> [T; 3] {
        let spectrums = XYZ_CURVES64.get_or_init(initialise_xyz_spectrums);

        let mut result = [T::ZERO; 3];
        for i in 0..Self::SAMPLE_COUNT {
            result[0] += T::cast(spectrums[0].samples[i]) * self.samples[i];
            result[1] += T::cast(spectrums[1].samples[i]) * self.samples[i];
            result[2] += T::cast(spectrums[2].samples[i]) * self.samples[i];
        }

        let scale = T::cast((Self::END - Self::START) as i64)
            / (T::cast(106.856895) * T::cast(Self::SAMPLE_COUNT as i64));

        result.map(|a| a * scale)
    }

    pub fn y(&self) -> T {
        let &[_, y, _] = XYZ_CURVES64.get_or_init(initialise_xyz_spectrums);

        let mut yy = T::ZERO;
        for i in 0..Self::SAMPLE_COUNT {
            yy += T::cast(y.samples[i]) * self.samples[i];
        }

        yy * T::cast((Self::END - Self::START) as i64) / T::cast(Self::SAMPLE_COUNT as i64)
    }
}

fn initialise_xyz_spectrums() -> [SampledSpectrum<f64>; 3] {
    let mut spectrums = [SampledSpectrum::default(); 3];

    type SS = SampledSpectrum<f64>;

    for i in 0..SS::SAMPLE_COUNT {
        let wl0 = lerp((i / SS::SAMPLE_COUNT) as _, SS::START as _, SS::END as _);
        let wl1 = lerp(
            ((i + 1) / SS::SAMPLE_COUNT) as _,
            SS::START as _,
            SS::END as _,
        );

        spectrums[0].samples[i] =
            SampledSpectrum::average_spectrum_samples(CIE_LAMBDA, CIE_X, wl0, wl1);
        spectrums[1].samples[i] =
            SampledSpectrum::average_spectrum_samples(CIE_LAMBDA, CIE_Y, wl0, wl1);
        spectrums[2].samples[i] =
            SampledSpectrum::average_spectrum_samples(CIE_LAMBDA, CIE_Z, wl0, wl1);
    }

    spectrums
}
