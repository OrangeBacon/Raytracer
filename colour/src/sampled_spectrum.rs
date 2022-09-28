use geometry::{lerp, Number};

use crate::CoefficientSpectrum;

/// uniformly sampled spectrum
pub type SampledSpectrum<T> = CoefficientSpectrum<60, T>;

impl<T: Number> SampledSpectrum<T> {
    /// Starting wavelength
    pub const START: usize = 400;

    /// Ending wavelength
    pub const END: usize = 700;

    /// Convert a list of (wavelength, value) samples into a uniform spectrum
    pub fn from_sampled(data: &[(T, T)]) -> Self {
        let mut data = data.to_vec();
        data.sort_by(|a, b| a.partial_cmp(b).unwrap());

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
            spectrum[idx] = Self::average_spectrum_samples(&data, lambda0, lambda1);
        }

        spectrum
    }

    /// Sample a spectrum between a given wavelengths to calculate the average value
    fn average_spectrum_samples(data: &[(T, T)], lambda_start: T, lambda_end: T) -> T {
        let n = data.len();

        if lambda_end <= data[0].0 {
            return data[0].1;
        }
        if lambda_start >= data[n - 1].0 {
            return data[n - 1].1;
        }
        if n == 1 {
            return data[0].1;
        }

        let mut sum = T::ZERO;
        if lambda_start < data[0].0 {
            sum += data[0].1 * (data[0].0 - lambda_start);
        }
        if lambda_end > data[n - 1].0 {
            sum += data[n - 1].0 * (lambda_end - data[n - 1].0);
        }

        let mut i = 0;
        while lambda_start > data[i + 1].0 {
            i += 1;
        }

        let interp = |w: T, i: usize| {
            lerp(
                (w - data[i].0) / (data[i + 1].0 - data[i].0),
                data[i].1,
                data[i + 1].1,
            )
        };

        loop {
            let start = lambda_start.max(data[i].0);
            let end = lambda_end.min(data[i + 1].0);

            sum += T::HALF * (interp(start, i) + interp(end, i)) * (end - start);

            if i + 1 >= n || lambda_end < data[i].0 {
                break;
            }
            i += 1;
        }

        sum / (lambda_end - lambda_start)
    }
}
