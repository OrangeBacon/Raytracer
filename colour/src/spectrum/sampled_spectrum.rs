use data::{cie::*, rgb_spectra::*};
use geometry::{lerp, Number};
use once_cell::sync::Lazy;

use super::{coefficient_spectrum::xyz_to_rgb, CoefficientSpectrum, RGBSpectrum, SpectrumType};

/// uniformly sampled spectrum
pub type SampledSpectrum<T> = CoefficientSpectrum<60, T>;

impl<T: Number> SampledSpectrum<T> {
    /// Starting wavelength
    pub const START: usize = 400;

    /// Ending wavelength
    pub const END: usize = 700;

    /// Convert a list of (wavelength, value) samples into a uniform spectrum
    pub fn from_sampled(lambda: &[T], value: &[T]) -> Self {
        let mut data: Vec<_> = lambda.iter().copied().zip(value.iter().copied()).collect();
        data.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
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

    /// Create a new sampled spectrum from
    pub fn from_rgb([r, g, b]: [T; 3], kind: SpectrumType) -> Self {
        let mut result = Self::default();

        let [white, cyan, yellow, magenta, red, green, blue] =
            STATIC_SPECTRA.spectra(kind == SpectrumType::Reflectance);

        if r <= g && r <= b {
            result += white * r;
            if g <= b {
                result += cyan * (g - r);
                result += blue * (b - g);
            } else {
                result += cyan * (b - r);
                result += green * (g - b);
            }
        } else if g <= r && g <= b {
            result += white * g;
            if r <= b {
                result += magenta * (r - g);
                result += blue * (b - r);
            } else {
                result += magenta * (b - g);
                result += red * (r - b);
            }
        } else {
            result += white * b;
            if r <= g {
                result += yellow * (r - b);
                result += green * (g - r);
            } else {
                result += yellow * (g - b);
                result += red * (r - g);
            }
        }

        if kind == SpectrumType::Reflectance {
            result *= T::cast(0.94);
        } else {
            result += T::cast(0.86445);
        }

        result.clamp(T::ZERO, T::INFINITY)
    }

    /// XYZ color components into a spectrum
    pub fn from_xyz(xyz: [T; 3], kind: SpectrumType) -> Self {
        let rgb = xyz_to_rgb(xyz);
        Self::from_rgb(rgb, kind)
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
        let mut result = [T::ZERO; 3];
        for i in 0..Self::SAMPLE_COUNT {
            result[0] += T::cast(STATIC_SPECTRA.x.samples[i]) * self.samples[i];
            result[1] += T::cast(STATIC_SPECTRA.y.samples[i]) * self.samples[i];
            result[2] += T::cast(STATIC_SPECTRA.z.samples[i]) * self.samples[i];
        }

        let scale = T::cast((Self::END - Self::START) as i64)
            / (T::cast(CIE_Y_INTEGRAL) * T::cast(Self::SAMPLE_COUNT as i64));

        result.map(|a| a * scale)
    }

    /// Calculate the XYZ luminance of the spectrum
    pub fn y(&self) -> T {
        let mut yy = T::ZERO;
        for i in 0..Self::SAMPLE_COUNT {
            yy += T::cast(STATIC_SPECTRA.y.samples[i]) * self.samples[i];
        }

        yy * T::cast((Self::END - Self::START) as i64) / T::cast(Self::SAMPLE_COUNT as i64)
    }

    /// Convert this spectrum into an RGB colour
    pub fn to_rgb(&self) -> [T; 3] {
        xyz_to_rgb(self.to_xyz())
    }

    /// Convert this spectrum into an RGB spectrum
    pub fn to_rgb_spectrum(&self) -> RGBSpectrum<T> {
        let rgb = self.to_rgb();
        RGBSpectrum::from_rgb(rgb, SpectrumType::Reflectance)
    }
}

#[derive(Default)]
struct StaticSpectra {
    x: SampledSpectrum<f64>,
    y: SampledSpectrum<f64>,
    z: SampledSpectrum<f64>,
    rgb_reflection_to_spectrum_white: SampledSpectrum<f64>,
    rgb_reflection_to_spectrum_cyan: SampledSpectrum<f64>,
    rgb_reflection_to_spectrum_magenta: SampledSpectrum<f64>,
    rgb_reflection_to_spectrum_yellow: SampledSpectrum<f64>,
    rgb_reflection_to_spectrum_red: SampledSpectrum<f64>,
    rgb_reflection_to_spectrum_green: SampledSpectrum<f64>,
    rgb_reflection_to_spectrum_blue: SampledSpectrum<f64>,
    rgb_illumination_to_spectrum_white: SampledSpectrum<f64>,
    rgb_illumination_to_spectrum_cyan: SampledSpectrum<f64>,
    rgb_illumination_to_spectrum_magenta: SampledSpectrum<f64>,
    rgb_illumination_to_spectrum_yellow: SampledSpectrum<f64>,
    rgb_illumination_to_spectrum_red: SampledSpectrum<f64>,
    rgb_illumination_to_spectrum_green: SampledSpectrum<f64>,
    rgb_illumination_to_spectrum_blue: SampledSpectrum<f64>,
}

static STATIC_SPECTRA: Lazy<StaticSpectra> = Lazy::new(StaticSpectra::new);

impl StaticSpectra {
    fn new() -> Self {
        let mut spectrums = StaticSpectra::default();

        type SS = SampledSpectrum<f64>;

        for i in 0..SS::SAMPLE_COUNT {
            let wl0 = lerp((i / SS::SAMPLE_COUNT) as _, SS::START as _, SS::END as _);
            let wl1 = lerp(
                ((i + 1) / SS::SAMPLE_COUNT) as _,
                SS::START as _,
                SS::END as _,
            );

            let s = |a, b| SampledSpectrum::average_spectrum_samples(a, b, wl0, wl1);

            spectrums.x.samples[i] = s(CIE_LAMBDA, CIE_X);
            spectrums.y.samples[i] = s(CIE_LAMBDA, CIE_Y);
            spectrums.z.samples[i] = s(CIE_LAMBDA, CIE_Z);

            spectrums.rgb_reflection_to_spectrum_white.samples[i] =
                s(RGB_TO_SPECTRUM_LAMBDA, RGB_REFLECTION_TO_SPECTRUM_WHITE);
            spectrums.rgb_reflection_to_spectrum_cyan.samples[i] =
                s(RGB_TO_SPECTRUM_LAMBDA, RGB_REFLECTION_TO_SPECTRUM_CYAN);
            spectrums.rgb_reflection_to_spectrum_magenta.samples[i] =
                s(RGB_TO_SPECTRUM_LAMBDA, RGB_REFLECTION_TO_SPECTRUM_MAGENTA);
            spectrums.rgb_reflection_to_spectrum_yellow.samples[i] =
                s(RGB_TO_SPECTRUM_LAMBDA, RGB_REFLECTION_TO_SPECTRUM_YELLOW);
            spectrums.rgb_reflection_to_spectrum_red.samples[i] =
                s(RGB_TO_SPECTRUM_LAMBDA, RGB_REFLECTION_TO_SPECTRUM_RED);
            spectrums.rgb_reflection_to_spectrum_green.samples[i] =
                s(RGB_TO_SPECTRUM_LAMBDA, RGB_REFLECTION_TO_SPECTRUM_GREEN);
            spectrums.rgb_reflection_to_spectrum_blue.samples[i] =
                s(RGB_TO_SPECTRUM_LAMBDA, RGB_REFLECTION_TO_SPECTRUM_BLUE);

            spectrums.rgb_illumination_to_spectrum_white.samples[i] =
                s(RGB_TO_SPECTRUM_LAMBDA, RGB_ILLUMINATION_TO_SPECTRUM_WHITE);
            spectrums.rgb_illumination_to_spectrum_cyan.samples[i] =
                s(RGB_TO_SPECTRUM_LAMBDA, RGB_ILLUMINATION_TO_SPECTRUM_CYAN);
            spectrums.rgb_illumination_to_spectrum_magenta.samples[i] =
                s(RGB_TO_SPECTRUM_LAMBDA, RGB_ILLUMINATION_TO_SPECTRUM_MAGENTA);
            spectrums.rgb_illumination_to_spectrum_yellow.samples[i] =
                s(RGB_TO_SPECTRUM_LAMBDA, RGB_ILLUMINATION_TO_SPECTRUM_YELLOW);
            spectrums.rgb_illumination_to_spectrum_red.samples[i] =
                s(RGB_TO_SPECTRUM_LAMBDA, RGB_ILLUMINATION_TO_SPECTRUM_RED);
            spectrums.rgb_illumination_to_spectrum_green.samples[i] =
                s(RGB_TO_SPECTRUM_LAMBDA, RGB_ILLUMINATION_TO_SPECTRUM_GREEN);
            spectrums.rgb_illumination_to_spectrum_blue.samples[i] =
                s(RGB_TO_SPECTRUM_LAMBDA, RGB_ILLUMINATION_TO_SPECTRUM_BLUE);
        }

        spectrums
    }

    fn spectra<T: Number>(&self, is_reflectance: bool) -> [SampledSpectrum<T>; 7] {
        let result = if is_reflectance {
            [
                self.rgb_reflection_to_spectrum_white,
                self.rgb_reflection_to_spectrum_cyan,
                self.rgb_reflection_to_spectrum_magenta,
                self.rgb_reflection_to_spectrum_yellow,
                self.rgb_reflection_to_spectrum_red,
                self.rgb_reflection_to_spectrum_green,
                self.rgb_reflection_to_spectrum_blue,
            ]
        } else {
            [
                self.rgb_illumination_to_spectrum_white,
                self.rgb_illumination_to_spectrum_cyan,
                self.rgb_illumination_to_spectrum_magenta,
                self.rgb_illumination_to_spectrum_yellow,
                self.rgb_illumination_to_spectrum_red,
                self.rgb_illumination_to_spectrum_green,
                self.rgb_illumination_to_spectrum_blue,
            ]
        };

        result.map(|s| s.cast())
    }
}
