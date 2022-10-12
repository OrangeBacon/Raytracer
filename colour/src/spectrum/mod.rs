mod coefficient_spectrum;
mod rgb_spectrum;
mod sampled_spectrum;

pub use coefficient_spectrum::CoefficientSpectrum;
pub use coefficient_spectrum::SpectrumType;

pub use rgb_spectrum::RGBSpectrum;

pub use sampled_spectrum::SampledSpectrum;

use geometry::Number;

pub trait Spectrum<T: Number> {
    // pub fn from_sampled(lambda: &[T], value: &[T]) -> Self {

    // pub fn to_xyz(&self) -> [T; 3] {

    // pub fn y(&self) -> T {

    // pub fn to_rgb(&self) -> [T; 3] {

    // pub fn from_rgb([r, g, b]: [T; 3], kind: SpectrumType) -> Self {

    // pub fn from_xyz(xyz: [T; 3], kind: SpectrumType) -> Self {

    // to_rgb_spectrum
}
