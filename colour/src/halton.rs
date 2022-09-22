use data::{RADICAL_INVERSE_PERMUTATIONS, PRIME_SUMS};
use geometry::{Bounds2i, Number, Point2i};

use crate::{
    global_sampler::{GlobalSampler, GlobalSamplerData, GlobalSamplerImpl},
    low_discrepancy::{inverse_radical_inverse, radical_inverse, scrambled_radical_inverse},
};

/// Maximum number of points in one direction in one dimension to generate so that
/// floating point precision is kept.
const MAX_RESOLUTION: i32 = 128;

/// Sampler using the halton sequence
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default, Hash)]
pub struct HaltonSampler {
    base_scales: Point2i,
    base_exponents: Point2i,
    sample_stride: usize,
    multiplicative_inverse: [i64; 2],
    pixel_for_offset: Point2i,
    offset_for_current_pixel: usize,
    sample_at_pixel_centre: bool,
}

impl HaltonSampler {
    /// Create a new halton sampler
    pub fn new<T: Number>(
        samples_per_pixel: usize,
        sample_bounds: Bounds2i,
    ) -> GlobalSampler<T, HaltonSampler> {
        let mut this = HaltonSampler::default();

        let res = sample_bounds.max - sample_bounds.min;
        for i in 0..2 {
            let base = i + 2; // (i == 0) ? 2 : 3
            let mut scale = 1;
            let mut exp = 0;
            while scale < res[i].min(MAX_RESOLUTION) {
                scale *= base as i32;
                exp += 1;
            }

            this.base_scales[i] = scale;
            this.base_exponents[i] = exp;
        }

        this.sample_stride = (this.base_scales.x * this.base_scales.y) as _;
        this.multiplicative_inverse[0] =
            multiplicative_inverse(this.base_scales[1] as _, this.base_scales[0] as _);
        this.multiplicative_inverse[1] =
            multiplicative_inverse(this.base_scales[0] as _, this.base_scales[1] as _);
        this.pixel_for_offset = Point2i::new(i32::MAX, i32::MAX);

        GlobalSampler::new(samples_per_pixel, this)
    }
}

impl<T: Number> GlobalSamplerImpl<T> for HaltonSampler {
    const ARRAY_START_DIM: usize = 5;

    fn get_index_for_sample(&mut self, data: &GlobalSamplerData<T>, sample_num: usize) -> usize {
        if data.current_pixel != self.pixel_for_offset {
            self.offset_for_current_pixel = 0;
            if self.sample_stride > 1 {
                let pm = Point2i::new(
                    data.current_pixel[0].rem_euclid(MAX_RESOLUTION),
                    data.current_pixel[1].rem_euclid(MAX_RESOLUTION),
                );

                for i in 0..2 {
                    let dim_offset = if i == 0 {
                        inverse_radical_inverse::<2>(pm[i] as _, self.base_exponents[i] as _)
                    } else {
                        inverse_radical_inverse::<3>(pm[i] as _, self.base_exponents[i] as _)
                    };

                    self.offset_for_current_pixel += dim_offset as usize
                        * (self.sample_stride / self.base_scales[i] as usize)
                        * self.multiplicative_inverse[i] as usize;
                }

                self.offset_for_current_pixel %= self.sample_stride;
            }

            self.pixel_for_offset = data.current_pixel;
        }

        self.offset_for_current_pixel + sample_num * self.sample_stride
    }

    fn sample_dimension(&mut self, _: &GlobalSamplerData<T>, index: usize, dimension: usize) -> T {
        if self.sample_at_pixel_centre && (dimension == 0 || dimension == 1) {
            return T::HALF;
        }

        if dimension == 0 {
            radical_inverse(dimension, (index >> self.base_exponents[0]) as _)
        } else if dimension == 1 {
            radical_inverse(dimension, (index >> self.base_scales[1]) as _)
        } else {
            scrambled_radical_inverse(
                dimension,
                index as _,
                &RADICAL_INVERSE_PERMUTATIONS[PRIME_SUMS[dimension]..],
            )
        }
    }

    fn clone_seed(&mut self, _: u64) -> Self {
        self.clone()
    }
}

fn multiplicative_inverse(a: i64, n: i64) -> i64 {
    let (x, _) = extended_gcd(a, n);
    x.rem_euclid(n)
}

fn extended_gcd(a: i64, b: i64) -> (i64, i64) {
    if b == 0 {
        return (1, 0);
    }

    let (xp, yp) = extended_gcd(b, a % b);
    (yp, xp - (a / b * yp))
}
