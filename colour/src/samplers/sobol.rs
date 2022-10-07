use geometry::{log2usize, round_up_pow_2, Bounds2i, Number};

use crate::samplers::global_sampler::{GlobalSampler, GlobalSamplerData, GlobalSamplerImpl};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default, Hash)]
pub struct Sobol {
    bounds: Bounds2i,
    resolution: usize,
    log2resolution: usize,
}

impl Sobol {
    /// Create a new sobol sampler
    pub fn new<T: Number>(samples_per_pixel: usize, bounds: Bounds2i) -> GlobalSampler<T, Sobol> {
        let diagonal = bounds.diagonal();
        let resolution = round_up_pow_2(diagonal.x.max(diagonal.y) as _);

        let sampler = Sobol {
            bounds,
            resolution,
            log2resolution: log2usize(resolution),
        };

        GlobalSampler::new(round_up_pow_2(samples_per_pixel), sampler)
    }
}

impl<T: Number> GlobalSamplerImpl<T> for Sobol {
    const ARRAY_START_DIM: usize = 5;

    fn get_index_for_sample(
        &mut self,
        data: &GlobalSamplerData<T>,
        mut sample_num: usize,
    ) -> usize {
        let point = (data.current_pixel - self.bounds.min).to_point();

        if self.log2resolution == 0 {
            return 0;
        }

        let m2 = self.log2resolution << 1;
        let mut index = sample_num << m2;

        let mut delta = 0;
        let mut c = 0;
        loop {
            if sample_num & 1 != 0 {
                delta ^= data::sobol::MATRICES[self.log2resolution - 1][c];
            }
            if sample_num == 0 {
                break;
            }
            sample_num >>= 1;
            c += 1;
        }

        let mut b =
            ((point.x as u32) << self.log2resolution) as u64 | ((point.y as u32 as u64) ^ delta);

        let mut c = 0;
        loop {
            if b & 1 != 0 {
                index ^= data::sobol::MATRICES_INV[self.log2resolution - 1][c] as usize;
            }
            if b == 0 {
                break;
            }
            b >>= 1;
            c += 1;
        }

        index
    }

    fn sample_dimension(
        &mut self,
        data: &GlobalSamplerData<T>,
        index: usize,
        dimension: usize,
    ) -> T {
        if dimension >= data::sobol::DIMENSIONS {
            panic!("Run out of sobol sampler dimensions")
        }

        let mut s: T = sobol_sample(index as _, dimension, 0);

        if dimension == 0 || dimension == 1 {
            s = s * T::cast(self.resolution as i32) + T::cast(self.bounds.min[dimension]);
            s -= T::cast(data.current_pixel[dimension]);
            s = s.clamp(T::ZERO, T::ONE_MINUS_EPSILON);
        }

        s
    }

    fn clone_seed(&mut self, _: u64) -> Self {
        *self
    }
}

/// Use more accurate implementation of sobol sampler for f64 if relevant
fn sobol_sample<T: Number>(index: u64, dimension: usize, scramble: u64) -> T {
    if std::mem::size_of::<T>() == std::mem::size_of::<f32>() {
        T::cast(sobol_sample_f32(index, dimension, scramble as _))
    } else {
        T::cast(sobol_sample_f64(index, dimension, scramble))
    }
}

fn sobol_sample_f32(mut index: u64, dimension: usize, mut scramble: u32) -> f32 {
    assert!(
        dimension < data::sobol::DIMENSIONS,
        "Too many sobol dimensions consumed"
    );

    for i in dimension * data::sobol::MATRIX_SIZE.. {
        if index & 1 != 0 {
            scramble ^= data::sobol::MATRICES_32[i]
        }

        index >>= 1;
        if index == 0 {
            break;
        }
    }

    (scramble as f32 * 2.0f32.powi(-32)).min(f32::ONE_MINUS_EPSILON)
}

fn sobol_sample_f64(mut index: u64, dimension: usize, scramble: u64) -> f64 {
    assert!(
        dimension < data::sobol::DIMENSIONS,
        "Too many sobol dimensions consumed"
    );

    let mut result = scramble & (!-(1i64 << data::sobol::MATRIX_SIZE)) as u64;

    for i in dimension * data::sobol::MATRIX_SIZE.. {
        if index & 1 != 0 {
            result ^= data::sobol::MATRICES_64[i];
        }
        index >>= 1;
        if index == 0 {
            break;
        }
    }

    ((result as f64) * (1.0 / ((1i64 << data::sobol::MATRIX_SIZE) as f64)))
        .min(f64::ONE_MINUS_EPSILON)
}
