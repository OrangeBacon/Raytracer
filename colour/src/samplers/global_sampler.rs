use std::ops::{Deref, DerefMut};

use geometry::{Number, Point2};

use crate::samplers::sampler::{
    set_sample_number, start_next_sample, start_pixel, Sampler, SamplerData,
};

/// A sampler that generates samples for multiple pixels at once
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Default, Hash)]
pub struct GlobalSampler<T: Number, U: GlobalSamplerImpl<T>> {
    data: GlobalSamplerData<T>,
    sampler: U,
}

/// Data required to implement a global sampler
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Default, Hash)]
pub struct GlobalSamplerData<T: Number> {
    data: SamplerData<T>,

    /// The next dimension the sampler will be asked for
    dimension: usize,

    /// Index of the sample in the current pixel
    interval_sample_index: usize,

    /// dimension at the end of the samples
    array_end_dim: usize,
}

/// Functions required to implement a global sampler
pub trait GlobalSamplerImpl<T: Number>: Copy {
    /// The number of dimensions that are skipped before generating array samples
    const ARRAY_START_DIM: usize = 5;

    /// Map the pixel-local sample index into a global sample index
    fn get_index_for_sample(&mut self, data: &GlobalSamplerData<T>, sample_num: usize) -> usize;

    /// Return the sample at the given dimension and global sample index
    fn sample_dimension(
        &mut self,
        data: &GlobalSamplerData<T>,
        index: usize,
        dimension: usize,
    ) -> T;

    /// Clone this sampler, but set the random seed to the given seed
    fn clone_seed(&mut self, seed: u64) -> Self;
}

impl<T: Number, U: GlobalSamplerImpl<T>> GlobalSampler<T, U> {
    /// Create a new global sampler data instance
    pub fn new(samples_per_pixel: usize, sampler: U) -> Self {
        GlobalSampler {
            data: GlobalSamplerData {
                data: SamplerData::new(samples_per_pixel),
                ..Default::default()
            },
            sampler,
        }
    }
}

impl<T: Number, U: GlobalSamplerImpl<T>> Sampler<T> for GlobalSampler<T, U> {
    fn start_pixel(&mut self, point: Point2<i32>) {
        start_pixel(self, point);
        self.dimension = 0;
        self.interval_sample_index = self.sampler.get_index_for_sample(&self.data, 0);
        self.array_end_dim =
            U::ARRAY_START_DIM + self.sample_array_one.len() + 2 * self.sample_array_two.len();
        let samples_per_pixel = self.samples_per_pixel;

        for i in 0..self.samples_one_array_sizes.len() {
            let n_samples = self.samples_one_array_sizes[i] * samples_per_pixel;
            for j in 0..n_samples {
                let index = self.sampler.get_index_for_sample(&self.data, j);
                self.sample_array_one[i][j] =
                    self.sampler
                        .sample_dimension(&self.data, index, U::ARRAY_START_DIM + i);
            }
        }

        let mut dim = U::ARRAY_START_DIM + self.samples_one_array_sizes.len();
        for i in 0..self.samples_two_array_sizes.len() {
            let n_samples = self.samples_two_array_sizes[i] * samples_per_pixel;
            for j in 0..n_samples {
                let index = self.sampler.get_index_for_sample(&self.data, j);
                self.sample_array_one[i][j] = self.sampler.sample_dimension(&self.data, index, dim);
                self.sample_array_one[i][j] =
                    self.sampler.sample_dimension(&self.data, index, dim + 1);
            }
            dim += 2;
            assert_eq!(dim, self.array_end_dim);
        }
    }

    fn start_next_sample(&mut self) -> bool {
        self.dimension = 0;
        let sample_num = self.current_pixel_sample_index + 1;
        self.interval_sample_index = self.sampler.get_index_for_sample(&self.data, sample_num);
        start_next_sample(self)
    }

    fn set_sample_number(&mut self, sample_num: usize) -> bool {
        self.dimension = 0;
        self.interval_sample_index = self.sampler.get_index_for_sample(&self.data, sample_num);
        set_sample_number(self, sample_num)
    }

    fn one(&mut self) -> T {
        if self.dimension >= U::ARRAY_START_DIM && self.dimension < self.array_end_dim {
            self.dimension = self.array_end_dim;
        }

        self.dimension += 1;
        let idx = self.interval_sample_index;
        let dimension = self.dimension - 1;
        self.sampler.sample_dimension(&self.data, idx, dimension)
    }

    fn two(&mut self) -> Point2<T> {
        if self.dimension + 1 >= U::ARRAY_START_DIM && self.dimension < self.array_end_dim {
            self.dimension = self.array_end_dim;
        }

        self.dimension += 2;
        let idx = self.interval_sample_index;
        let dimension = self.dimension - 1;
        Point2::new(
            self.sampler.sample_dimension(&self.data, idx, dimension),
            self.sampler
                .sample_dimension(&self.data, idx, dimension + 1),
        )
    }

    fn clone_seed(&mut self, seed: u64) -> Self {
        Self {
            sampler: self.sampler.clone_seed(seed),
            data: self.data.clone(),
        }
    }

    fn sample_data(&mut self) -> &mut SamplerData<T> {
        &mut self.data
    }
}

impl<T: Number, U: GlobalSamplerImpl<T>> Deref for GlobalSampler<T, U> {
    type Target = GlobalSamplerData<T>;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl<T: Number, U: GlobalSamplerImpl<T>> DerefMut for GlobalSampler<T, U> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}

impl<T: Number> Deref for GlobalSamplerData<T> {
    type Target = SamplerData<T>;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl<T: Number> DerefMut for GlobalSamplerData<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}
