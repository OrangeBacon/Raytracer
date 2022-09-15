use std::ops::{Deref, DerefMut};

use geometry::{ConstZero, Number, Point2, Point2i, Rng};

use crate::{
    sampler::{set_sample_number, start_next_sample, start_pixel, SamplerData},
    Sampler,
};

/// Generic sampler that generates all its data at the same time
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Default, Hash)]
pub struct PixelSampler<T: Number, U: PixelSamplerImpl<T>> {
    data: PixelSamplerData<T>,
    implementation: U,
}

/// Inner data accessible by pixel sampler implementations
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Default, Hash)]
pub struct PixelSamplerData<T: Number> {
    pub data: SamplerData<T>,
    pub samples_one: Vec<Vec<T>>,
    pub samples_two: Vec<Vec<Point2<T>>>,
    pub current_one_dimension: usize,
    pub current_two_dimension: usize,
    pub rng: Rng,
}

/// Required functions for a pixel sampler implementation
pub trait PixelSamplerImpl<T: Number>: Copy {
    fn start_pixel(&mut self, data: &mut PixelSamplerData<T>, point: Point2i);
}

impl<T: Number, U: PixelSamplerImpl<T>> PixelSampler<T, U> {
    /// Create a new pixel sampler
    pub fn new(samples_per_pixel: usize, sampled_dimensions: usize, sampler: U) -> Self {
        let mut samples_one = vec![];
        let mut samples_two = vec![];

        for _ in 0..sampled_dimensions {
            samples_one.push(vec![T::ZERO; samples_per_pixel]);
            samples_two.push(vec![Point2::ZERO; samples_per_pixel]);
        }

        Self {
            data: PixelSamplerData {
                data: SamplerData::new(samples_per_pixel),
                samples_one,
                samples_two,
                current_one_dimension: 0,
                current_two_dimension: 0,
                rng: Rng::default(),
            },
            implementation: sampler,
        }
    }
}

impl<T: Number, U: PixelSamplerImpl<T>> Sampler<T> for PixelSampler<T, U> {
    fn one(&mut self) -> T {
        if self.current_one_dimension < self.samples_one.len() {
            self.current_one_dimension += 1;
            self.samples_one[self.current_one_dimension - 1][self.current_pixel_sample_index]
        } else {
            self.rng.uniform_float()
        }
    }

    fn two(&mut self) -> Point2<T> {
        if self.current_two_dimension < self.samples_two.len() {
            self.current_two_dimension += 1;
            self.samples_two[self.current_two_dimension - 1][self.current_pixel_sample_index]
        } else {
            Point2::new(self.rng.uniform_float(), self.rng.uniform_float())
        }
    }

    fn start_pixel(&mut self, point: Point2i) {
        self.implementation.start_pixel(&mut self.data, point);
        start_pixel(self, point);
    }

    fn start_next_sample(&mut self) -> bool {
        self.current_one_dimension = 0;
        self.current_two_dimension = 0;
        start_next_sample(self)
    }

    fn set_sample_number(&mut self, sample_num: usize) -> bool {
        self.current_one_dimension = 0;
        self.current_two_dimension = 0;
        set_sample_number(self, sample_num)
    }

    fn clone_seed(&mut self, seed: u64) -> Self {
        let mut new = self.clone();
        new.rng.set_sequence(seed);
        new
    }

    fn sample_data(&mut self) -> &mut SamplerData<T> {
        &mut self.data.data
    }
}

impl<T: Number, U: PixelSamplerImpl<T>> Deref for PixelSampler<T, U> {
    type Target = PixelSamplerData<T>;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl<T: Number, U: PixelSamplerImpl<T>> DerefMut for PixelSampler<T, U> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}

impl<T: Number> Deref for PixelSamplerData<T> {
    type Target = SamplerData<T>;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl<T: Number> DerefMut for PixelSamplerData<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}
