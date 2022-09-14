// use std::ops::{Deref, DerefMut};

// use geometry::{ConstZero, Number, Point2};

// use crate::{
//     sampler::{set_sample_number, start_next_sample, SamplerData},
//     Sampler,
// };

// /// Generic sampler that generates all its data at the same time
// pub struct PixelSampler<T: Number> {
//     data: SamplerData<T>,
//     samples_one: Vec<Vec<T>>,
//     samples_two: Vec<Vec<Point2<T>>>,
//     current_one_dimension: usize,
//     current_two_dimension: usize,
// }

// impl<T: Number> PixelSampler<T> {
//     pub fn new(samples_per_pixel: usize, sampled_dimensions: usize) -> Self {
//         let mut samples_one = vec![];
//         let mut samples_two = vec![];

//         for _ in 0..sampled_dimensions {
//             samples_one.push(vec![T::ZERO; samples_per_pixel]);
//             samples_two.push(vec![Point2::ZERO; samples_per_pixel]);
//         }

//         Self {
//             data: SamplerData::new(samples_per_pixel),
//             samples_one,
//             samples_two,
//             current_one_dimension: 0,
//             current_two_dimension: 0,
//         }
//     }
// }

// impl<T: Number> Sampler<T> for PixelSampler<T> {
//     fn one(&mut self) -> T {
//         if self.current_one_dimension < self.samples_one.len() {
//             self.current_one_dimension += 1;
//             self.samples_one[self.current_one_dimension - 1][self.current_pixel_sample_index]
//         } else {
//             todo!()
//         }
//     }

//     fn two(&mut self) -> Point2<T> {
//         if self.current_two_dimension < self.samples_two.len() {
//             self.current_two_dimension += 1;
//             self.samples_two[self.current_two_dimension - 1][self.current_pixel_sample_index]
//         } else {
//             todo!()
//         }
//     }

//     fn start_next_sample(&mut self) -> bool {
//         self.current_one_dimension = 0;
//         self.current_two_dimension = 0;
//         start_next_sample(self)
//     }

//     fn set_sample_number(&mut self, sample_num: usize) -> bool {
//         self.current_one_dimension = 0;
//         self.current_two_dimension = 0;
//         set_sample_number(self, sample_num)
//     }

//     fn clone_seed(&mut self, seed: i32) -> Self {
//         todo!()
//     }
// }

// impl<T: Number> Deref for PixelSampler<T> {
//     type Target = SamplerData<T>;

//     fn deref(&self) -> &Self::Target {
//         &self.data
//     }
// }

// impl<T: Number> DerefMut for PixelSampler<T> {
//     fn deref_mut(&mut self) -> &mut Self::Target {
//         &mut self.data
//     }
// }
