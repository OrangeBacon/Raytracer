use geometry::{log2usize, Number, Point2};

use crate::{
    pixel_sampler::{PixelSampler, PixelSamplerImpl},
    zero_two::{sobol_2d, van_der_corput},
};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default, Hash)]
pub struct MaxMinDistSampler {
    matrix_index: usize,
}

impl MaxMinDistSampler {
    pub fn new<T: Number>(
        samples_per_pixel: usize,
        sampled_dimensions: usize,
    ) -> PixelSampler<T, MaxMinDistSampler> {
        let sampler = MaxMinDistSampler {
            matrix_index: log2usize(samples_per_pixel) as usize,
        };

        PixelSampler::new(samples_per_pixel, sampled_dimensions, sampler)
    }
}

impl<T: Number> PixelSamplerImpl<T> for MaxMinDistSampler {
    fn start_pixel(
        &mut self,
        data: &mut crate::pixel_sampler::PixelSamplerData<T>,
        _point: geometry::Point2i,
    ) {
        let inv = T::ONE / T::cast(data.samples_per_pixel as i32);

        for i in 0..data.samples_per_pixel {
            data.samples_two[0][i] = Point2::new(
                T::cast(i as i32) * inv,
                sample_generator_matrix(&MAX_MIN_MATRIX[self.matrix_index], i as _, 0),
            );
        }
        let samples = data.samples_per_pixel;

        data.rng
            .shuffle_dims(&mut data.samples_two[0][0..samples], 1);

        for i in 0..data.samples_one.len() {
            van_der_corput(1, samples, &mut data.samples_one[i], &mut data.rng);
        }

        for i in 0..data.samples_two.len() {
            sobol_2d(1, samples, &mut data.samples_two[i], &mut data.rng);
        }

        for i in 0..data.samples_one_array_sizes.len() {
            van_der_corput(
                data.samples_one_array_sizes[i],
                samples,
                &mut data.data.sample_array_one[i],
                &mut data.rng,
            );
        }

        for i in 0..data.samples_two_array_sizes.len() {
            sobol_2d(
                data.samples_two_array_sizes[i],
                samples,
                &mut data.data.sample_array_two[i],
                &mut data.rng,
            );
        }
    }
}

/// Multiply a generator matrix by a bit vector a
fn multiply_generator(lines: &[u32], mut a: u32) -> u32 {
    let mut res = 0;

    for &line in lines {
        if a & 1 != 0 {
            res ^= line
        }

        if a == 0 {
            break;
        }

        a >>= 1;
    }

    res
}

/// Generate sample value from a generator matrix
fn sample_generator_matrix<T: Number>(lines: &[u32], a: u32, scramble: u32) -> T {
    T::ONE_MINUS_EPSILON.min(T::cast(
        (multiply_generator(lines, a) ^ scramble) as f32 * 2.0f32.powi(-32),
    ))
}

const MAX_MIN_MATRIX: [[u32; 32]; 17] = [
    [
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x80000000,
    ],
    [
        0xc0000000, 0x80000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0xa0000000, 0x40000000, 0x80000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0xc0000000, 0x50000000, 0x20000000, 0x30000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0x88000000, 0x58000000, 0x20000000, 0x40000000, 0x80000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0xe0000000, 0x60000000, 0x28000000, 0x10000000, 0x18000000, 0x04000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0x82000000, 0x44000000, 0x2c000000, 0x10000000, 0x20000000, 0x40000000, 0x80000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0x80000000, 0x70000000, 0x30000000, 0x14000000, 0x08000000, 0x0c000000, 0x02000000,
        0x01000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0xc0000000, 0x41000000, 0x22000000, 0x16000000, 0x08000000, 0x10000000, 0x20000000,
        0x40800000, 0x80000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0x80000000, 0x40000000, 0x20000000, 0x18000000, 0x08000000, 0x1c000000, 0x1e000000,
        0x03000000, 0x00800000, 0x00400000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0x80000000, 0x60400000, 0x20800000, 0x11000000, 0x0b000000, 0x04000000, 0x08000000,
        0x10000000, 0x20000000, 0x40000000, 0x00200000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0x80000000, 0x40000000, 0x20000000, 0x1c000000, 0x0c000000, 0x05000000, 0x02000000,
        0x03000000, 0x00800000, 0x00400000, 0x00200000, 0x00100000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0x80000000, 0x40000000, 0x30200000, 0x10400000, 0x08800000, 0x05800000, 0x02000000,
        0x04000000, 0x08000000, 0x10000000, 0x20000000, 0x00100000, 0x00080000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0x80000000, 0x40000000, 0x20000000, 0x10000000, 0x0e000000, 0x06000000, 0x02800000,
        0x01000000, 0x01800000, 0x00400000, 0x00200000, 0x00100000, 0x00080000, 0x00040000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0x80000000, 0x40000000, 0x20000000, 0x18100000, 0x08200000, 0x04400000, 0x02c00000,
        0x01000000, 0x02000000, 0x04000000, 0x08000000, 0x10000000, 0x00080000, 0x00040000,
        0x00020000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
    [
        0x80000000, 0x40000000, 0x20000000, 0x10000000, 0x08000000, 0x07000000, 0x03000000,
        0x01400000, 0x00800000, 0x00c00000, 0x00200000, 0x00100000, 0x00080000, 0x00040000,
        0x00020000, 0x00010000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000,
        0x00000000, 0x00000000, 0x00000000, 0x00000000,
    ],
];
