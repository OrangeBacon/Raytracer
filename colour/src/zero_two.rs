use geometry::{Number, Point2, Point2i, Rng};

use crate::pixel_sampler::{PixelSampler, PixelSamplerData, PixelSamplerImpl};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default, Hash)]
pub struct ZeroTwoSampler;

impl ZeroTwoSampler {
    /// Create a new (0, 2)-sequence sampler
    pub fn new<T: Number>(
        samples_per_pixel: usize,
        sampled_dimensions: usize,
    ) -> PixelSampler<T, ZeroTwoSampler> {
        PixelSampler::new(
            round_up_pow_2(samples_per_pixel),
            sampled_dimensions,
            ZeroTwoSampler,
        )
    }
}

impl<T: Number> PixelSamplerImpl<T> for ZeroTwoSampler {
    fn start_pixel(&mut self, data: &mut PixelSamplerData<T>, _point: Point2i) {
        // pixel samples
        for i in 0..data.samples_one.len() {
            van_der_corput(
                1,
                data.samples_per_pixel,
                &mut data.samples_one[i],
                &mut data.rng,
            );
        }
        for i in 0..data.samples_two.len() {
            sobol_2d(
                1,
                data.samples_per_pixel,
                &mut data.samples_two[i],
                &mut data.rng,
            );
        }

        // array samples
        for i in 0..data.samples_one_array_sizes.len() {
            van_der_corput(
                data.samples_one_array_sizes[i],
                data.samples_per_pixel,
                &mut data.data.sample_array_one[i],
                &mut data.rng,
            )
        }

        for i in 0..data.samples_two_array_sizes.len() {
            sobol_2d(
                data.samples_two_array_sizes[i],
                data.samples_per_pixel,
                &mut data.data.sample_array_two[i],
                &mut data.rng,
            )
        }
    }

    fn round_count(&mut self, count: usize) -> usize {
        round_up_pow_2(count)
    }
}

/// Round a number up to the nearest power of 2.  Returns the input if it is
/// a power of two.
fn round_up_pow_2(x: usize) -> usize {
    1 << (std::mem::size_of::<usize>() * 8 - (x - 1).leading_zeros() as usize)
}

/// Create scrambled 1D sample values using grey code
fn van_der_corput<T: Number>(
    samples_per_pixel: usize,
    pixel_samples: usize,
    samples: &mut [T],
    rng: &mut Rng,
) {
    const MATRIX: [u32; 32] = [
        0x80000000, 0x40000000, 0x20000000, 0x10000000, 0x8000000, 0x4000000, 0x2000000, 0x1000000,
        0x800000, 0x400000, 0x200000, 0x100000, 0x80000, 0x40000, 0x20000, 0x10000, 0x8000, 0x4000,
        0x2000, 0x1000, 0x800, 0x400, 0x200, 0x100, 0x80, 0x40, 0x20, 0x10, 0x8, 0x4, 0x2, 0x1,
    ];

    let scramble = rng.uniform_u32();
    let total_samples = samples_per_pixel * pixel_samples;
    assert_eq!(total_samples, samples.len());

    sample_gray_code(&MATRIX, scramble, samples);

    for i in 0..pixel_samples {
        let idx = i * samples_per_pixel;
        rng.shuffle_dims(&mut samples[idx..idx + samples_per_pixel], 1);
    }

    rng.shuffle_dims(samples, samples_per_pixel);
}

/// Create 2D samples using sobol2d sampling
fn sobol_2d<T: Number>(
    samples_per_pixel: usize,
    pixel_samples: usize,
    samples: &mut [Point2<T>],
    rng: &mut Rng,
) {
    const MATRIX_1: [u32; 32] = [
        0x80000000, 0x40000000, 0x20000000, 0x10000000, 0x8000000, 0x4000000, 0x2000000, 0x1000000,
        0x800000, 0x400000, 0x200000, 0x100000, 0x80000, 0x40000, 0x20000, 0x10000, 0x8000, 0x4000,
        0x2000, 0x1000, 0x800, 0x400, 0x200, 0x100, 0x80, 0x40, 0x20, 0x10, 0x8, 0x4, 0x2, 0x1,
    ];
    const MATRIX_2: [u32; 32] = [
        0x80000000, 0xc0000000, 0xa0000000, 0xf0000000, 0x88000000, 0xcc000000, 0xaa000000,
        0xff000000, 0x80800000, 0xc0c00000, 0xa0a00000, 0xf0f00000, 0x88880000, 0xcccc0000,
        0xaaaa0000, 0xffff0000, 0x80008000, 0xc000c000, 0xa000a000, 0xf000f000, 0x88008800,
        0xcc00cc00, 0xaa00aa00, 0xff00ff00, 0x80808080, 0xc0c0c0c0, 0xa0a0a0a0, 0xf0f0f0f0,
        0x88888888, 0xcccccccc, 0xaaaaaaaa, 0xffffffff,
    ];

    let scramble = Point2::new(rng.uniform_u32() as _, rng.uniform_u32() as _);
    let total_samples = samples_per_pixel * pixel_samples;
    assert_eq!(total_samples, samples.len());

    sample_gray_code2(&MATRIX_1, &MATRIX_2, scramble, samples);

    for i in 0..pixel_samples {
        let idx = i * samples_per_pixel;
        rng.shuffle_dims(&mut samples[idx..idx + samples_per_pixel], 1);
    }

    rng.shuffle_dims(samples, samples_per_pixel);
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

/// Compute the nth gray code number
fn gray_code(n: u32) -> u32 {
    (n >> 1) ^ n
}

/// Generate n sample values from a generator matrix in gray code ordering
fn sample_gray_code<T: Number>(lines: &[u32], scramble: u32, output: &mut [T]) {
    let mut v = scramble;

    for i in 0..output.len() {
        output[i] = T::ONE_MINUS_EPSILON.min(T::cast((v as f32) * 2.0f32.powi(-32)));
        v ^= lines[(i + 1).trailing_zeros() as usize];
    }
}

/// Generate n sample values from a pair of generator matrixes in gray code ordering.
/// 2D version of sample_gray_code.
fn sample_gray_code2<T: Number>(
    mat1: &[u32],
    mat2: &[u32],
    scramble: Point2i,
    output: &mut [Point2<T>],
) {
    let mut v = [scramble.x as u32, scramble.y as u32];

    for i in 0..output.len() {
        output[i].x = T::ONE_MINUS_EPSILON.min(T::cast((v[0] as f32) * 2.0f32.powi(-32)));
        output[i].y = T::ONE_MINUS_EPSILON.min(T::cast((v[1] as f32) * 2.0f32.powi(-32)));

        v[0] ^= mat1[(i + 1).trailing_zeros() as usize];
        v[1] ^= mat2[(i + 1).trailing_zeros() as usize];
    }
}
