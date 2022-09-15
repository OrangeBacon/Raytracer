use geometry::{Number, Point2, Point2i, Rng};

use crate::pixel_sampler::{PixelSampler, PixelSamplerData, PixelSamplerImpl};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default, Hash)]
pub struct StratifiedSampler {
    x_pixel_samples: usize,
    y_pixel_samples: usize,
    jitter_samples: bool,
}

impl StratifiedSampler {
    /// Create a new stratified sampler
    pub fn new<T: Number>(
        x_pixel_samples: usize,
        y_pixel_samples: usize,
        jitter_samples: bool,
        sampled_dimensions: usize,
    ) -> PixelSampler<T, StratifiedSampler> {
        let sampler = StratifiedSampler {
            x_pixel_samples,
            y_pixel_samples,
            jitter_samples,
        };
        PixelSampler::new(
            x_pixel_samples * y_pixel_samples,
            sampled_dimensions,
            sampler,
        )
    }
}

impl<T: Number> PixelSamplerImpl<T> for StratifiedSampler {
    fn start_pixel(&mut self, data: &mut PixelSamplerData<T>, _point: Point2i) {
        for sample in &mut data.samples_one {
            stratified_sample_one(sample, &mut data.rng, self.jitter_samples);
            shuffle(sample, 1, &mut data.rng)
        }

        for sample in &mut data.samples_two {
            stratified_sample_two(
                sample,
                self.x_pixel_samples,
                self.y_pixel_samples,
                &mut data.rng,
                self.jitter_samples,
            );
            shuffle(sample, 1, &mut data.rng)
        }

        for i in 0..data.samples_one_array_sizes.len() {
            for j in 0..data.samples_per_pixel {
                let count = data.samples_one_array_sizes[i];
                stratified_sample_one(
                    &mut data.data.sample_array_one[i][j * count..count],
                    &mut data.rng,
                    self.jitter_samples,
                )
            }
        }

        for i in 0..data.samples_two_array_sizes.len() {
            for j in 0..data.samples_per_pixel {
                let count = data.samples_two_array_sizes[i];
                latin_hypercube(
                    &mut data.data.sample_array_two[i][j * count..count],
                    count,
                    2,
                    &mut data.rng,
                )
            }
        }
    }
}

/// Create 1D stratified samples
fn stratified_sample_one<T: Number>(samples: &mut [T], rng: &mut Rng, jitter: bool) {
    let inv = T::ONE / T::cast(samples.len() as i32);

    for (idx, sample) in samples.iter_mut().enumerate() {
        let delta = if jitter { rng.uniform_float() } else { T::HALF };
        *sample = T::ONE_MINUS_EPSILON.min((T::cast(idx as i32) + delta) * inv);
    }
}

/// Create 2D stratified samples
fn stratified_sample_two<T: Number>(
    samples: &mut [Point2<T>],
    x: usize,
    y: usize,
    rng: &mut Rng,
    jitter: bool,
) {
    assert_eq!(
        x * y,
        samples.len(),
        "Stratified sample dimensions are inconsistent"
    );

    let dx = T::ONE / T::cast(x as i32);
    let dy = T::ONE / T::cast(y as i32);

    for (idx, sample) in samples.iter_mut().enumerate() {
        let jx = if jitter { rng.uniform_float() } else { T::HALF };
        let jy = if jitter { rng.uniform_float() } else { T::HALF };

        *sample = Point2::new(
            T::ONE_MINUS_EPSILON.min((T::cast(idx as i32) + jx) * dx),
            T::ONE_MINUS_EPSILON.min((T::cast(idx as i32) + jy) * dy),
        );
    }
}

/// Shuffle blocks of data of length n_dimensions within the given data
fn shuffle<T: Copy>(data: &mut [T], n_dimensions: usize, rng: &mut Rng) {
    for i in 0..data.len() {
        let other = i + rng.uniform_u32_limit((data.len() - i) as _) as usize;
        for j in 0..n_dimensions {
            data.swap(n_dimensions * i + j, n_dimensions * other + j);
        }
    }
}

/// Create a number of samples of arbitrary dimension using latin hypercube sampling
fn latin_hypercube<T: Number>(
    samples: &mut [Point2<T>],
    n_samples: usize,
    n_dimensions: usize,
    rng: &mut Rng,
) {
    let inv = T::ONE / T::cast(n_samples as i32);
    for i in 0..n_samples {
        for j in 0..n_dimensions {
            let sj = (T::cast(i as i32) + rng.uniform_float()) * inv;
            let idx = n_dimensions * i + j;
            let div = idx / 2;
            let rem = idx % 2;
            *(if rem == 0 {
                &mut samples[div].x
            } else {
                &mut samples[div].y
            }) = T::ONE_MINUS_EPSILON.min(sj);
        }
    }

    for i in 0..n_dimensions {
        for j in 0..n_samples {
            let other = j + rng.uniform_u32_limit((n_samples - j) as _) as usize;

            // All this below is in reality just a mem::swap if you treat 
            // [Point2<T>; N] == [[T; 2]; N] == [T; N * 2], but below does it without
            // having to assume anything about the layout of Point2<T>
            let idx = n_dimensions * j + i;
            let div_a = idx / 2;
            let rem_a = idx % 2;
            let tmp = if rem_a == 0 {
                samples[div_a].x
            } else {
                samples[div_a].y
            };

            let idx = n_dimensions * other + i;
            let div_b = idx / 2;
            let rem_b = idx % 2;
            let b = if rem_b == 0 {
                samples[div_b].x
            } else {
                samples[div_b].y
            };

            *(if rem_a == 0 {
                &mut samples[div_a].x
            } else {
                &mut samples[div_a].y
            }) = b;

            *(if rem_b == 0 {
                &mut samples[div_b].x
            } else {
                &mut samples[div_b].y
            }) = tmp;
        }
    }
}
