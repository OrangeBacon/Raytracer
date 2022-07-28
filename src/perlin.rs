use glam::{dvec3, DVec2, DVec3};
use rand::Rng;

use crate::texture::Texture;

#[derive(Debug)]
struct Perlin {
    values: Vec<DVec3>,
    x: Vec<usize>,
    y: Vec<usize>,
    z: Vec<usize>,
}

impl Perlin {
    const POINT_COUNT: usize = 256;

    pub fn new() -> Self {
        let mut values = vec![DVec3::ZERO; Self::POINT_COUNT];

        let mut rng = rand::thread_rng();
        for value in &mut values {
            *value = DVec3::from_array([
                rng.gen_range(-1.0..=1.0),
                rng.gen_range(-1.0..=1.0),
                rng.gen_range(-1.0..=1.0),
            ])
            .normalize();
        }

        Self {
            values,
            x: Self::generate_permutation(),
            y: Self::generate_permutation(),
            z: Self::generate_permutation(),
        }
    }

    fn generate_permutation() -> Vec<usize> {
        let mut data = vec![0; Self::POINT_COUNT];

        for (idx, data) in data.iter_mut().enumerate() {
            *data = idx;
        }

        let mut rng = rand::thread_rng();
        for idx in (1..data.len()).rev() {
            let target = rng.gen_range(0..=idx);
            let tmp = data[idx];
            data[idx] = data[target];
            data[target] = tmp;
        }

        data
    }

    pub fn noise(&self, point: DVec3) -> f64 {
        let uvw = point - point.floor();
        let ijk = point.floor().as_ivec3();

        let mut data = [[[DVec3::ZERO; 2]; 2]; 2];

        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    data[i][j][k] = self.values[self.x[((i as i32 + ijk.x) & 255) as usize]
                        ^ self.y[((j as i32 + ijk.y) & 255) as usize]
                        ^ self.z[((k as i32 + ijk.z) & 255) as usize]]
                }
            }
        }

        Self::perlin_interp(data, uvw)
    }

    fn perlin_interp(data: [[[DVec3; 2]; 2]; 2], uvw: DVec3) -> f64 {
        let m_uvw = uvw * uvw * (3.0 - 2.0 * uvw);

        let mut accum = 0.0;
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    let weight = uvw - dvec3(i as f64, j as f64, k as f64);

                    accum += (i as f64 * m_uvw.x + (1.0 - i as f64) * (1.0 - m_uvw.x))
                        * (j as f64 * m_uvw.y + (1.0 - j as f64) * (1.0 - m_uvw.y))
                        * (k as f64 * m_uvw.z + (1.0 - k as f64) * (1.0 - m_uvw.z))
                        * data[i][j][k].dot(weight);
                }
            }
        }

        accum
    }

    fn turbulence(&self, point: DVec3, depth: usize) -> f64 {
        let mut accum = 0.0;
        let mut temp = point;
        let mut weight = 1.0;

        for _ in 0..depth {
            accum += weight * self.noise(temp);
            weight *= 0.5;
            temp *= 2.0;
        }

        accum.abs()
    }
}

#[derive(Debug)]
pub struct NoiseTexture {
    noise: Perlin,
    scale: f64,
}

impl NoiseTexture {
    pub fn new(scale: f64) -> Self {
        Self {
            noise: Perlin::new(),
            scale,
        }
    }
}

impl Texture for NoiseTexture {
    fn value(&self, _: DVec2, point: DVec3) -> DVec3 {
        let noise = self.noise.turbulence(self.scale * point, 7);
        DVec3::ONE * 0.5 * (1.0 + (self.scale * point.z + 10.0 * noise).sin())
    }
}
