use std::{fmt::Debug, path::Path, sync::Arc};

use anyhow::Result;
use glam::{dvec3, DVec2, DVec3};
use image::RgbImage;

pub trait Texture: Send + Sync + Debug {
    fn value(&self, uv: DVec2, point: DVec3) -> DVec3;
}

#[derive(Debug)]
pub struct SolidColour {
    pub colour: DVec3,
}

impl Texture for SolidColour {
    fn value(&self, _: DVec2, _: DVec3) -> DVec3 {
        self.colour
    }
}

#[derive(Debug)]
pub struct Checker {
    pub odd: Arc<dyn Texture>,
    pub even: Arc<dyn Texture>,
}

impl Texture for Checker {
    fn value(&self, uv: DVec2, point: DVec3) -> DVec3 {
        let sines = point * 10.0;
        let sines = sines.x.sin() * sines.y.sin() * sines.z.sin();

        if sines < 0.0 {
            self.odd.value(uv, point)
        } else {
            self.even.value(uv, point)
        }
    }
}

#[derive(Debug)]
pub struct ImageTexture {
    image: RgbImage,
}

impl ImageTexture {
    pub fn new(file_name: impl AsRef<Path>) -> Result<Self> {
        let image = image::open(file_name)?;

        Ok(ImageTexture {
            image: image.into_rgb8(),
        })
    }
}

impl Texture for ImageTexture {
    fn value(&self, uv: DVec2, _: DVec3) -> DVec3 {
        let u = uv.x.clamp(0.0, 1.0);
        let v = 1.0 - uv.y.clamp(0.0, 1.0);
        let mut x = (u * self.image.width() as f64) as u32;
        let mut y = (v * self.image.height() as f64) as u32;

        if x >= self.image.width() {
            x = self.image.width() - 1
        }
        if y >= self.image.height() {
            y = self.image.height() - 1
        }

        let data = self.image.get_pixel(x, y);
        dvec3(data[0] as f64, data[1] as f64, data[2] as f64) / 255.0
    }
}
