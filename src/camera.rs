use glam::{dvec3, DVec3};

use crate::ray::Ray;

pub struct Camera {
    origin: DVec3,
    lower_left: DVec3,
    horizontal: DVec3,
    vertical: DVec3,
}

impl Camera {
    pub fn new(width: u32, height: u32) -> Camera {
        let aspect_ratio = width as f64 / height as f64;
        let viewport_height = 2.0;
        let viewport_width = viewport_height * aspect_ratio;
        let focal_length = 1.0;

        let origin = DVec3::ZERO;
        let horizontal = dvec3(viewport_width, 0.0, 0.0);
        let vertical = dvec3(0.0, viewport_height, 0.0);
        let lower_left = origin - horizontal / 2.0 - vertical / 2.0 - dvec3(0.0, 0.0, focal_length);

        Camera {
            origin,
            lower_left,
            horizontal,
            vertical,
        }
    }

    pub fn ray(&self, u: f64, v: f64) -> Ray {
        Ray {
            origin: self.origin,
            direction: self.lower_left + u * self.horizontal + v * self.vertical - self.origin,
        }
    }
}
