use glam::DVec3;

use crate::{ray::Ray, scene::Scene, utils::rand_unit_disk};

pub struct Camera {
    origin: DVec3,
    lower_left: DVec3,
    horizontal: DVec3,
    vertical: DVec3,
    pub width: u32,
    pub height: u32,
    u: DVec3,
    v: DVec3,
    lens_radius: Option<f64>,
}

impl Camera {
    pub fn new(scene: &Scene) -> Camera {
        let theta = scene.vfov.to_radians();
        let h = (theta / 2.0).tan();

        let viewport_height = 2.0 * h;
        let viewport_width = viewport_height * scene.aspect_ratio;

        let w = (scene.origin - scene.look_at).normalize();
        let u = scene.up_direction.cross(w).normalize();
        let v = w.cross(u);

        let focal_distance = scene.focal_distance.unwrap_or(1.0);
        let horizontal = focal_distance * viewport_width * u;
        let vertical = focal_distance * viewport_height * v;
        let lower_left = scene.origin - horizontal / 2.0 - vertical / 2.0 - focal_distance * w;

        let lens_radius = if scene.focal_distance.is_some() {
            Some(scene.aperture / 2.0)
        } else {
            None
        };

        Camera {
            origin: scene.origin,
            lower_left,
            horizontal,
            vertical,
            width: viewport_width.floor() as _,
            height: viewport_height.floor() as _,
            u,
            v,
            lens_radius,
        }
    }

    pub fn ray(&self, s: f64, t: f64) -> Ray {
        let offset = if let Some(radius) = self.lens_radius {
            let rd = radius * rand_unit_disk();
            self.u * rd.x + self.v * rd.y
        } else {
            DVec3::ZERO
        };

        Ray {
            origin: self.origin + offset,
            direction: self.lower_left + s * self.horizontal + t * self.vertical
                - self.origin
                - offset,
        }
    }
}
