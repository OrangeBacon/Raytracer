use glam::DVec3;

use crate::ray::Ray;

pub struct Camera {
    origin: DVec3,
    lower_left: DVec3,
    horizontal: DVec3,
    vertical: DVec3,
    pub width: u32,
    pub height: u32,
}

impl Camera {
    pub fn new(
        look_from: DVec3,
        look_at: DVec3,
        v_up: DVec3,
        vfov: f64,
        aspect_ratio: f64,
    ) -> Camera {
        let theta = vfov.to_radians();
        let h = (theta / 2.0).tan();

        let viewport_height = 2.0 * h;
        let viewport_width = viewport_height * aspect_ratio;

        let w = (look_from - look_at).normalize();
        let u = v_up.cross(w).normalize();
        let v = w.cross(u);

        let horizontal = viewport_width * u;
        let vertical = viewport_height * v;
        let lower_left = look_from - horizontal / 2.0 - vertical / 2.0 - w;

        Camera {
            origin: look_from,
            lower_left,
            horizontal,
            vertical,
            width: viewport_width.floor() as _,
            height: viewport_height.floor() as _,
        }
    }

    pub fn ray(&self, s: f64, t: f64) -> Ray {
        Ray {
            origin: self.origin,
            direction: self.lower_left + s * self.horizontal + t * self.vertical - self.origin,
        }
    }
}
