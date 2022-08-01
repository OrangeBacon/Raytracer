use std::sync::Arc;

use glam::{DVec2, DVec3};

use crate::{aabb::Aabb, material::Material, ray::Ray};

use super::{HitRecord, Hittable};

/// Fog/mist etc volumetric object
#[derive(Debug)]
pub struct Volume {
    boundary: Arc<dyn Hittable>,
    phase: Arc<dyn Material>,
    neg_inv_density: f64,
}

impl Volume {
    pub fn new(boundary: Arc<dyn Hittable>, density: f64, texture: Arc<dyn Material>) -> Self {
        Self {
            boundary,
            phase: texture,
            neg_inv_density: -1.0 / density,
        }
    }
}

impl Hittable for Volume {
    fn hit(&self, ray: Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut rec1 = self.boundary.hit(ray, f64::NEG_INFINITY, f64::INFINITY)?;
        let mut rec2 = self.boundary.hit(ray, rec1.t + 0.0001, f64::INFINITY)?;

        if rec1.t < t_min {
            rec1.t = t_min
        }
        if rec2.t > t_max {
            rec2.t = t_max
        }

        if rec1.t >= rec2.t {
            return None;
        }

        if rec1.t < 0.0 {
            rec1.t = 0.0
        }

        let ray_len = ray.direction.length();
        let distance_inside = (rec2.t - rec1.t) * ray_len;
        let hit_distance = self.neg_inv_density * rand::random::<f64>().ln();

        if hit_distance > distance_inside {
            return None;
        }

        let t = rec1.t + hit_distance / ray_len;
        let point = ray.at(t);

        Some(HitRecord::new(
            ray,
            point,
            DVec3::X,
            t,
            Arc::clone(&self.phase),
            DVec2::ZERO,
        ))
    }

    fn bounding_box(&self, time0: f64, time1: f64) -> Option<Aabb> {
        self.boundary.bounding_box(time0, time1)
    }
}
