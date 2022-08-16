use glam::DVec3;

use crate::ray::Ray;

#[derive(Debug, Clone, Copy)]
pub struct Aabb {
    pub minimum: DVec3,
    pub maximum: DVec3,
}

impl Aabb {
    pub fn hit(&self, ray: Ray, mut t_min: f64, mut t_max: f64) -> bool {
        for a in 0..3 {
            let inv_d = 1.0 / ray.direction[a];

            let mut t0 = (self.minimum[a] - ray.origin[a]) * inv_d;
            let mut t1 = (self.maximum[a] - ray.origin[a]) * inv_d;
            if inv_d < 0.0 {
                std::mem::swap(&mut t0, &mut t1);
            }

            t_min = if t0 > t_min { t0 } else { t_min };
            t_max = if t1 < t_max { t1 } else { t_max };
            if t_max <= t_min {
                return false;
            }
        }

        true
    }

    pub fn surrounds(&self, other: Aabb) -> Aabb {
        Aabb {
            minimum: self.minimum.min(other.minimum),
            maximum: self.maximum.max(other.maximum),
        }
    }
}
