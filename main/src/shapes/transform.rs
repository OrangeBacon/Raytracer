use std::sync::Arc;

use glam::{dvec3, DVec3};

use crate::{aabb::Aabb, ray::Ray};

use super::{HitRecord, Hittable};

#[derive(Debug)]
pub struct Translate {
    pub child: Arc<dyn Hittable>,
    pub offset: DVec3,
}

impl Hittable for Translate {
    fn hit(&self, ray: Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let ray = Ray {
            origin: ray.origin - self.offset,
            direction: ray.direction,
            time: ray.time,
        };

        let mut hit = self.child.hit(ray, t_min, t_max)?;
        hit.point += self.offset;
        hit.set_face_normal(ray, hit.normal);

        Some(hit)
    }

    fn bounding_box(&self, time0: f64, time1: f64) -> Option<Aabb> {
        let inner = self.child.bounding_box(time0, time1)?;

        Some(Aabb {
            minimum: inner.minimum + self.offset,
            maximum: inner.maximum + self.offset,
        })
    }
}

#[derive(Debug)]
pub struct RotateY {
    child: Arc<dyn Hittable>,
    bounding_box: Option<Aabb>,
    sin_theta: f64,
    cos_theta: f64,
}

impl RotateY {
    pub fn new(child: Arc<dyn Hittable>, angle: f64) -> Self {
        let radians = angle.to_radians();
        let sin_theta = radians.sin();
        let cos_theta = radians.cos();
        let bounding_box = child.bounding_box(0.0, 1.0).map(|bbox| {
            let mut minimum = DVec3::splat(f64::INFINITY);
            let mut maximum = DVec3::splat(f64::NEG_INFINITY);

            for i in 0..2 {
                for j in 0..2 {
                    for k in 0..2 {
                        let ijk = dvec3(i as f64, j as f64, k as f64);
                        let xyz = ijk * bbox.maximum + (1.0 - ijk) * bbox.minimum;
                        let newx = cos_theta * xyz.x + sin_theta * xyz.z;
                        let newz = -sin_theta * xyz.x + cos_theta * xyz.z;
                        let tester = dvec3(newx, xyz.y, newz);

                        for c in 0..3 {
                            minimum[c] = minimum[c].min(tester[c]);
                            maximum[c] = maximum[c].max(tester[c]);
                        }
                    }
                }
            }

            Aabb { minimum, maximum }
        });

        Self {
            child,
            sin_theta,
            cos_theta,
            bounding_box,
        }
    }
}

impl Hittable for RotateY {
    fn hit(&self, ray: Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let origin = dvec3(
            self.cos_theta * ray.origin.x - self.sin_theta * ray.origin.z,
            ray.origin.y,
            self.sin_theta * ray.origin.x + self.cos_theta * ray.origin.z,
        );

        let direction = dvec3(
            self.cos_theta * ray.direction.x - self.sin_theta * ray.direction.z,
            ray.direction.y,
            self.sin_theta * ray.direction.x + self.cos_theta * ray.direction.z,
        );

        let rotated = Ray {
            origin,
            direction,
            time: ray.time,
        };

        if let Some(hit) = self.child.hit(rotated, t_min, t_max) {
            let point = dvec3(
                self.cos_theta * hit.point.x + self.sin_theta * hit.point.z,
                hit.point.y,
                -self.sin_theta * hit.point.x + self.cos_theta * hit.point.z,
            );

            let normal = dvec3(
                self.cos_theta * hit.normal.x + self.sin_theta * hit.normal.z,
                hit.normal.y,
                -self.sin_theta * hit.normal.x + self.cos_theta * hit.normal.z,
            );

            let hit = HitRecord::new(rotated, point, normal, hit.t, hit.material, hit.uv);

            return Some(hit);
        }

        None
    }

    fn bounding_box(&self, _time0: f64, _time1: f64) -> Option<Aabb> {
        self.bounding_box
    }
}
