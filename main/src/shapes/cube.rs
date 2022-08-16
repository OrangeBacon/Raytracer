use std::sync::Arc;

use glam::DVec3;

use crate::{aabb::Aabb, material::Material, ray::Ray};

use super::{
    aa_rect::{XyRect, XzRect, YzRect},
    HitRecord, Hittable, HittableVec,
};

#[derive(Debug)]
pub struct Cube {
    min: DVec3,
    max: DVec3,
    sides: HittableVec,
}

impl Cube {
    pub fn new(min: DVec3, max: DVec3, mat: Arc<dyn Material>) -> Self {
        let sides = HittableVec {
            objects: vec![
                XyRect::new(&mat, min.x, max.x, min.y, max.y, min.z),
                XyRect::new(&mat, min.x, max.x, min.y, max.y, max.z),
                XzRect::new(&mat, min.x, max.x, min.z, max.z, min.y),
                XzRect::new(&mat, min.x, max.x, min.z, max.z, max.y),
                YzRect::new(&mat, min.y, max.y, min.z, max.z, min.x),
                YzRect::new(&mat, min.y, max.y, min.z, max.z, max.x),
            ],
        };

        Self { min, max, sides }
    }
}

impl Hittable for Cube {
    fn hit(&self, ray: Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        self.sides.hit(ray, t_min, t_max)
    }

    fn bounding_box(&self, _: f64, _: f64) -> Option<Aabb> {
        Some(Aabb {
            minimum: self.min,
            maximum: self.max,
        })
    }
}
