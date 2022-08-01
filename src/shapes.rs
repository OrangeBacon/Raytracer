pub mod aa_rect;
pub mod sphere;
pub mod cube;
pub mod transform;

use std::{fmt::Debug, sync::Arc};

use glam::{DVec2, DVec3};

use crate::{aabb::Aabb, material::Material, ray::Ray};

pub struct HitRecord {
    pub point: DVec3,
    pub normal: DVec3,
    pub t: f64,
    pub front_face: bool,
    pub material: Arc<dyn Material>,
    pub uv: DVec2,
}

pub trait Hittable: Send + Sync + Debug {
    fn hit(&self, ray: Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;

    fn bounding_box(&self, time0: f64, time1: f64) -> Option<Aabb>;
}

impl HitRecord {
    pub fn new(
        ray: Ray,
        point: DVec3,
        outward_normal: DVec3,
        t: f64,
        material: Arc<dyn Material>,
        uv: DVec2,
    ) -> Self {
        let mut record = HitRecord {
            point,
            material,
            normal: DVec3::ZERO,
            t,
            front_face: false,
            uv,
        };
        record.set_face_normal(ray, outward_normal);
        record
    }

    fn set_face_normal(&mut self, ray: Ray, outward_normal: DVec3) {
        self.front_face = ray.direction.dot(outward_normal) < 0.0;
        self.normal = if self.front_face {
            outward_normal
        } else {
            -outward_normal
        }
    }
}

#[derive(Debug)]
pub struct HittableVec {
    pub objects: Vec<Arc<dyn Hittable>>,
}

impl Hittable for HittableVec {
    fn hit(&self, ray: Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut closest = t_max;
        let mut output = None;

        for object in &self.objects {
            if let Some(record) = object.hit(ray, t_min, closest) {
                closest = record.t;
                output = Some(record);
            }
        }

        output
    }

    fn bounding_box(&self, time0: f64, time1: f64) -> Option<Aabb> {
        let mut result = None;

        for object in &self.objects {
            if let Some(aabb) = object.bounding_box(time0, time1) {
                match result {
                    None => result = Some(aabb),
                    Some(a) => result = Some(a.surrounds(aabb)),
                }
            } else {
                return None;
            }
        }

        result
    }
}


