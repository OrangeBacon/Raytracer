//! Axis aligned rectangles

use std::sync::Arc;

use glam::{dvec2, dvec3, DVec3};

use crate::{aabb::Aabb, material::Material, ray::Ray};

use super::{HitRecord, Hittable};

#[derive(Debug)]
pub struct XyRect {
    pub mat: Arc<dyn Material>,
    pub x0: f64,
    pub x1: f64,
    pub y0: f64,
    pub y1: f64,
    pub k: f64,
}

impl Hittable for XyRect {
    fn hit(&self, ray: Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let t = (self.k - ray.origin.z) / ray.direction.z;
        if t < t_min || t > t_max {
            return None;
        }

        let x = ray.origin.x + t * ray.direction.x;
        let y = ray.origin.y + t * ray.direction.y;
        if x < self.x0 || x > self.x1 || y < self.y0 || y > self.y1 {
            return None;
        }

        Some(HitRecord::new(
            ray,
            ray.at(t),
            DVec3::Z,
            t,
            Arc::clone(&self.mat),
            dvec2(
                (x - self.x0) / (self.x1 - self.x0),
                (y - self.y0) / (self.y1 - self.y0),
            ),
        ))
    }

    fn bounding_box(&self, _: f64, _: f64) -> Option<Aabb> {
        Some(Aabb {
            minimum: dvec3(self.x0, self.y0, self.k - 0.0001),
            maximum: dvec3(self.x1, self.y1, self.k + 0.0001),
        })
    }
}

#[derive(Debug)]
pub struct XzRect {
    pub mat: Arc<dyn Material>,
    pub x0: f64,
    pub x1: f64,
    pub z0: f64,
    pub z1: f64,
    pub k: f64,
}

impl Hittable for XzRect {
    fn hit(&self, ray: Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let t = (self.k - ray.origin.y) / ray.direction.y;
        if t < t_min || t > t_max {
            return None;
        }

        let x = ray.origin.x + t * ray.direction.x;
        let z = ray.origin.z + t * ray.direction.z;
        if x < self.x0 || x > self.x1 || z < self.z0 || z > self.z1 {
            return None;
        }

        Some(HitRecord::new(
            ray,
            ray.at(t),
            DVec3::Y,
            t,
            Arc::clone(&self.mat),
            dvec2(
                (x - self.x0) / (self.x1 - self.x0),
                (z - self.z0) / (self.z1 - self.z0),
            ),
        ))
    }

    fn bounding_box(&self, _: f64, _: f64) -> Option<Aabb> {
        Some(Aabb {
            minimum: dvec3(self.x0, self.k - 0.0001, self.z0),
            maximum: dvec3(self.x1, self.k + 0.0001, self.z1),
        })
    }
}

#[derive(Debug)]
pub struct YzRect {
    pub mat: Arc<dyn Material>,
    pub y0: f64,
    pub y1: f64,
    pub z0: f64,
    pub z1: f64,
    pub k: f64,
}

impl Hittable for YzRect {
    fn hit(&self, ray: Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let t = (self.k - ray.origin.x) / ray.direction.x;
        if t < t_min || t > t_max {
            return None;
        }

        let y = ray.origin.y + t * ray.direction.y;
        let z = ray.origin.z + t * ray.direction.z;
        if y < self.y0 || y > self.y1 || z < self.z0 || z > self.z1 {
            return None;
        }

        Some(HitRecord::new(
            ray,
            ray.at(t),
            DVec3::X,
            t,
            Arc::clone(&self.mat),
            dvec2(
                (y - self.y0) / (self.y1 - self.y0),
                (z - self.z0) / (self.z1 - self.z0),
            ),
        ))
    }

    fn bounding_box(&self, _: f64, _: f64) -> Option<Aabb> {
        Some(Aabb {
            minimum: dvec3(self.k - 0.0001, self.y0, self.z0),
            maximum: dvec3(self.k + 0.0001, self.y1, self.z1),
        })
    }
}
