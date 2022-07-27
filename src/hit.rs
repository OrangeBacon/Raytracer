use std::{fmt::Debug, sync::Arc};

use glam::DVec3;

use crate::{aabb::AABB, material::Material, ray::Ray};

pub struct HitRecord {
    pub point: DVec3,
    pub normal: DVec3,
    pub t: f64,
    pub front_face: bool,
    pub material: Arc<dyn Material>,
}

pub trait Hittable: Send + Sync + Debug {
    fn hit(&self, ray: Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;

    fn bounding_box(&self, time0: f64, time1: f64) -> Option<AABB>;
}

impl HitRecord {
    pub fn new(
        ray: Ray,
        point: DVec3,
        outward_normal: DVec3,
        t: f64,
        material: Arc<dyn Material>,
    ) -> Self {
        let mut record = HitRecord {
            point,
            material,
            normal: DVec3::ZERO,
            t,
            front_face: false,
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

#[derive(Clone, Debug)]
pub struct Sphere {
    pub centre: DVec3,
    pub radius: f64,
    pub material: Arc<dyn Material>,
}

impl Hittable for Sphere {
    fn hit(&self, ray: Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let oc = ray.origin - self.centre;
        let a = ray.direction.length_squared();
        let half_b = oc.dot(ray.direction);
        let c = oc.length_squared() - self.radius * self.radius;
        let discriminant = half_b * half_b - a * c;
        if discriminant < 0.0 {
            return None;
        }

        let sqrt_discriminant = discriminant.sqrt();

        let mut root = (-half_b - sqrt_discriminant) / a;
        if root < t_min || t_max < root {
            root = (-half_b + sqrt_discriminant) / a;
            if root < t_min || t_max < root {
                return None;
            }
        }

        let point = ray.at(root);
        return Some(HitRecord::new(
            ray,
            point,
            (point - self.centre) / self.radius,
            root,
            Arc::clone(&self.material),
        ));
    }

    fn bounding_box(&self, _time0: f64, _time1: f64) -> Option<AABB> {
        Some(AABB {
            minimum: self.centre - DVec3::splat(self.radius),
            maximum: self.centre + DVec3::splat(self.radius),
        })
    }
}

#[derive(Debug)]
pub struct MovingSphere {
    pub centre0: DVec3,
    pub centre1: DVec3,
    pub time0: f64,
    pub time1: f64,
    pub radius: f64,
    pub material: Arc<dyn Material>,
}

impl MovingSphere {
    fn centre(&self, time: f64) -> DVec3 {
        self.centre0
            + ((time - self.time0) / (self.time1 - self.time0)) * (self.centre1 - self.centre0)
    }
}

impl Hittable for MovingSphere {
    fn hit(&self, ray: Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let oc = ray.origin - self.centre(ray.time);
        let a = ray.direction.length_squared();
        let half_b = oc.dot(ray.direction);
        let c = oc.length_squared() - self.radius * self.radius;
        let discriminant = half_b * half_b - a * c;
        if discriminant < 0.0 {
            return None;
        }

        let sqrt_discriminant = discriminant.sqrt();

        let mut root = (-half_b - sqrt_discriminant) / a;
        if root < t_min || t_max < root {
            root = (-half_b + sqrt_discriminant) / a;
            if root < t_min || t_max < root {
                return None;
            }
        }

        let point = ray.at(root);
        return Some(HitRecord::new(
            ray,
            point,
            (point - self.centre(ray.time)) / self.radius,
            root,
            Arc::clone(&self.material),
        ));
    }

    fn bounding_box(&self, time0: f64, time1: f64) -> Option<AABB> {
        let box0 = AABB {
            minimum: self.centre(time0) - DVec3::splat(self.radius),
            maximum: self.centre(time0) + DVec3::splat(self.radius),
        };

        let box1 = AABB {
            minimum: self.centre(time1) - DVec3::splat(self.radius),
            maximum: self.centre(time1) + DVec3::splat(self.radius),
        };

        Some(box0.surrounds(box1))
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

    fn bounding_box(&self, time0: f64, time1: f64) -> Option<AABB> {
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
