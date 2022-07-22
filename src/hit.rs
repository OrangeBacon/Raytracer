use glam::Vec3;

use crate::ray::Ray;

pub struct HitRecord {
    pub point: Vec3,
    pub normal: Vec3,
    pub t: f32,
    pub front_face: bool,
}

pub trait Hittable: Send + Sync {
    fn hit(&self, ray: Ray, t_min: f32, t_max: f32) -> Option<HitRecord>;
}

impl HitRecord {
    pub fn new(ray: Ray, point: Vec3, outward_normal: Vec3, t: f32) -> Self {
        let mut record = HitRecord {
            point,
            normal: Vec3::ZERO,
            t,
            front_face: false,
        };
        record.set_face_normal(ray, outward_normal);
        record
    }

    fn set_face_normal(&mut self, ray: Ray, outward_normal: Vec3) {
        self.front_face = ray.direction.dot(outward_normal) < 0.0;
        self.normal = if self.front_face {
            outward_normal
        } else {
            -outward_normal
        }
    }
}

#[derive(Clone, Copy)]
pub struct Sphere {
    pub centre: Vec3,
    pub radius: f32,
}

impl Hittable for Sphere {
    fn hit(&self, ray: Ray, t_min: f32, t_max: f32) -> Option<HitRecord> {
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
        ));
    }
}

impl Hittable for [Box<dyn Hittable>] {
    fn hit(&self, ray: Ray, t_min: f32, t_max: f32) -> Option<HitRecord> {
        let mut closest = t_max;
        let mut output = None;

        for object in self {
            if let Some(record) = object.hit(ray, t_min, closest) {
                closest = record.t;
                output = Some(record);
            }
        }

        output
    }
}