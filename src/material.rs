use std::fmt::Debug;

use glam::DVec3;

use crate::{
    hit::HitRecord,
    ray::Ray,
    scene_format::RandomKind,
    utils::{near_zero, rand_hemisphere_point},
};

pub struct ScatterRecord {
    pub attenuation: DVec3,
    pub scattered: Ray,
}
pub trait Material: Send + Sync + Debug {
    fn scatter(&self, ray: Ray, hit_record: &HitRecord) -> Option<ScatterRecord>;
}

#[derive(Debug, Default, Clone, Copy)]
pub struct Lambertian {
    pub albedo: DVec3,
    pub random_kind: RandomKind,
}

impl Material for Lambertian {
    fn scatter(&self, _ray: Ray, hit_record: &HitRecord) -> Option<ScatterRecord> {
        let mut direction = hit_record.normal + rand_hemisphere_point(hit_record.normal);

        if near_zero(direction) {
            direction = hit_record.normal;
        }

        let scattered = Ray {
            origin: hit_record.point,
            direction,
        };
        let attenuation = self.albedo;

        Some(ScatterRecord {
            attenuation,
            scattered,
        })
    }
}

fn reflect(vec: DVec3, normal: DVec3) -> DVec3 {
    vec - 2.0 * vec.dot(normal) * normal
}

#[derive(Debug, Default, Clone, Copy)]
pub struct Metal {
    pub albedo: DVec3,
    pub fuzz: f64,
    pub random_kind: RandomKind,
}

impl Material for Metal {
    fn scatter(&self, ray: Ray, hit_record: &HitRecord) -> Option<ScatterRecord> {
        let reflected = reflect(ray.direction.normalize(), hit_record.normal);
        let scattered = Ray {
            origin: hit_record.point,
            direction: reflected + self.fuzz * rand_hemisphere_point(hit_record.normal),
        };
        let attenuation = self.albedo;

        if scattered.direction.dot(hit_record.normal) > 0.0 {
            Some(ScatterRecord {
                attenuation,
                scattered,
            })
        } else {
            None
        }
    }
}
