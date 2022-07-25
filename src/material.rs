use glam::DVec3;

use crate::{hit::HitRecord, rand_sphere_point, ray::Ray};

pub struct ScatterRecord {
    pub attenuation: DVec3,
    pub scattered: Ray,
}
pub trait Material: Send + Sync {
    fn scatter(&self, ray: Ray, hit_record: &HitRecord) -> Option<ScatterRecord>;
}

pub struct Lambertian {
    pub albedo: DVec3,
}

impl Material for Lambertian {
    fn scatter(&self, ray: Ray, hit_record: &HitRecord) -> Option<ScatterRecord> {
        let mut direction = hit_record.normal + rand_sphere_point();

        let e = 1.0e-8;
        if direction.x.abs() < e && direction.y.abs() < e && direction.x.abs() < e {
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

pub struct Metal {
    pub albedo: DVec3,
}

impl Material for Metal {
    fn scatter(&self, ray: Ray, hit_record: &HitRecord) -> Option<ScatterRecord> {
        let reflected = reflect(ray.direction.normalize(), hit_record.normal);
        let scattered = Ray {
            origin: hit_record.point,
            direction: reflected,
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
