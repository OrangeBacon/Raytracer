use std::{fmt::Debug, sync::Arc};

use glam::DVec3;
use rand::Rng;

use crate::{
    hit::HitRecord,
    ray::Ray,
    scene_format::RandomKind,
    texture::Texture,
    utils::{near_zero, rand_hemisphere_point, rand_sphere_point},
};

pub struct ScatterRecord {
    pub attenuation: DVec3,
    pub scattered: Ray,
}
pub trait Material: Send + Sync + Debug {
    fn scatter(&self, ray: Ray, hit_record: &HitRecord) -> Option<ScatterRecord>;
}

fn rand_sphere(normal: DVec3, kind: RandomKind) -> DVec3 {
    match kind {
        RandomKind::Sphere => rand_sphere_point(),
        RandomKind::Hemisphere => rand_hemisphere_point(normal),
    }
}

#[derive(Debug, Clone)]
pub struct Lambertian {
    pub albedo: Arc<dyn Texture>,
    pub random_kind: RandomKind,
}

impl Material for Lambertian {
    fn scatter(&self, ray: Ray, hit_record: &HitRecord) -> Option<ScatterRecord> {
        let mut direction = hit_record.normal + rand_sphere(hit_record.normal, self.random_kind);

        if near_zero(direction) {
            direction = hit_record.normal;
        }

        let scattered = Ray {
            origin: hit_record.point,
            direction,
            time: ray.time,
        };
        let attenuation = self.albedo.value(hit_record.uv, hit_record.point);

        Some(ScatterRecord {
            attenuation,
            scattered,
        })
    }
}

fn reflect(vec: DVec3, normal: DVec3) -> DVec3 {
    vec - 2.0 * vec.dot(normal) * normal
}

#[derive(Debug, Clone)]
pub struct Metal {
    pub albedo: Arc<dyn Texture>,
    pub fuzz: f64,
    pub random_kind: RandomKind,
}

impl Material for Metal {
    fn scatter(&self, ray: Ray, hit_record: &HitRecord) -> Option<ScatterRecord> {
        let reflected = reflect(ray.direction.normalize(), hit_record.normal);
        let scattered = Ray {
            origin: hit_record.point,
            direction: reflected + self.fuzz * rand_sphere(hit_record.normal, self.random_kind),
            time: ray.time,
        };
        let attenuation = self.albedo.value(hit_record.uv, hit_record.point);

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

fn refract(uv: DVec3, normal: DVec3, ratio: f64) -> DVec3 {
    let cos = (-uv).dot(normal).min(1.0);
    let perpendicular = ratio * (uv + cos * normal);
    let parallel = -(1.0 - perpendicular.length_squared()).abs().sqrt() * normal;

    perpendicular + parallel
}

/// Schlick's reflectance approximation
fn reflectance(cos: f64, refractive_index: f64) -> f64 {
    let r0 = (1.0 - refractive_index) / (1.0 + refractive_index);
    let r0 = r0 * r0;
    r0 + (1.0 - r0) * (1.0 - cos).powi(5)
}

#[derive(Debug, Default)]
pub struct Dielectric {
    pub refractive_index: f64,
}

impl Material for Dielectric {
    fn scatter(&self, ray: Ray, hit_record: &HitRecord) -> Option<ScatterRecord> {
        let attenuation = DVec3::ONE;
        let refraction_ratio = if hit_record.front_face {
            1.0 / self.refractive_index
        } else {
            self.refractive_index
        };

        let direction = ray.direction.normalize();
        let cos = (-direction).dot(hit_record.normal).min(1.0);
        let sin = (1.0 - cos * cos).sqrt();

        let mut rng = rand::thread_rng();
        let direction = if refraction_ratio * sin > 1.0
            || reflectance(cos, self.refractive_index) > rng.gen()
        {
            reflect(direction, hit_record.normal)
        } else {
            refract(direction, hit_record.normal, refraction_ratio)
        };

        let scattered = Ray {
            origin: hit_record.point,
            direction,
            time: ray.time,
        };

        Some(ScatterRecord {
            attenuation,
            scattered,
        })
    }
}
