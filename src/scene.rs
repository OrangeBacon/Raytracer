//! Convert from the scene file format into the internal rendering formats

use std::{collections::HashMap, sync::Arc};

use glam::DVec3;

use crate::{
    hit::{Hittable, Sphere},
    material::{Dielectric, Lambertian, Material, Metal},
    scene_format,
};

/// Scene description
#[derive(Debug)]
pub struct Scene {
    /// Vertical field of view in degrees
    pub vfov: f64,

    /// Aspect ratio (aspect width / aspect height)
    pub aspect_ratio: f64,

    /// Height in pixels of the output image
    pub height: u32,

    /// Width in pixels of the output image
    pub width: u32,

    /// The number of rays sent for each pixel
    pub samples_per_pixel: u32,

    /// The maximum recursive depth of each ray
    pub recursive_depth: u32,

    /// output = pow(output, 1.0 / gamma), gamma correction/tone mapping
    pub gamma: f64,

    /// Where the camera is located
    pub origin: DVec3,

    /// Direction the camera is looking at
    pub look_at: DVec3,

    /// The vertical upwards direction of the camera
    pub up_direction: DVec3,

    /// All the objects in the scene
    pub world: Vec<Box<dyn Hittable>>,
}

impl Scene {
    pub fn from_file(scene: scene_format::Scene) -> Scene {
        let materials: HashMap<_, _> = scene
            .materials
            .iter()
            .map(|(k, v)| (k.to_string(), Scene::material(&scene.settings, &v)))
            .collect();

        let mut world = vec![];
        for object in scene.world {
            world.push(Scene::object(&scene.settings, object, &materials))
        }

        let aspect_ratio = match (scene.settings.aspect_height, scene.settings.aspect_width) {
            (Some(h), Some(w)) => w / h,
            (_, _) => scene.settings.width as f64 / scene.settings.height as f64,
        };

        Scene {
            origin: scene.settings.origin,
            look_at: scene.settings.look_at.unwrap_or(DVec3::NEG_Z),
            up_direction: scene.settings.up_direction.unwrap_or(DVec3::Y),
            height: scene.settings.height,
            width: scene.settings.width,
            vfov: scene.settings.vfov.unwrap_or(90.0),
            aspect_ratio,
            samples_per_pixel: scene.settings.samples_per_pixel.unwrap_or(100),
            recursive_depth: scene.settings.recursive_depth.unwrap_or(50),
            gamma: scene.settings.gamma.unwrap_or(2.2),
            world,
        }
    }

    /// Convert material to internal format
    fn material(
        settings: &scene_format::Settings,
        material: &scene_format::Material,
    ) -> Arc<dyn Material> {
        match material {
            scene_format::Material::Lambertian(mat) => Arc::new(Lambertian {
                albedo: mat.albedo,
                random_kind: settings.scattering_mode.unwrap_or_default(),
            }),
            scene_format::Material::Metallic(mat) => Arc::new(Metal {
                albedo: mat.albedo,
                fuzz: mat.fuzz,
                random_kind: settings.scattering_mode.unwrap_or_default(),
            }),
            scene_format::Material::Dielectric(mat) => Arc::new(Dielectric {
                refractive_index: mat.refractive_index,
            }),
        }
    }

    fn object(
        settings: &scene_format::Settings,
        object: scene_format::Object,
        materials: &HashMap<String, Arc<dyn Material>>,
    ) -> Box<dyn Hittable> {
        match object {
            scene_format::Object::Sphere(obj) => Box::new(Sphere {
                centre: obj.centre,
                radius: obj.radius,
                material: match obj.material {
                    scene_format::MaterialReference::Named(name) => Arc::clone(&materials[&name]),
                    scene_format::MaterialReference::Material(mat) => {
                        Scene::material(settings, &mat)
                    }
                },
            }),
        }
    }
}
