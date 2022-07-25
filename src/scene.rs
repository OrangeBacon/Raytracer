//! Convert from the scene file format into the internal rendering formats

use std::{collections::HashMap, sync::Arc};

use crate::{
    hit::{Hittable, Sphere},
    material::{Dielectric, Lambertian, Material, Metal},
    scene_format,
};

/// Scene description
#[derive(Debug)]
pub struct Scene {
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

    /// All the objects in the scene
    pub world: Vec<Box<dyn Hittable>>,
}

#[derive(Default)]
struct Defaults {
    materials: DefaultMaterials,
}

#[derive(Default)]
struct DefaultMaterials {
    lambertian: Lambertian,
    metallic: Metal,
    dielectric: Dielectric,
}

impl Scene {
    pub fn from_file(scene: scene_format::Scene) -> Scene {
        let defaults = Scene::defaults(&scene);
        let materials: HashMap<_, _> = scene
            .materials
            .into_iter()
            .map(|(k, v)| (k, Scene::material(v, &defaults, &Default::default())))
            .collect();

        let mut world = vec![];
        for object in scene.world {
            world.push(Scene::object(object, &defaults, &materials))
        }

        Scene {
            height: scene.settings.height,
            width: scene.settings.width,
            samples_per_pixel: scene.settings.samples_per_pixel,
            recursive_depth: scene.settings.recursive_depth,
            gamma: scene.settings.gamma,
            world,
        }
    }

    /// Get the defaults specified by the scene
    fn defaults(scene: &scene_format::Scene) -> Defaults {
        macro_rules! set {
            ($scene:expr, $($path:tt).+) => {
                if let Some(value) = $scene.$($path).+ {
                    $($path).+ = value;
                }
            };
        }

        let mut defaults = Defaults::default();
        set!(scene, defaults.materials.lambertian.albedo);
        set!(scene, defaults.materials.lambertian.random_kind);
        set!(scene, defaults.materials.metallic.albedo);
        set!(scene, defaults.materials.metallic.fuzz);
        set!(scene, defaults.materials.metallic.random_kind);
        set!(scene, defaults.materials.dielectric.refractive_index);

        defaults
    }

    /// Convert material to internal format
    fn material(
        material: scene_format::Material,
        defaults: &Defaults,
        materials: &HashMap<String, Arc<dyn Material>>,
    ) -> Arc<dyn Material> {
        match material {
            scene_format::Material::Reference { name } => Arc::clone(&materials[&name]),
            scene_format::Material::Lambertian(mat) => Arc::new(Lambertian {
                albedo: mat.albedo.unwrap_or(defaults.materials.lambertian.albedo),
                random_kind: mat
                    .random_kind
                    .unwrap_or(defaults.materials.lambertian.random_kind),
            }),
            scene_format::Material::Metallic(mat) => Arc::new(Metal {
                albedo: mat.albedo.unwrap_or(defaults.materials.metallic.albedo),
                fuzz: mat.fuzz.unwrap_or(defaults.materials.metallic.fuzz),
                random_kind: mat
                    .random_kind
                    .unwrap_or(defaults.materials.metallic.random_kind),
            }),
            scene_format::Material::Dielectric(mat) => Arc::new(Dielectric {
                refractive_index: mat
                    .refractive_index
                    .unwrap_or(defaults.materials.dielectric.refractive_index),
            }),
        }
    }

    fn object(
        object: scene_format::Object,
        defaults: &Defaults,
        materials: &HashMap<String, Arc<dyn Material>>,
    ) -> Box<dyn Hittable> {
        match object {
            scene_format::Object::Sphere(obj) => Box::new(Sphere {
                centre: obj.centre,
                radius: obj.radius,
                material: Scene::material(obj.material, defaults, materials),
            }),
        }
    }
}
