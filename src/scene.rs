//! Convert from the scene file format into the internal rendering formats

use std::{collections::HashMap, sync::Arc};

use glam::{dvec3, DVec3};
use rand::Rng;

use crate::{
    bvh::BVHNode,
    hit::{Hittable, MovingSphere, Sphere},
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

    /// The size of the aperture
    pub aperture: f64,

    /// The distance of the focal plane from the origin
    /// If not specified, no depth of field effect will be used.
    pub focal_distance: Option<f64>,

    /// All the objects in the scene
    pub world: Arc<dyn Hittable>,

    /// Start time of the camera aperture opening
    pub time0: f64,

    /// End time of the camera aperture opening
    pub time1: f64,
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
            aperture: scene.settings.aperture.unwrap_or(0.5),
            focal_distance: scene.settings.focus_distance,
            time0: scene.settings.time0,
            time1: scene.settings.time1.unwrap_or(1.0),
            world: BVHNode::new(
                &world,
                scene.settings.time0,
                scene.settings.time1.unwrap_or(1.0),
            )
            .unwrap(),
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
    ) -> Arc<dyn Hittable> {
        match object {
            scene_format::Object::Sphere(obj) => Arc::new(Sphere {
                centre: obj.centre,
                radius: obj.radius,
                material: match obj.material {
                    scene_format::MaterialReference::Named(name) => Arc::clone(&materials[&name]),
                    scene_format::MaterialReference::Material(mat) => {
                        Scene::material(settings, &mat)
                    }
                },
            }),
            scene_format::Object::SphereField => gen_sphere_field(),
        }
    }
}

fn gen_sphere_field() -> Arc<dyn Hittable> {
    let mut rng = rand::thread_rng();

    let mut objects: Vec<Arc<dyn Hittable>> = vec![];

    for a in -11..11 {
        for b in -11..11 {
            let radius = 0.2;
            let centre = dvec3(
                a as f64 + 0.9 * rng.gen::<f64>(),
                radius,
                b as f64 + 0.9 * rng.gen::<f64>(),
            );

            if (centre - dvec3(4.0, 0.2, 0.0)).length() <= 0.9 {
                continue;
            }

            let mat: f64 = rng.gen();

            if mat < 0.8 {
                let albedo = DVec3::from_array(rng.gen()) * DVec3::from_array(rng.gen());
                let centre2 = centre + dvec3(0.0, rng.gen_range(0.0..=0.5), 0.0);

                objects.push(Arc::new(MovingSphere {
                    centre0: centre,
                    centre1: centre2,
                    time0: 0.0,
                    time1: 1.0,
                    radius,
                    material: Arc::new(Lambertian {
                        albedo,
                        random_kind: scene_format::RandomKind::Hemisphere,
                    }),
                }));
            } else if mat < 0.95 {
                let albedo = dvec3(
                    rng.gen_range(0.5..1.0),
                    rng.gen_range(0.5..1.0),
                    rng.gen_range(0.5..1.0),
                );
                let fuzz = rng.gen_range(0.0..0.5);
                objects.push(Arc::new(Sphere {
                    centre,
                    radius,
                    material: Arc::new(Metal {
                        albedo,
                        fuzz,
                        random_kind: scene_format::RandomKind::Hemisphere,
                    }),
                }));
            } else {
                objects.push(Arc::new(Sphere {
                    centre,
                    radius,
                    material: Arc::new(Dielectric {
                        refractive_index: 1.5,
                    }),
                }))
            }
        }
    }

    BVHNode::new(&objects, 0.0, 1.0).unwrap()
}
