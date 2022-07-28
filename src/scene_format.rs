//! Serde description of the scene data structure as read by this program.
//! Should not be used for rendering, should be converted to internal formats

use std::collections::HashMap;

use glam::DVec3;
use serde::Deserialize;

/// Top level scene file
#[derive(Deserialize, Debug, Clone)]
pub struct Scene {
    /// Required settings to render a scene
    pub settings: Settings,

    /// Globally available materials
    #[serde(default)]
    pub materials: HashMap<String, Material>,

    #[serde(default)]
    pub textures: HashMap<String, Texture>,

    /// All objects present in the world
    #[serde(default)]
    pub world: Vec<Object>,
}

/// Required settings to render a scene
#[derive(Deserialize, Debug, Clone, Copy, Default)]
#[serde(default)]
pub struct Settings {
    /// Vertical field of view in degrees
    pub vfov: Option<f64>,

    /// Aspect ratio numerator
    pub aspect_width: Option<f64>,

    /// Aspect ratio denominator
    pub aspect_height: Option<f64>,

    /// Height in pixels of the output image
    pub height: u32,

    /// Width in pixels of the output image
    pub width: u32,

    /// The number of rays sent for each pixel
    pub samples_per_pixel: Option<u32>,

    /// The maximum recursive depth of each ray
    pub recursive_depth: Option<u32>,

    /// output = pow(output, 1.0 / gamma), gamma correction/tone mapping
    pub gamma: Option<f64>,

    /// The random number generation to use during scattering
    pub scattering_mode: Option<RandomKind>,

    /// Where the camera is located
    pub origin: DVec3,

    /// Direction the camera is looking at
    pub look_at: Option<DVec3>,

    /// The vertical upwards direction of the camera
    pub up_direction: Option<DVec3>,

    /// The size of the camera's aperture
    pub aperture: Option<f64>,

    /// The distance of the focal plane from the origin
    /// If not specified, no depth of field effect will be used.
    pub focus_distance: Option<f64>,

    /// Start time of the camera aperture opening
    pub time0: f64,

    /// End time of the camera aperture opening
    pub time1: Option<f64>,
}

/// All possible materials to use for an object
#[derive(Deserialize, Debug, Clone)]
#[serde(tag = "kind")]
pub enum Material {
    Lambertian(Lambertian),
    Metallic(Metal),
    Dielectric(Dielectric),
}

/// Lambertian (diffuse) material
#[derive(Deserialize, Debug, Default, Clone)]
#[serde(default)]
pub struct Lambertian {
    /// The base colour of the material
    pub albedo: TextureReference,
}

/// Metallic materials
#[derive(Deserialize, Debug, Default, Clone)]
#[serde(default)]
pub struct Metal {
    /// The base colour of the material
    pub albedo: TextureReference,

    /// How rough the metallic object is
    pub fuzz: f64,
}

/// Dielectric Material
#[derive(Deserialize, Debug, Default, Clone, Copy)]
#[serde(default)]
pub struct Dielectric {
    pub refractive_index: f64,
}

/// The random number generation to use during scattering
#[derive(Deserialize, Debug, Default, Clone, Copy)]
pub enum RandomKind {
    /// Scatter in all directions around a collision
    #[default]
    Sphere,

    /// Only scatter in the hemisphere facing outwards from a collision
    Hemisphere,
}

#[derive(Deserialize, Debug, Clone)]
#[serde(untagged)]
pub enum TextureReference {
    Solid(DVec3),
    Named(String),
    Texture(Texture),
}
impl Default for TextureReference {
    fn default() -> Self {
        Self::Solid(DVec3::default())
    }
}

#[derive(Deserialize, Debug, Clone)]
#[serde(tag = "kind")]
pub enum Texture {
    Solid {
        albedo: DVec3,
    },
    Checker {
        odd: Box<Texture>,
        even: Box<Texture>,
    },
}

#[derive(Deserialize, Debug, Clone)]
#[serde(tag = "shape")]
pub enum Object {
    Sphere(Sphere),
    SphereField,
}

/// A Single sphere
#[derive(Deserialize, Debug, Clone)]
pub struct Sphere {
    /// Radius of the sphere
    pub radius: f64,

    /// The centre location
    pub centre: DVec3,

    /// Material the object is made from
    pub material: MaterialReference,
}

#[derive(Deserialize, Debug, Clone)]
#[serde(untagged)]
pub enum MaterialReference {
    Named(String),
    Material(Material),
}
