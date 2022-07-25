//! Serde description of the scene data structure as read by this program.
//! Should not be used for rendering, should be converted to internal formats

use std::collections::HashMap;

use glam::DVec3;
use serde::Deserialize;

/// Top level scene file
#[derive(Deserialize, Debug)]
pub struct Scene {
    /// Required settings to render a scene
    pub settings: Settings,

    /// Values to override the defaults of some types
    #[serde(default)]
    pub defaults: Defaults,

    /// Globally available materials
    #[serde(default)]
    pub materials: HashMap<String, Material>,

    /// All objects present in the world
    #[serde(default)]
    pub world: Vec<Object>,
}

/// Required settings to render a scene
#[derive(Deserialize, Debug)]
pub struct Settings {
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
}

/// Values to override the defaults of some types
#[derive(Deserialize, Debug, Default)]
#[serde(default)]
pub struct Defaults {
    /// Default material settings
    pub materials: DefaultMaterials,
}

/// Default material settings
#[derive(Deserialize, Debug, Default)]
#[serde(default)]
pub struct DefaultMaterials {
    /// Default lambertian material settings
    pub lambertian: Lambertian,

    /// Default metallic material settings
    pub metallic: Metal,

    /// Default dielectric material settings
    pub dielectric: Dielectric,
}

/// All possible materials to use for an object
#[derive(Deserialize, Debug, Clone)]
#[serde(tag = "kind")]
pub enum Material {
    Reference { name: String },
    Lambertian(Lambertian),
    Metallic(Metal),
    Dielectric(Dielectric),
}

/// Lambertian (diffuse) material
#[derive(Deserialize, Debug, Default, Clone, Copy)]
#[serde(default)]
pub struct Lambertian {
    /// The base colour of the material
    pub albedo: Option<DVec3>,

    /// The random number generation to use during scattering
    pub random_kind: Option<RandomKind>,
}

/// Metallic materials
#[derive(Deserialize, Debug, Default, Clone, Copy)]
#[serde(default)]
pub struct Metal {
    /// The base colour of the material
    pub albedo: Option<DVec3>,

    /// How rough the metallic object is
    pub fuzz: Option<f64>,

    /// The random number generation to use during scattering
    pub random_kind: Option<RandomKind>,
}

/// Dielectric Material
#[derive(Deserialize, Debug, Default, Clone, Copy)]
#[serde(default)]
pub struct Dielectric {
    pub refractive_index: Option<f64>,
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

#[derive(Deserialize, Debug)]
#[serde(tag = "shape")]
pub enum Object {
    Sphere(Sphere),
}

/// A Single sphere
#[derive(Deserialize, Debug)]
pub struct Sphere {
    /// Radius of the sphere
    pub radius: f64,

    /// The centre location
    pub centre: DVec3,

    /// Material the object is made from
    pub material: Material,
}
