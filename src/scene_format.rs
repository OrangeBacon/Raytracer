//! Serde description of the scene data structure as read by this program.
//! Should not be used for rendering, should be converted to internal formats

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

    /// Globally avaliable materials
    #[serde(default)]
    pub materials: Vec<Material>,

    /// All objects present in the world
    #[serde(default)]
    pub world: Vec<Object>,
}

#[derive(Deserialize, Debug)]
pub struct Settings {
    pub height: u32,
    pub width: u32,
    pub samples_per_pixel: u32,
    pub recursive_depth: u32,
}

#[derive(Deserialize, Debug, Default)]
#[serde(default)]
pub struct Defaults {
    pub materials: DefaultMaterials,
}

#[derive(Deserialize, Debug, Default)]
#[serde(default)]
pub struct DefaultMaterials {
    pub lambertian: Lambertian,
    pub metal: Metal,
}

#[derive(Deserialize, Debug)]
#[serde(tag = "kind")]
pub enum Material {
    Lambertian(Lambertian),
    Metal(Metal),
}

#[derive(Deserialize, Debug, Default)]
#[serde(default)]
pub struct Lambertian {
    pub albedo: Option<DVec3>,
    pub random_kind: Option<RandomKind>,
}

#[derive(Deserialize, Debug, Default)]
#[serde(default)]
pub struct Metal {
    pub albedo: Option<DVec3>,
    pub fuzz: Option<f64>,
    pub random_kind: Option<RandomKind>,
}

#[derive(Deserialize, Debug)]
pub enum RandomKind {
    Hemisphere,
    Sphere,
}

impl Default for RandomKind {
    fn default() -> Self {
        RandomKind::Sphere
    }
}

#[derive(Deserialize, Debug)]
#[serde(tag = "shape")]
pub enum Object {
    Sphere(Sphere),
}

#[derive(Deserialize, Debug)]
pub struct Sphere {
    pub radius: f64,
    pub centre: DVec3,
    pub material: MaterialReference,
}

#[derive(Deserialize, Debug)]
#[serde(untagged)]
pub enum MaterialReference {
    Name(String),
    Material(Material),
}
