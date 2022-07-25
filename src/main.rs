mod camera;
mod hit;
mod material;
mod ray;
mod ray_trace;
mod scene;
mod scene_format;
mod utils;

use std::path::PathBuf;

use anyhow::Result;
use clap::Parser;
use glam::{dvec3, DVec3};
use hit::Hittable;
use rand::Rng;
use ray::Ray;

use crate::{camera::Camera, scene::Scene};

#[derive(Debug, Parser)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// File name to write the image to
    #[clap(short, long, value_parser, default_value = "out.png")]
    output: PathBuf,

    /// File containing the scene description
    #[clap(short, long, value_parser)]
    scene: PathBuf,
}

fn main() -> Result<()> {
    let args = Args::parse();

    let scene = std::fs::read_to_string(&args.scene)?;
    let scene = Scene::from_file(toml::from_str(&scene)?);

    let camera = Camera::new(scene.width, scene.height);

    ray_trace::ray_trace(scene.width, scene.height, scene.gamma, |x, y| {
        let mut colour = DVec3::ZERO;
        let mut rng = rand::thread_rng();
        for _ in 0..scene.samples_per_pixel {
            let u = (x as f64 + rng.gen::<f64>()) / scene.width as f64;
            let v = (y as f64 + rng.gen::<f64>()) / scene.height as f64;

            let ray = camera.ray(u, v);
            colour += ray_colour(ray, &scene.world[..], scene.recursive_depth);
        }

        colour / scene.samples_per_pixel as f64
    })
    .save(&args.output)?;

    Ok(())
}

fn ray_colour(ray: Ray, world: &[Box<dyn Hittable>], depth: u32) -> DVec3 {
    if depth <= 0 {
        return DVec3::ZERO;
    }

    if let Some(record) = world.hit(ray, 0.001, f64::INFINITY) {
        if let Some(scatter) = record.material.scatter(ray, &record) {
            return scatter.attenuation * ray_colour(scatter.scattered, world, depth - 1);
        }
        return DVec3::ZERO;
    }

    let direction = ray.direction.normalize();
    let t = 0.5 * (direction.y + 1.0);
    (1.0 - t) * DVec3::ONE + t * dvec3(0.5, 0.7, 1.0)
}
