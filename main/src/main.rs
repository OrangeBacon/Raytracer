mod aabb;
mod bvh;
mod camera;
mod material;
mod pbrt_parser;
mod perlin;
mod ray;
mod ray_trace;
mod scene;
mod scene_format;
mod shapes;
mod texture;
mod utils;

use std::{path::PathBuf, sync::Arc};

use anyhow::Result;
use clap::Parser;
use glam::DVec3;
use rand::Rng;

use crate::{camera::Camera, pbrt_parser::parse_pbrt, ray::Ray, scene::Scene, shapes::Hittable};

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

    if args.scene.extension().map(|s| s == "pbrt").unwrap_or(false) {
        parse_pbrt(&args.scene)?;
        return Ok(());
    }

    let scene = std::fs::read_to_string(&args.scene)?;
    let scene = Scene::from_file(toml::from_str(&scene)?)?;

    let camera = Camera::new(&scene);

    ray_trace::ray_trace(scene.width, scene.height, scene.gamma, |x, y| {
        let mut colour = DVec3::ZERO;
        let mut rng = rand::thread_rng();
        for _ in 0..scene.samples_per_pixel {
            let u = (x as f64 + rng.gen::<f64>()) / scene.width as f64;
            let v = (y as f64 + rng.gen::<f64>()) / scene.height as f64;

            let ray = camera.ray(u, v);
            colour += ray_colour(
                ray,
                scene.background_colour,
                &scene.world,
                scene.recursive_depth,
            );
        }

        colour / scene.samples_per_pixel as f64
    })
    .save(&args.output)?;

    Ok(())
}

fn ray_colour(ray: Ray, background: DVec3, world: &Arc<dyn Hittable>, depth: u32) -> DVec3 {
    if depth == 0 {
        return DVec3::ZERO;
    }

    if let Some(record) = world.hit(ray, 0.001, f64::INFINITY) {
        let emitted = record.material.emitted(record.uv, record.point);

        if let Some(scatter) = record.material.scatter(ray, &record) {
            return emitted
                + scatter.attenuation
                    * ray_colour(scatter.scattered, background, world, depth - 1);
        }
        return emitted;
    }

    background
}
