mod camera;
mod hit;
mod material;
mod ray;
mod ray_trace;
mod utils;
mod scene_format;

use std::{path::PathBuf, sync::Arc};

use anyhow::Result;
use camera::Camera;
use clap::Parser;
use glam::{dvec3, DVec3};
use hit::{Hittable, Sphere};
use material::{Lambertian, Metal};
use rand::Rng;
use ray::Ray;

use crate::scene_format::Scene;

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
    let scene: Scene = toml::from_str(&scene)?;
    println!("{scene:?}");
    /*
    let camera = Camera::new(args.width, args.height);

    ray_trace::ray_trace(args.width, args.height, |x, y| {
        let mut color = DVec3::ZERO;
        let mut rng = rand::thread_rng();
        for _ in 0..args.samples_per_pixel {
            let u = (x as f64 + rng.gen::<f64>()) / args.width as f64;
            let v = (y as f64 + rng.gen::<f64>()) / args.height as f64;

            let ray = camera.ray(u, v);
            color += ray_color(ray, &world[..], args.max_depth);
        }

        color / args.samples_per_pixel as f64
    })
    .save(&args.output)?;*/

    Ok(())
}

fn ray_color(ray: Ray, world: &[Box<dyn Hittable>], depth: u32) -> DVec3 {
    if depth <= 0 {
        return DVec3::ZERO;
    }

    if let Some(record) = world.hit(ray, 0.001, f64::INFINITY) {
        if let Some(scatter) = record.material.scatter(ray, &record) {
            return scatter.attenuation * ray_color(scatter.scattered, world, depth - 1);
        }
        return DVec3::ZERO;
    }

    let direction = ray.direction.normalize();
    let t = 0.5 * (direction.y + 1.0);
    (1.0 - t) * DVec3::ONE + t * dvec3(0.5, 0.7, 1.0)
}
