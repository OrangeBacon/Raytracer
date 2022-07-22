mod ray;
mod ray_trace;
mod hit;

use std::path::PathBuf;

use anyhow::Result;
use clap::Parser;
use glam::{vec3, Vec3};
use hit::{Sphere, Hittable};
use ray::Ray;

#[derive(Debug, Parser)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Width in pixels of the generated image
    #[clap(short, long, value_parser = clap::value_parser!(u32).range(1..), default_value_t = 640)]
    width: u32,

    /// Height in pixels of the generated image
    #[clap(short, long, value_parser = clap::value_parser!(u32).range(1..), default_value_t = 480)]
    height: u32,

    #[clap(short, long, value_parser, default_value = "out.png")]
    output: PathBuf,
}

fn main() -> Result<()> {
    let args = Args::parse();

    let aspect_ratio = args.width as f32 / args.height as f32;
    let viewport_height = 2.0;
    let viewport_width = viewport_height * aspect_ratio;
    let focal_length = 1.0;

    let origin = glam::Vec3::ZERO;
    let horizontal = vec3(viewport_width, 0.0, 0.0);
    let vertical = vec3(0.0, viewport_height, 0.0);
    let lower_left = origin - horizontal / 2.0 - vertical / 2.0 - vec3(0.0, 0.0, focal_length);

    let world: Vec<Box<dyn Hittable>> = vec![
        Box::new(Sphere { centre: Vec3::NEG_Z, radius: 0.5 }),
        Box::new(Sphere { centre: vec3(0.0, -100.5, -1.0), radius: 100.0 }),
    ];

    ray_trace::ray_trace(args.width, args.height, |x, y| {
        let u = x as f32 / args.width as f32;
        let v = y as f32 / args.height as f32;

        let ray = Ray {
            origin,
            direction: lower_left + u * horizontal + v * vertical - origin,
        };

        ray_color(ray, &world[..])
    })
    .save(&args.output)?;

    Ok(())
}

fn ray_color(ray: Ray, world: &[Box<dyn Hittable>]) -> Vec3 {
    if let Some(record) = world.hit(ray, 0.0, f32::INFINITY) {
        return 0.5 * (record.normal + 1.0);
    }
    
    let direction = ray.direction.normalize();
    let t = 0.5 * (direction.y + 1.0);
    (1.0 - t) * Vec3::ONE + t * vec3(0.5, 0.7, 1.0)
}
