mod camera;
mod hit;
mod material;
mod ray;
mod ray_trace;

use std::{path::PathBuf, sync::Arc};

use anyhow::Result;
use camera::Camera;
use clap::Parser;
use glam::{dvec3, DVec3};
use hit::{Hittable, Sphere};
use material::{Lambertian, Metal};
use rand::Rng;
use ray::Ray;

#[derive(Debug, Parser)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Width in pixels of the generated image
    #[clap(short, long, value_parser = clap::value_parser!(u32).range(1..), default_value_t = 400)]
    width: u32,

    /// Height in pixels of the generated image
    #[clap(short, long, value_parser = clap::value_parser!(u32).range(1..), default_value_t = 225)]
    height: u32,

    /// File name to write the image to
    #[clap(short, long, value_parser, default_value = "out.png")]
    output: PathBuf,

    /// How many rays should be fired per pixel
    #[clap(short, long, value_parser = clap::value_parser!(u32).range(1..), default_value_t = 100)]
    samples_per_pixel: u32,

    /// Maximum recursive depth of every ray cast
    #[clap(short, long, value_parser = clap::value_parser!(u32).range(1..), default_value_t = 50)]
    max_depth: u32,
}

fn main() -> Result<()> {
    let args = Args::parse();

    let camera = Camera::new(args.width, args.height);

    let ground = Arc::new(Lambertian {
        albedo: dvec3(0.8, 0.8, 0.0),
    });
    let centre = Arc::new(Lambertian {
        albedo: dvec3(0.7, 0.3, 0.3),
    });
    let left = Arc::new(Metal {
        albedo: dvec3(0.8, 0.8, 0.8),
    });
    let right = Arc::new(Metal {
        albedo: dvec3(0.8, 0.6, 0.2),
    });

    let world: Vec<Box<dyn Hittable>> = vec![
        Box::new(Sphere {
            centre: dvec3(0.0, -100.5, -1.0),
            radius: 100.0,
            material: ground,
        }),
        Box::new(Sphere {
            centre: dvec3(0.0, 0.0, -1.0),
            radius: 0.5,
            material: centre,
        }),
        Box::new(Sphere {
            centre: dvec3(-1.0, 0.0, -1.0),
            radius: 0.5,
            material: left,
        }),
        Box::new(Sphere {
            centre: dvec3(1.0, 0.0, -1.0),
            radius: 0.5,
            material: right,
        }),
    ];

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
    .save(&args.output)?;

    Ok(())
}

pub fn rand_sphere_point() -> DVec3 {
    let mut rng = rand::thread_rng();

    let point = DVec3::from_array(rng.gen());

    point.length_recip() * point
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
