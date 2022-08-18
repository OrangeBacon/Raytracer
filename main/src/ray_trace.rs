use std::time::Instant;

use image::{ImageBuffer, Rgb, RgbImage};
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressDrawTarget};
use rayon::prelude::*;

/// Create a ray-traced image in parallel and display the progress
/// to the terminal.  Assumes the shading function returns colours as
/// vec3(r, g, b) where all components 0..=1
pub fn ray_trace(
    width: u32,
    height: u32,
    gamma: f64,
    f: impl Send + Sync + Fn(u32, u32) -> glam::DVec3,
) -> RgbImage {
    let progress = ProgressBar::with_draw_target(
        (height * width).into(),
        ProgressDrawTarget::stderr_with_hz(5),
    );

    let now = Instant::now();
    let pixels: Vec<_> = (0..height)
        .into_par_iter()
        .flat_map(|y| (0..width).into_par_iter().map(move |x| (x, y)))
        .progress_with(progress)
        .flat_map(|(x, y)| {
            let [r, g, b] = (f(x, y).powf(1.0 / gamma) * 255.99).to_array();
            [r as u8, g as u8, b as u8]
        })
        .collect();

    let mut image = ImageBuffer::<Rgb<_>, _>::from_raw(width, height, pixels).unwrap();

    image::imageops::flip_vertical_in_place(&mut image);

    let time = now.elapsed();
    print!("Image rendered: ");
    let min = time.as_secs() / 60;
    let hours = min / 60;

    match hours {
        2.. => print!("{} hours ", hours),
        1 => print!("{} hour ", hours),
        _ => (),
    }
    match min % 60 {
        min @ 2.. => print!("{} minutes ", min),
        1 => print!("{} minute ", min),
        _ => (),
    }

    print!("{:.2} seconds", time.as_secs_f64() % 60.0);

    image
}