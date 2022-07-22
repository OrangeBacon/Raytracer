use image::{ImageBuffer, Rgb, RgbImage};
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressDrawTarget};
use rayon::prelude::*;

/// Create a raytraced image in parallel and display the progress
/// to the terminal.  Assumes the shading function returns colors as
/// vec3(r, g, b) where all components 0..=1
pub fn ray_trace(
    width: u32,
    height: u32,
    f: impl Send + Sync + Fn(u32, u32) -> glam::Vec3,
) -> RgbImage {
    let progress = ProgressBar::with_draw_target(
        height.into(),
        ProgressDrawTarget::stderr_with_hz(5),
    );

    let pixels: Vec<_> = (0..height)
        .into_par_iter()
        .progress_with(progress)
        .flat_map(|y| (0..width).into_par_iter().map(move |x| (x, y)))
        .flat_map(|(x, y)| {
            let [r, g, b] = (f(x, y) * 255.99).to_array();
            [r as u8, g as u8, b as u8]
        })
        .collect();

    let mut image = ImageBuffer::<Rgb<_>, _>::from_raw(width, height, pixels).unwrap();

    image::imageops::flip_vertical_in_place(&mut image);

    image
}
