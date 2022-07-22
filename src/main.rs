use std::path::PathBuf;

use anyhow::Result;
use clap::Parser;
use image::{ImageBuffer, Rgb};
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressDrawTarget};
use rayon::prelude::*;

#[derive(Debug, Parser)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Width in pixels of the generated image
    #[clap(short, long, value_parser = clap::value_parser!(u32).range(1..), default_value_t = 256)]
    width: u32,

    /// Height in pixels of the generated image
    #[clap(short, long, value_parser = clap::value_parser!(u32).range(1..), default_value_t = 256)]
    height: u32,

    #[clap(short, long, value_parser, default_value = "out.png")]
    output: PathBuf,
}

fn main() -> Result<()> {
    let args = Args::parse();

    let progress = ProgressBar::with_draw_target(
        (args.width * args.height).into(),
        ProgressDrawTarget::stderr_with_hz(5),
    );

    let pixels: Vec<_> = (0..args.height)
        .into_par_iter()
        .flat_map(|y| (0..args.width).into_par_iter().map(move |x| (x, y)))
        .progress_with(progress)
        .flat_map(|(x, y)| {
            let color = glam::vec3(
                x as f32 / args.width as f32,
                y as f32 / args.height as f32,
                0.25,
            );

            let [r,g,b] = (color * 255.99).to_array();
            [r as u8, g as u8, b as u8]
        })
        .collect();

    let mut image = ImageBuffer::<Rgb<_>, _>::from_raw(args.width, args.height, pixels).unwrap();

    image::imageops::flip_vertical_in_place(&mut image);
    image.save(&args.output)?;

    Ok(())
}
