use std::{iter::FromIterator, path::PathBuf};

use clap::Parser;
use geometry::{Number, Rng};
use parsers::{
    pbrt::PbrtFile,
    png::{Png, PngCreateInfo},
};

#[derive(Debug, Parser)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Uses 64-bit double precision floats instead of 32 bit single precision.
    #[clap(short = 'd', long = "double", value_parser, default_value = "false")]
    use_double: bool,

    /// Print the internal representation from all specified components
    #[clap(short, long, value_parser, value_delimiter = ',')]
    print: Vec<String>,

    /// File name to write the image to.
    #[clap(short, long, value_parser, default_value = "out.png")]
    output: PathBuf,

    /// File containing the scene description.
    #[clap(value_parser)]
    scene: PathBuf,
}

type Result<T = ()> = std::result::Result<T, Box<dyn std::error::Error>>;

fn main() -> Result {
    let args = Args::parse();

    if args.scene.extension().map(|s| s == "pbrt").unwrap_or(false) {
        if args.use_double {
            run_pbrt::<f64>(&args)?;
        } else {
            run_pbrt::<f32>(&args)?;
        }
    } else {
        panic!("Unknown file type")
    }

    Ok(())
}

fn run_pbrt<T: Number>(args: &Args) -> Result {
    let file: PbrtFile<T> = PbrtFile::parse(&args.scene)?;

    if args.print.iter().any(|x| x == "pbrt_ast") {
        println!("{:#?}", file)
    }

    // generate some random noise because haven't finished implementing the raytracer
    let mut rng = Rng::new(0);
    let image = Vec::from_iter(
        std::iter::from_fn(|| Some(rng.uniform_u32_limit(256) as _))
            .take(file.film.x_resolution as usize * file.film.y_resolution as usize * 3),
    );
    let image = Png::new(PngCreateInfo {
        image_data: image,
        width: file.film.x_resolution as usize,
        height: file.film.y_resolution as usize,
        colour_type: parsers::png::InputColourType::RGB,
        u16: false,
    })?;

    // write the noise out as png
    let mut name = file.film.file_name;
    name.set_extension("png");
    std::fs::write(name, image.write())?;

    Ok(())
}
