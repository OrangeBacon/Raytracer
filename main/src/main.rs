use std::path::PathBuf;

use clap::Parser;
use geometry::Number;
use parsers::pbrt::PbrtFile;

#[derive(Debug, Parser)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Uses 64-bit double precision floats instead of 32 bit single precision.
    #[clap(short = 'd', long = "double", value_parser, default_value = "false")]
    use_double: bool,

    /// Print the internal representation of all parsed input files.
    #[clap(short, long, value_parser, default_value = "false")]
    print_files: bool,

    /// File name to write the image to.
    #[clap(short, long, value_parser, default_value = "out.png")]
    output: PathBuf,

    /// File containing the scene description.
    #[clap(value_parser)]
    scene: PathBuf,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    if args.scene.extension().map(|s| s == "pbrt").unwrap_or(false) {
        if args.use_double {
            run_pbrt(&args, PbrtFile::<f64>::parse(&args.scene)?);
        } else {
            run_pbrt(&args, PbrtFile::<f32>::parse(&args.scene)?);
        }
    } else {
        panic!("Unknown file type")
    }

    Ok(())
}

fn run_pbrt<T: Number>(args: &Args, file: PbrtFile<T>) {
    if args.print_files {
        println!("{:#?}", file)
    }
}
