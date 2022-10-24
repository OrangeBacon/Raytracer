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

    Ok(())
}
