use std::path::PathBuf;

use anyhow::Result;
use clap::Parser;

use parsers::parse_pbrt;

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

    Ok(())
}
