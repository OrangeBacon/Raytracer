use api::{Options, Pbrt};

fn main() -> Result<(), ()> {
    let options = Options {
        thread_count: 1,
        quick_render: true,
        quiet: false,
        verbose: true,
        image_file: "image.png".into(),
    };

    let mut renderer = Pbrt::new(options);

    for arg in std::env::args() {
        renderer.parse_pbrt(arg)?;
    }

    renderer.render();

    Ok(())
}
