use std::path::{PathBuf, Path};

/// A collection of configuration settings for rendering with
pub struct Pbrt {
    options: Options,
}

/// Configuration settings to initialise pbrt with
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Default, Hash)]
pub struct Options {
    /// Number of threads to use while rendering
    pub thread_count: u32,

    /// Should this render be done quickly (or not, do it high quality)
    pub quick_render: bool,

    /// Print no information while rendering
    pub quiet: bool,

    /// Print extra debug information while rendering, will be slower.
    pub verbose: bool,

    /// File path to write the render to
    pub image_file: PathBuf,
}

impl Pbrt {
    /// Start a new pbrt instance
    pub fn new(options: Options) -> Self {
        Self {
            options: options
        }
    }

    /// Load a new pbrt file
    pub fn parse_pbrt(&mut self, file_name: impl AsRef<Path>) -> Result<(), ()> {
        Ok(())
    }

    /// Render the requested image
    pub fn render(self) -> Result<(), ()> {
        Ok(())
    }
}
