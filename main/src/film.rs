use std::path::PathBuf;

use colour::{filters::Filter, spectrum::RGBSpectrum};
use geometry::{Bounds2, Bounds2i, Number, Point2, Point2i, Vector2};

/// A surface to render an image onto
pub struct Film<T: Number> {
    /// The width and height of the surface in pixels
    full_resolution: Point2i,

    /// The physical size of the film, stored in metres
    diagonal: T,

    /// Filtering function in use
    filter: Box<dyn Filter<T>>,

    /// The path to write the file to
    file_name: PathBuf,

    /// pixel scale factor to multiply every color by once the image has been
    /// generated, useful for maximising dynamic range in 8 bit files
    scale: T,

    /// The bounds of the area that will actually be rendered to, not all of the
    /// film has to be used for rendering.
    crop_bounds: Bounds2i,

    /// The actual image data
    pixels: Box<[Pixel<T>]>,

    /// Pre-computed output values from the filter
    filter_table: [T; FILTER_TABLE_WIDTH * FILTER_TABLE_WIDTH],
}

/// Size of pre-computed filter evaluations
const FILTER_TABLE_WIDTH: usize = 16;

/// A single pixel in the final image
#[derive(Debug, Clone, Copy, Default)]
struct Pixel<T: Number> {
    xyz: [T; 3],
    filter_weight_sum: T,
    splat_xyz: [T; 3],
    _padding: T,
}
const _: () = assert!(std::mem::size_of::<Pixel<f32>>() == 32);

/// A small region of an image that can be rendered to
pub struct FilmTile<'a, T: Number> {
    bounds: Bounds2i,
    filter_radius: Vector2<T>,
    inv_filter_radius: Vector2<T>,
    filter_table: &'a [T],
    pixels: Vec<FilmTilePixel<T>>,
}

struct FilmTilePixel<T: Number> {
    contribution_sum: RGBSpectrum<T>,
    filter_weight_sum: T,
}

impl<T: Number> Film<T> {
    /// Create a new film area to render to.
    /// Parameters:
    /// - Resolution: pixel width and height of the image
    /// - crop window: how much of the film is going to be rendered to, specified
    ///   in normalised coordinates from (0,0) top left to (1,1) bottom right.
    /// - filter: The filter used when combining pixel samples
    /// - diagonal: How large the film is, distance specified in mm.
    /// - file name: The path to save the final image to.
    /// - scale: multiplication factor to use once the image has been created
    pub fn new<F: Filter<T> + 'static>(
        resolution: Point2i,
        crop_window: Bounds2<T>,
        filter: F,
        diagonal: T,
        file_name: impl Into<PathBuf>,
        scale: T,
    ) -> Self {
        let crop_bounds = Bounds2::new(
            Point2i::new(
                (T::cast(resolution.x) * crop_window.min.x).ceil().i32(),
                (T::cast(resolution.y) * crop_window.min.y).ceil().i32(),
            ),
            Point2i::new(
                (T::cast(resolution.x) * crop_window.max.x).ceil().i32(),
                (T::cast(resolution.y) * crop_window.max.y).ceil().i32(),
            ),
        );

        let mut filter_table = [T::ZERO; FILTER_TABLE_WIDTH * FILTER_TABLE_WIDTH];
        let mut offset = 0;
        for y in 0..FILTER_TABLE_WIDTH {
            for x in 0..FILTER_TABLE_WIDTH {
                let point = Point2::new(
                    (T::cast(x as i64) + T::HALF) * filter.radius().x
                        / T::cast(FILTER_TABLE_WIDTH as i64),
                    (T::cast(y as i64) + T::HALF) * filter.radius().y
                        / T::cast(FILTER_TABLE_WIDTH as i64),
                );
                filter_table[offset] = filter.eval(point);
                offset += 1;
            }
        }

        Self {
            full_resolution: resolution,
            diagonal: diagonal * T::cast(0.001),
            filter: Box::new(filter),
            file_name: file_name.into(),
            scale,
            crop_bounds,
            pixels: vec![Pixel::default(); crop_bounds.surface_area() as _].into_boxed_slice(),
            filter_table,
        }
    }

    /// Get the bounds of a single sample in the film
    pub fn sample_bounds(&self) -> Bounds2i {
        let half = Vector2::new(T::HALF, T::HALF);

        Bounds2::new(
            (self.crop_bounds.min.cast() + half - self.filter.radius()).floor(),
            (self.crop_bounds.max.cast() - half + self.filter.radius()).ceil(),
        )
        .cast()
    }

    /// Get the size of the image in the scene
    pub fn physical_extent(&self) -> Bounds2<T> {
        let aspect = T::cast(self.full_resolution.y) / T::cast(self.full_resolution.x);
        let x = (self.diagonal * self.diagonal / (T::ONE + aspect * aspect)).sqrt();
        let y = aspect * x;
        Bounds2::new(
            Point2::new(-x / T::TWO, -y / T::TWO),
            Point2::new(x / T::TWO, y / T::TWO),
        )
    }

    /// Get a film tile covering the given bounds
    pub fn film_tile(&self, bounds: Bounds2i) -> FilmTile<T> {
        let half = Vector2::new(T::HALF, T::HALF);
        let bounds = bounds.cast();
        let p1 = (bounds.min - half - self.filter.radius()).ceil();
        let p2 = (bounds.max - half + self.filter.radius()).floor() + Point2::new(T::ONE, T::ONE);

        let tile_bounds = Bounds2::new(p1.cast(), p2.cast()).intersect(self.crop_bounds);

        FilmTile::new(tile_bounds, self.filter.radius(), &self.filter_table)
    }
}

impl<'a, T: Number> FilmTile<'a, T> {
    fn new(bounds: Bounds2i, filter_radius: Vector2<T>, filter_table: &'a [T]) -> Self {
        Self {
            bounds,
            filter_radius,
            inv_filter_radius: Vector2::new(T::ONE / filter_radius.x, T::ONE / filter_radius.y),
            filter_table,
            pixels: Vec::with_capacity(bounds.surface_area().max(0) as _),
        }
    }

    /// Add a radiance value to the film tile at a given location
    pub fn add_sample(&mut self, point: Point2<T>, radiance: RGBSpectrum<T>) {
        self.add_sample_weighted(point, radiance, T::ONE);
    }

    /// Add a radiance value to the film tile at a given location where the value
    /// has a given weight
    pub fn add_sample_weighted(&mut self, point: Point2<T>, radiance: RGBSpectrum<T>, weight: T) {
        let discrete_point = point - Vector2::splat(T::HALF);
        let p0 = (discrete_point - self.filter_radius).cast();
        let p1 = (discrete_point + self.filter_radius).cast() + Point2::new(1, 1);
        let p0 = p0.max(self.bounds.min);
        let p1 = p1.min(self.bounds.max);

        for y in p0.y..p1.y {
            for x in p0.x..p1.x {}
        }
    }
}
