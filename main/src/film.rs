use std::{
    array,
    path::PathBuf,
    sync::{Arc, Mutex},
};

use colour::{
    filters::Filter,
    spectrum::{xyz_to_rgb, RGBSpectrum},
};
use geometry::{AtomicNumber, Bounds2, Bounds2i, Number, Point2, Point2i, Vector2};

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
    pixels: Arc<Mutex<Vec<Pixel<T>>>>,

    /// Pre-computed output values from the filter
    filter_table: [T; FILTER_TABLE_WIDTH * FILTER_TABLE_WIDTH],
}

/// Size of pre-computed filter evaluations
const FILTER_TABLE_WIDTH: usize = 16;

/// A single pixel in the final image
#[derive(Debug, Clone, Copy, Default)]
pub struct Pixel<T: Number> {
    xyz: [T; 3],
    filter_weight_sum: T,
    splat_xyz: [T::Atomic; 3],
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

pub struct FilmTilePixel<T: Number> {
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

        let mut pixels = Vec::with_capacity(crop_bounds.surface_area() as _);
        pixels.fill_with(Default::default);
        let pixels = Arc::new(Mutex::new(pixels));

        Self {
            full_resolution: resolution,
            diagonal: diagonal * T::cast(0.001),
            filter: Box::new(filter),
            file_name: file_name.into(),
            scale,
            crop_bounds,
            pixels,
            filter_table,
        }
    }

    /// Set all pixels in the image to the value stored in spectrum
    fn set_image(&mut self, spectrum: RGBSpectrum<T>) {
        let mut pixels = self.pixels.lock().unwrap();

        for pixel in pixels.iter_mut() {
            pixel.xyz = spectrum.to_xyz();
            pixel.filter_weight_sum = T::ONE;
            pixel.splat_xyz = array::from_fn(|_| T::Atomic::from_float(T::ZERO));
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

    /// Add a film tile to this film
    pub fn merge_film_tile(&self, tile: FilmTile<T>) {
        let mut pixels = self.pixels.lock().unwrap();

        for pixel in tile.bounds {
            let tile_pixel = tile.pixel(pixel);
            let merge_pixel = self.pixel_mut(&mut pixels, pixel);
            let xyz = tile_pixel.contribution_sum.to_xyz();
            for (&xyz, merge) in xyz.iter().zip(merge_pixel.xyz.iter_mut()) {
                *merge += xyz;
            }
            merge_pixel.filter_weight_sum += tile_pixel.filter_weight_sum;
        }
    }

    /// Get a reference to a pixel within the referenced unlocked pixel
    /// buffer from calling `self.pixels.lock()`
    fn pixel<'a>(&self, pixels: &'a [Pixel<T>], p: Point2i) -> &'a Pixel<T> {
        let width = self.crop_bounds.max.x - self.crop_bounds.min.x;
        let offset = (p.x - self.crop_bounds.min.x) + (p.y - self.crop_bounds.min.y) * width;

        &pixels[offset as usize]
    }

    /// Get a mutable reference to a pixel within the referenced unlocked pixel
    /// buffer from calling `self.pixels.lock()`
    fn pixel_mut<'a>(&self, pixels: &'a mut [Pixel<T>], p: Point2i) -> &'a mut Pixel<T> {
        let width = self.crop_bounds.max.x - self.crop_bounds.min.x;
        let offset = (p.x - self.crop_bounds.min.x) + (p.y - self.crop_bounds.min.y) * width;

        &mut pixels[offset as usize]
    }

    /// Add a contribution to an arbitrary pixel without using weighted sums
    pub fn add_splat(&self, point: Point2i, v: RGBSpectrum<T>) {
        if !self.crop_bounds.inside_exclusive(point) {
            return;
        }

        let xyz = v.to_xyz();
        let mut pixels = self.pixels.lock().unwrap();
        let pixel = self.pixel_mut(&mut pixels, point);
        for (&xyz, splat) in xyz.iter().zip(&pixel.splat_xyz) {
            splat.add(xyz);
        }
    }

    /// Write the image data stored in the film out to a file
    pub fn write_image(&self, splat_scale: T) {
        let mut image = Vec::with_capacity(3 * self.crop_bounds.surface_area() as usize);
        let pixels = self.pixels.lock().unwrap();
        for point in self.crop_bounds {
            let pixel = self.pixel(&pixels, point);
            let mut rgb = xyz_to_rgb(pixel.xyz);

            // normalise pixel value
            let filter_weight_sum = pixel.filter_weight_sum;
            if filter_weight_sum != T::ZERO {
                let inv = T::ONE / filter_weight_sum;
                rgb[0] = T::ZERO.max(rgb[0] * inv);
                rgb[1] = T::ZERO.max(rgb[1] * inv);
                rgb[2] = T::ZERO.max(rgb[2] * inv);
            }

            // add pixel splat values
            let splat = [
                pixel.splat_xyz[0].to_float(),
                pixel.splat_xyz[1].to_float(),
                pixel.splat_xyz[2].to_float(),
            ];
            let splat = xyz_to_rgb(splat);

            rgb[0] += splat_scale * splat[0];
            rgb[1] += splat_scale * splat[1];
            rgb[2] += splat_scale * splat[2];

            rgb[0] *= self.scale;
            rgb[1] *= self.scale;
            rgb[2] *= self.scale;

            image.extend_from_slice(&rgb);
        }

        todo!("Write image data to file: {:?}", image);
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
    pub fn add_sample_weighted(
        &mut self,
        point: Point2<T>,
        radiance: RGBSpectrum<T>,
        sample_weight: T,
    ) {
        let discrete_point = point - Vector2::splat(T::HALF);
        let p0 = (discrete_point - self.filter_radius).cast();
        let p1 = (discrete_point + self.filter_radius).cast() + Point2::new(1, 1);
        let p0 = p0.max(self.bounds.min);
        let p1 = p1.min(self.bounds.max);

        let mut ifx = vec![0; (p1.x - p0.x) as _];
        for x in p0.x..p1.x {
            let fx = ((T::cast(x) - discrete_point.x)
                * self.inv_filter_radius.x
                * T::cast(self.filter_table.len() as i64))
            .abs();
            ifx[(x - p0.x) as usize] = (fx.floor().i32()).min(self.filter_table.len() as _);
        }

        let mut ify = vec![0; (p1.y - p0.y) as _];
        for y in p0.y..p1.y {
            let fy = ((T::cast(y) - discrete_point.y)
                * self.inv_filter_radius.y
                * T::cast(self.filter_table.len() as i64))
            .abs();
            ify[(y - p0.y) as usize] = (fy.floor().i32()).min(self.filter_table.len() as _);
        }

        for y in p0.y..p1.y {
            for x in p0.x..p1.x {
                let offset = ify[(y - p0.y) as usize] * self.filter_table.len() as i32
                    + ifx[(x - p0.x) as usize];
                let filter_weight = self.filter_table[offset as usize];

                let pixel = self.pixel_mut(Point2i::new(x, y));
                pixel.contribution_sum += radiance * sample_weight * filter_weight;
                pixel.filter_weight_sum += filter_weight;
            }
        }
    }

    /// Get a pixel at a given location
    pub fn pixel(&self, p: Point2i) -> &FilmTilePixel<T> {
        let width = self.bounds.max.x - self.bounds.min.x;
        let offset = (p.x - self.bounds.min.x) + (p.y - self.bounds.min.y) * width;

        &self.pixels[offset as usize]
    }

    /// Get a pixel at a given location
    pub fn pixel_mut(&mut self, p: Point2i) -> &mut FilmTilePixel<T> {
        let width = self.bounds.max.x - self.bounds.min.x;
        let offset = (p.x - self.bounds.min.x) + (p.y - self.bounds.min.y) * width;

        &mut self.pixels[offset as usize]
    }
}
