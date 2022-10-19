mod orthographic;
mod projective;

pub use orthographic::Orthographic;
pub use projective::Projective;

use geometry::{AnimatedTransform, Number, Point2, Ray, RayDifferential, RayDifferentials};

use crate::film::Film;

/// Data that every camera should store
pub struct CameraData<'a, T: Number> {
    pub camera_to_world: AnimatedTransform<T>,
    pub shutter_open: T,
    pub shutter_close: T,
    pub film: &'a Film<T>,
}

#[derive(Debug, Clone, Copy)]
pub struct CameraSample<T: Number> {
    pub film: Point2<T>,
    pub lens: Point2<T>,
    pub time: T,
}

pub trait Camera<'a, T: Number> {
    /// Generate a ray for a given sample.  Must return normalised ray.
    fn generate_ray(&self, sample: CameraSample<T>) -> (T, Ray<T>);

    fn generate_ray_differential(&self, sample: CameraSample<T>) -> (T, RayDifferential<T>) {
        let (wt, ray) = self.generate_ray(sample);
        let mut diff = RayDifferential::default();

        let mut shift = sample;
        shift.film.x += T::ONE;
        let (wtx, rx) = self.generate_ray(shift);
        if wtx == T::ZERO {
            return (T::ZERO, RayDifferential::from_ray(ray));
        }

        shift.film.x -= T::ONE;
        shift.film.y += T::ONE;
        let (wty, ry) = self.generate_ray(shift);
        if wty == T::ZERO {
            return (T::ZERO, RayDifferential::from_ray(ray));
        }

        diff.differentials = Some(RayDifferentials {
            rx_origin: rx.origin,
            ry_origin: ry.origin,
            rx_direction: rx.direction,
            ry_direction: ry.direction,
        });

        (wt, diff)
    }

    fn camera_data(&self) -> &'a CameraData<T>;

    fn camera_data_mut(&mut self) -> &'a mut CameraData<T>;
}
