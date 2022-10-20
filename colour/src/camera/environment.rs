use geometry::{lerp, AnimatedTransform, ConstZero, Number, Point3, Ray, Vector3};

use crate::film::Film;

use super::{Camera, CameraData};

/// 360 degree camera for recording all light reaching a point in the scene
pub struct Environment<'a, T: Number> {
    camera_data: CameraData<'a, T>,
}

impl<'a, T: Number> Environment<'a, T> {
    /// Create a new environment camera
    pub fn new(
        camera_to_world: AnimatedTransform<T>,
        shutter_open: T,
        shutter_close: T,
        film: &'a Film<T>,
    ) -> Self {
        Self {
            camera_data: CameraData {
                camera_to_world,
                shutter_open,
                shutter_close,
                film,
            },
        }
    }
}

impl<'a, T: Number> Camera<'a, T> for Environment<'a, T> {
    fn generate_ray(&self, sample: super::CameraSample<T>) -> (T, Ray<T>) {
        let film_resolution = self.camera_data.film.full_resolution().cast();

        let theta = T::PI * sample.film.y / film_resolution.y;
        let phi = T::TWO * T::PI * sample.film.x / film_resolution.x;
        let dir = Vector3::new(
            theta.sin() * phi.cos(),
            theta.cos(),
            theta.sin() * phi.sin(),
        );

        let mut ray = Ray::new(Point3::ZERO, dir);
        ray.t_max = T::INFINITY;
        ray.time = lerp(
            sample.time,
            self.camera_data.shutter_open,
            self.camera_data.shutter_close,
        );

        (T::ONE, ray * self.camera_data.camera_to_world)
    }

    fn camera_data(&self) -> &'a CameraData<T> {
        &self.camera_data
    }

    fn camera_data_mut(&mut self) -> &'a mut CameraData<T> {
        &mut self.camera_data
    }
}
