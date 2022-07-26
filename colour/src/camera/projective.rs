use geometry::{AnimatedTransform, Bounds2, Number, Transform, Vector3};

use crate::film::Film;

use super::CameraData;

pub struct Projective<'a, T: Number> {
    pub camera_data: CameraData<'a, T>,
    pub camera_to_screen: Transform<T>,
    pub raster_to_camera: Transform<T>,
    pub screen_to_raster: Transform<T>,
    pub raster_to_screen: Transform<T>,
    pub lens_radius: T,
    pub focal_distance: T,
}

impl<'a, T: Number> Projective<'a, T> {
    pub fn new(
        camera_to_world: AnimatedTransform<T>,
        camera_to_screen: Transform<T>,
        screen_window: Bounds2<T>,
        shutter_open: T,
        shutter_close: T,
        lens_radius: T,
        focal_distance: T,
        film: &'a Film<T>,
    ) -> Self {
        let film_resolution = film.full_resolution().cast();
        let screen_to_raster =
            Transform::scale(Vector3::new(film_resolution.x, film_resolution.y, T::ONE))
                * Transform::scale(Vector3::new(
                    T::ONE / (screen_window.max.x - screen_window.min.x),
                    T::ONE / (screen_window.min.y - screen_window.max.y),
                    T::ONE,
                ))
                * Transform::translation(Vector3::new(
                    -screen_window.min.x,
                    -screen_window.max.y,
                    T::ZERO,
                ));
        Self {
            camera_data: CameraData {
                camera_to_world,
                shutter_open,
                shutter_close,
                film,
            },
            screen_to_raster,
            camera_to_screen,
            raster_to_screen: screen_to_raster.inverse(),
            raster_to_camera: camera_to_screen.inverse() * screen_to_raster.inverse(),
            lens_radius,
            focal_distance,
        }
    }
}
