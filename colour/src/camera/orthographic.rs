use geometry::{
    concentric_sample_disk, lerp, AnimatedTransform, Bounds2, ConstZero, Number, Point3, Ray,
    RayDifferential, RayDifferentials, Transform, Vector3,
};

use crate::film::Film;

use super::{Camera, CameraData, CameraSample, Projective};

pub struct Orthographic<'a, T: Number> {
    projective: Projective<'a, T>,
    dx_camera: Vector3<T>,
    dy_camera: Vector3<T>,
}

impl<'a, T: Number> Orthographic<'a, T> {
    pub fn new(
        camera_to_world: AnimatedTransform<T>,
        screen_window: Bounds2<T>,
        shutter_open: T,
        shutter_close: T,
        lens_radius: T,
        focal_distance: T,
        film: &'a Film<T>,
    ) -> Self {
        let projective = Projective::new(
            camera_to_world,
            Transform::orthographic(T::ZERO, T::ONE),
            screen_window,
            shutter_open,
            shutter_close,
            lens_radius,
            focal_distance,
            film,
        );

        let dx_camera = projective
            .raster_to_camera
            .apply(Vector3::new(T::ONE, T::ZERO, T::ZERO));
        let dy_camera = projective
            .raster_to_camera
            .apply(Vector3::new(T::ONE, T::ZERO, T::ZERO));

        Self {
            projective,
            dx_camera,
            dy_camera,
        }
    }
}

impl<'a, T: Number> Camera<'a, T> for Orthographic<'a, T> {
    fn generate_ray(&self, sample: CameraSample<T>) -> (T, geometry::Ray<T>) {
        let film = Point3::new(sample.film.x, sample.film.y, T::ZERO);
        let camera = self.projective.raster_to_camera.apply(film);
        let mut ray = Ray::new(camera, Vector3::Z);

        // Modify ray for depth of field effect
        if self.projective.lens_radius > T::ZERO {
            let lens = concentric_sample_disk(sample.lens) * self.projective.lens_radius;

            // compute point on plane of focus
            let ft = self.projective.focal_distance / ray.direction.z;
            let focus = ray.at(ft);

            ray.origin = Point3::new(lens.x, lens.y, T::ZERO);
            ray.direction = (focus - ray.origin).normalise()
        }

        ray.time = lerp(
            sample.time,
            self.projective.camera_data.shutter_open,
            self.projective.camera_data.shutter_close,
        );

        (T::ONE, ray * self.projective.camera_data.camera_to_world)
    }

    fn generate_ray_differential(
        &self,
        sample: CameraSample<T>,
    ) -> (T, geometry::RayDifferential<T>) {
        let film = Point3::new(sample.film.x, sample.film.y, T::ZERO);
        let camera = self.projective.raster_to_camera.apply(film);
        let mut ray = RayDifferential::new(camera, Vector3::Z);

        // Modify ray for depth of field effect
        if self.projective.lens_radius > T::ZERO {
            let lens = concentric_sample_disk(sample.lens) * self.projective.lens_radius;

            // compute point on plane of focus
            let ft = self.projective.focal_distance / ray.direction.z;
            let focus = ray.at(ft);

            ray.main.origin = Point3::new(lens.x, lens.y, T::ZERO);
            ray.main.direction = (focus - ray.origin).normalise();

            let dx = (camera + self.dx_camera).to_vec().normalise();
            let ft = self.projective.focal_distance / dx.z;
            let focus = Point3::ZERO + (dx * ft);
            let rx_origin = Point3::new(lens.x, lens.y, T::ZERO);
            let rx_direction = (focus - rx_origin).normalise();

            let dy = (camera + self.dy_camera).to_vec().normalise();
            let ft = self.projective.focal_distance / dy.z;
            let focus = Point3::ZERO + (dy * ft);
            let ry_origin = Point3::new(lens.x, lens.y, T::ZERO);
            let ry_direction = (focus - ry_origin).normalise();

            ray.differentials = Some(RayDifferentials {
                rx_origin,
                ry_origin,
                rx_direction,
                ry_direction,
            });
        } else {
            ray.differentials = Some(RayDifferentials {
                rx_origin: ray.origin + self.dx_camera,
                ry_origin: ray.origin + self.dy_camera,
                rx_direction: ray.direction,
                ry_direction: ray.direction,
            });
        }

        ray.main.time = lerp(
            sample.time,
            self.projective.camera_data.shutter_open,
            self.projective.camera_data.shutter_close,
        );

        (T::ONE, ray * self.projective.camera_data.camera_to_world)
    }

    fn camera_data(&self) -> &'a CameraData<T> {
        &self.projective.camera_data
    }

    fn camera_data_mut(&mut self) -> &'a mut CameraData<T> {
        &mut self.projective.camera_data
    }
}
