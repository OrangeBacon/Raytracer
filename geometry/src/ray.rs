use std::ops::Deref;

use crate::{ConstZero, Float, Number, Point3, Vector3};

/// Ray with an origin and a direction from the origin
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Ray<T, F: Number> {
    pub origin: Point3<F>,
    pub direction: Vector3<F>,
    pub t_max: F,
    pub time: F,
    pub material: T,
}

impl<T: ConstZero, F: Number> Default for Ray<T, F> {
    fn default() -> Self {
        Self::ZERO
    }
}

impl<T: ConstZero, F: Number> ConstZero for Ray<T, F> {
    /// Ray from (0,0,0) towards (0,0,0)
    const ZERO: Self = Self {
        origin: Point3::ZERO,
        direction: Vector3::ZERO,
        t_max: F::INFINITY,
        time: F::ZERO,
        material: T::ZERO,
    };
}

impl<T: Default, F: Number> Ray<T, F> {
    /// Create a new ray with a given origin and direction
    pub fn new(origin: Point3<F>, direction: Vector3<F>) -> Self {
        Self {
            origin,
            direction,
            t_max: F::INFINITY,
            time: F::ZERO,
            material: Default::default(),
        }
    }
}

impl<T, F: Number> Ray<T, F> {
    /// Create a new ray with a given origin and direction
    pub fn new_with(origin: Point3<F>, direction: Vector3<F>, data: T) -> Self {
        Self {
            origin,
            direction,
            t_max: F::INFINITY,
            time: F::ZERO,
            material: data,
        }
    }

    /// Calculate a point at t distance along the ray
    pub fn at(&self, t: F) -> Point3<F> {
        self.origin + self.direction * t
    }
}

/// Information about two rays
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default, Hash)]
pub struct RayDifferentials<F: Number = Float> {
    pub rx_origin: Point3<F>,
    pub ry_origin: Point3<F>,
    pub rx_direction: Vector3<F>,
    pub ry_direction: Vector3<F>,
}

/// Ray with differential information
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct RayDifferential<T, F: Number> {
    pub main: Ray<T, F>,
    pub differentials: Option<RayDifferentials<F>>,
}

impl<T: ConstZero, F: Number> RayDifferential<T, F> {
    /// Ray::Zero, differentials = None
    pub const ZERO: Self = Self {
        main: Ray::ZERO,
        differentials: None,
    };
}

impl<T: Default, F: Number> RayDifferential<T, F> {
    /// Create a new differential ray with no recorded differentials
    pub fn new(origin: Point3<F>, direction: Vector3<F>) -> Self {
        Self {
            main: Ray::new(origin, direction),
            differentials: None,
        }
    }
}

impl<T, F: Number> RayDifferential<T, F> {
    /// Add differential storage to a given ray
    pub fn from_ray(ray: Ray<T, F>) -> Self {
        Self {
            main: ray,
            differentials: None,
        }
    }
}

impl<T: Copy, F: Number> RayDifferential<T, F> {
    /// Scale the differentials by a given spacing
    pub fn scale(&self, scale: F) -> Self {
        Self {
            main: self.main,
            differentials: self.differentials.map(|mut diff| {
                diff.rx_origin = self.origin + (diff.rx_origin - self.origin) * scale;
                diff.ry_origin = self.origin + (diff.ry_origin - self.origin) * scale;
                diff.rx_direction = self.direction + (diff.rx_direction - self.direction) * scale;
                diff.ry_direction = self.direction + (diff.ry_direction - self.direction) * scale;
                diff
            }),
        }
    }
}

impl<T, F: Number> Deref for RayDifferential<T, F> {
    type Target = Ray<T, F>;

    fn deref(&self) -> &Self::Target {
        &self.main
    }
}
