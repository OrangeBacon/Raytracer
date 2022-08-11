use std::ops::Deref;

use crate::{Float, Point3f, Vector3f, ConstZero};

/// Ray with an origin and a direction from the origin
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Ray<T = ()> {
    pub origin: Point3f,
    pub direction: Vector3f,
    pub t_max: Float,
    pub time: Float,
    pub material: T,
}

impl<T: ConstZero> Default for Ray<T> {
    fn default() -> Self {
        Self::ZERO
    }
}

impl<T: ConstZero> ConstZero for Ray<T> {
    /// Ray from (0,0,0) towards (0,0,0)
    const ZERO: Self = Self {
        origin: Point3f::ZERO,
        direction: Vector3f::ZERO,
        t_max: Float::INFINITY,
        time: 0.0,
        material: T::ZERO,
    };
}

impl<T: Default> Ray<T> {
    /// Create a new ray with a given origin and direction
    pub fn new(origin: Point3f, direction: Vector3f) -> Self {
        Self {
            origin,
            direction,
            t_max: Float::INFINITY,
            time: 0.0,
            material: Default::default(),
        }
    }
}

impl<T> Ray<T> {
    /// Create a new ray with a given origin and direction
    pub fn new_with(origin: Point3f, direction: Vector3f, data: T) -> Self {
        Self {
            origin,
            direction,
            t_max: Float::INFINITY,
            time: 0.0,
            material: data,
        }
    }

    /// Calculate a point at t distance along the ray
    pub fn at(&self, t: Float) -> Point3f {
        self.origin + self.direction * t
    }
}

/// Information about two rays
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RayDifferentials {
    pub rx_origin: Point3f,
    pub ry_origin: Point3f,
    pub rx_direction: Vector3f,
    pub ry_direction: Vector3f,
}

/// Ray with differential information
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RayDifferential<T = ()> {
    pub main: Ray<T>,
    pub differentials: Option<RayDifferentials>,
}

impl<T: ConstZero> RayDifferential<T> {
    /// Ray::Zero, differentials = None
    pub const ZERO: Self = Self {
        main: Ray::ZERO,
        differentials: None,
    };
}

impl<T: Default> RayDifferential<T> {
    /// Create a new differential ray with no recorded differentials
    pub fn new(origin: Point3f, direction: Vector3f) -> Self {
        Self {
            main: Ray::new(origin, direction),
            differentials: None,
        }
    }
}

impl<T> RayDifferential<T> {
    /// Add differential storage to a given ray
    pub fn from_ray(ray: Ray<T>) -> Self {
        Self {
            main: ray,
            differentials: None,
        }
    }
}

impl<T: Copy> RayDifferential<T> {
    /// Scale the differentials by a given spacing
    pub fn scale(&self, scale: Float) -> Self {
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

impl<T> Deref for RayDifferential<T> {
    type Target = Ray<T>;

    fn deref(&self) -> &Self::Target {
        &self.main
    }
}
