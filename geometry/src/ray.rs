use std::ops::Deref;

use crate::{ConstZero, Number, Point3, Vector3};

/// Ray with an origin and a direction from the origin
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Ray<T: Number> {
    pub origin: Point3<T>,
    pub direction: Vector3<T>,
    pub t_max: T,
    pub time: T,
}

impl<T: Number> Default for Ray<T> {
    fn default() -> Self {
        Self::ZERO
    }
}

impl<T: Number> ConstZero for Ray<T> {
    /// Ray from (0,0,0) towards (0,0,0)
    const ZERO: Self = Self {
        origin: Point3::ZERO,
        direction: Vector3::ZERO,
        t_max: T::INFINITY,
        time: T::ZERO,
    };
}

impl<T: Number> Ray<T> {
    /// Create a new ray with a given origin and direction
    pub fn new(origin: Point3<T>, direction: Vector3<T>) -> Self {
        Self {
            origin,
            direction,
            t_max: T::INFINITY,
            time: T::ZERO,
        }
    }

    /// Create a new ray with a given origin and direction
    pub fn new_with(origin: Point3<T>, direction: Vector3<T>) -> Self {
        Self {
            origin,
            direction,
            t_max: T::INFINITY,
            time: T::ZERO,
        }
    }

    /// Calculate a point at t distance along the ray
    pub fn at(&self, t: T) -> Point3<T> {
        self.origin + self.direction * t
    }
}

/// Information about two rays
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default, Hash)]
pub struct RayDifferentials<T: Number> {
    pub rx_origin: Point3<T>,
    pub ry_origin: Point3<T>,
    pub rx_direction: Vector3<T>,
    pub ry_direction: Vector3<T>,
}

/// Ray with differential information
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct RayDifferential<T: Number> {
    pub main: Ray<T>,
    pub differentials: Option<RayDifferentials<T>>,
}

impl<T: Number> RayDifferential<T> {
    /// Ray::Zero, differentials = None
    pub const ZERO: Self = Self {
        main: Ray::ZERO,
        differentials: None,
    };

    /// Create a new differential ray with no recorded differentials
    pub fn new(origin: Point3<T>, direction: Vector3<T>) -> Self {
        Self {
            main: Ray::new(origin, direction),
            differentials: None,
        }
    }

    /// Add differential storage to a given ray
    pub fn from_ray(ray: Ray<T>) -> Self {
        Self {
            main: ray,
            differentials: None,
        }
    }

    /// Scale the differentials by a given spacing
    pub fn scale(&self, scale: T) -> Self {
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

impl<T: Number> Deref for RayDifferential<T> {
    type Target = Ray<T>;

    fn deref(&self) -> &Self::Target {
        &self.main
    }
}
