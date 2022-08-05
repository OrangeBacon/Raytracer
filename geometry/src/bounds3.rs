use std::ops::Index;

use crate::{number::Number, Float, Point3, Vector3};

/// 3D Axis aligned bounding box
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Bounds3<T: Number> {
    pub min: Point3<T>,
    pub max: Point3<T>,
}

pub type Bounds3f = Bounds3<Float>;
pub type Bounds3i = Bounds3<i32>;

impl<T: Number> Default for Bounds3<T> {
    fn default() -> Self {
        Self::ZERO
    }
}

impl<T: Number> Bounds3<T> {
    /// Bounds that take up no area
    pub const ZERO: Self = Self {
        min: Point3::MAX,
        max: Point3::MIN,
    };

    /// Create a bounds that contains both the given points
    pub fn new(p1: Point3<T>, p2: Point3<T>) -> Self {
        Self {
            min: p1.min(p2),
            max: p1.max(p2),
        }
    }

    /// Create a bounds that contains only the given point
    pub fn at_point(point: Point3<T>) -> Self {
        Self {
            min: point,
            max: point,
        }
    }

    /// Get the coordinates of one corner of the bounds
    pub fn corner(&self, corner: usize) -> Point3<T> {
        Point3::new(
            self[corner & 1].x,
            self[(corner & 2) >> 1].y,
            self[(corner & 4) >> 2].z,
        )
    }

    /// create bounding box containing self and the provided point
    pub fn union_point(&self, point: Point3<T>) -> Self {
        Self {
            min: self.min.min(point),
            max: self.max.max(point),
        }
    }

    /// create bounding box containing self and the provided bounding box
    pub fn union_box(&self, other: Self) -> Self {
        Self {
            min: self.min.min(other.min),
            max: self.max.max(other.max),
        }
    }

    /// create bounding box for the area that self and the other both overlap
    pub fn intersect(&self, other: Self) -> Self {
        Self {
            min: self.min.max(other.min),
            max: self.max.min(other.max),
        }
    }

    /// Does bound intersect with self
    pub fn overlap(&self, other: Self) -> bool {
        self.max >= other.min && self.min <= other.max
    }

    /// Is the given point inside self
    pub fn inside(&self, point: Point3<T>) -> bool {
        point >= self.min && point <= self.max
    }

    /// Is the given point inside self, excluding all boundary walls
    pub fn inside_exclusive(&self, point: Point3<T>) -> bool {
        point > self.min && point < self.max
    }

    /// Expand the bounds by a constant factor in all directions
    pub fn expand<U: Number>(&self, delta: U) -> Self {
        Self {
            min: self.min - Vector3::splat(delta).cast(),
            max: self.max + Vector3::splat(delta).cast(),
        }
    }

    /// Get a vector from the minimum to the maximum points of the bounds
    pub fn diagonal(&self) -> Vector3<T> {
        self.max - self.min
    }

    /// Surface area of the box created by the bounding box
    pub fn surface_area(&self) -> T {
        let diagonal = self.diagonal();
        T::TWO * diagonal.dot(diagonal)
    }

    /// The volume enclosed by the bounding box
    pub fn volume(&self) -> T {
        let diagonal = self.diagonal();
        diagonal.x * diagonal.y * diagonal.z
    }

    /// Index of the longest axis (x = 0; y = 1; z = 2)
    pub fn maximum_extent(&self) -> usize {
        let diagonal = self.diagonal();
        if diagonal.x > diagonal.y && diagonal.x > diagonal.z {
            0
        } else if diagonal.y > diagonal.z {
            1
        } else {
            2
        }
    }

    /// Linear interpolation between the minimum and maximum coordinates
    pub fn lerp(&self, t: T) -> Point3<T> {
        self.min.lerp(self.max, t)
    }

    /// Vector relative to the coordinates of the box given a location in the box
    /// where (0, 0, 0) is the minimum bound and (1, 1, 1) is the maximum
    pub fn offset(&self, point: Point3<T>) -> Vector3<T> {
        let mut o = point - self.min;
        if self.max.x > self.min.x {
            o.x = o.x / (self.max.x - self.min.x)
        }
        if self.max.y > self.min.y {
            o.y = o.y / (self.max.y - self.min.y)
        }
        if self.max.z > self.min.z {
            o.z = o.z / (self.max.z - self.min.z)
        }

        o
    }

    /// Get a sphere that contains this bounding box
    pub fn bounding_sphere(&self) -> (Point3<T>, T) {
        let centre = (self.min + self.max) / T::TWO;
        let radius = if self.inside(centre) {
            centre.distance(self.max)
        } else {
            T::ZERO
        };

        (centre, radius)
    }
}

impl<T: Number> Index<usize> for Bounds3<T> {
    type Output = Point3<T>;

    fn index(&self, index: usize) -> &Self::Output {
        assert!(index == 0 || index == 1);
        if index == 0 {
            &self.min
        } else {
            &self.max
        }
    }
}
