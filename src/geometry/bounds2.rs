use std::ops::Index;

use crate::geometry::{number::Number, Float, Point2, Vector2};

/// 3D Axis aligned bounding box
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Bounds2<T: Number> {
    pub min: Point2<T>,
    pub max: Point2<T>,
}

pub type Bounds2f = Bounds2<Float>;
pub type Bounds2i = Bounds2<i32>;

impl<T: Number> Default for Bounds2<T> {
    fn default() -> Self {
        Self::ZERO
    }
}

impl<T: Number> Bounds2<T> {
    /// Bounds that take up no area
    pub const ZERO: Self = Self {
        min: Point2::MAX,
        max: Point2::MIN,
    };

    /// Create a bounds that contains both the given points
    pub fn new(p1: Point2<T>, p2: Point2<T>) -> Self {
        Self {
            min: p1.min(p2),
            max: p1.max(p2),
        }
    }

    /// Create a bounds that contains only the given point
    pub fn at_point(point: Point2<T>) -> Self {
        Self {
            min: point,
            max: point,
        }
    }

    /// Get the coordinates of one corner of the bounds
    pub fn corner(&self, corner: usize) -> Point2<T> {
        Point2::new(self[corner & 1].x, self[(corner & 2) >> 1].y)
    }

    /// create bounding box containing self and the provided point
    pub fn union_point(&self, point: Point2<T>) -> Self {
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
    pub fn inside(&self, point: Point2<T>) -> bool {
        point >= self.min && point <= self.max
    }

    /// Is the given point inside self, excluding all boundary walls
    pub fn inside_exclusive(&self, point: Point2<T>) -> bool {
        point > self.min && point < self.max
    }

    /// Expand the bounds by a constant factor in all directions
    pub fn expand<U: Number>(&self, delta: U) -> Self {
        Self {
            min: self.min - Vector2::splat(delta).cast(),
            max: self.max + Vector2::splat(delta).cast(),
        }
    }

    /// Get a vector from the minimum to the maximum points of the bounds
    pub fn diagonal(&self) -> Vector2<T> {
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
        diagonal.x * diagonal.y
    }

    /// Index of the longest axis (x = 0; y = 1)
    pub fn maximum_extent(&self) -> usize {
        let diagonal = self.diagonal();
        if diagonal.x > diagonal.y {
            0
        } else {
            1
        }
    }

    /// Linear interpolation between the minimum and maximum coordinates
    pub fn lerp(&self, t: T) -> Point2<T> {
        self.min.lerp(self.max, t)
    }

    /// Vector relative to the coordinates of the box given a location in the box
    /// where (0, 0, 0) is the minimum bound and (1, 1, 1) is the maximum
    pub fn offset(&self, point: Point2<T>) -> Vector2<T> {
        let mut o = point - self.min;
        if self.max.x > self.min.x {
            o.x = o.x / (self.max.x - self.min.x)
        }
        if self.max.y > self.min.y {
            o.y = o.y / (self.max.y - self.min.y)
        }

        o
    }

    /// Get a circle that contains this bounding box
    pub fn bounding_circle(&self) -> (Point2<T>, T) {
        let centre = (self.min + self.max) / T::TWO;
        let radius = if self.inside(centre) {
            centre.distance(self.max)
        } else {
            T::ZERO
        };

        (centre, radius)
    }
}

impl<T: Number> Index<usize> for Bounds2<T> {
    type Output = Point2<T>;

    fn index(&self, index: usize) -> &Self::Output {
        assert!(index == 0 || index == 1);
        if index == 0 {
            &self.min
        } else {
            &self.max
        }
    }
}
