use std::ops::{Index, Range};

use crate::{float::gamma, number::Number, ConstZero, Float, Point3, Ray, Vector3};

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

    /// Return an iterator over all the corners of the bounds
    pub fn iter_corners(&self) -> impl Iterator<Item = Point3<T>> {
        CornerIter {
            bound: *self,
            iter: 0..8,
        }
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
            o.x /= self.max.x - self.min.x
        }
        if self.max.y > self.min.y {
            o.y /= self.max.y - self.min.y
        }
        if self.max.z > self.min.z {
            o.z /= self.max.z - self.min.z
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

    /// Get the smallest squared distance between the bounds and a point.
    /// Returns 0 if the point is inside the bounds
    pub fn distance_squared(&self, point: Point3<T>) -> T {
        let d = (self.min - point).max(point - self.max).max(Vector3::ZERO);

        d.dot(d)
    }

    /// Get the smallest distance between the bounds and a point.
    /// Returns 0 if the point is inside the bounds
    pub fn distance(&self, point: Point3<T>) -> T {
        self.distance_squared(point).sqrt()
    }

    /// Calculate the intersection point between the bounds and a given ray
    pub fn intersect_p<U>(&self, ray: Ray<U, T>) -> Option<(T, T)> {
        let mut t0 = T::ZERO;
        let mut t1 = ray.t_max;

        for i in 0..3 {
            let inv_ray_dir = T::ONE / ray.direction[i];
            let mut t_near = (self.min[i] - ray.origin[i]) / inv_ray_dir;
            let mut t_far = (self.max[i] - ray.origin[i]) / inv_ray_dir;

            if t_near > t_far {
                std::mem::swap(&mut t_near, &mut t_far);
            }

            t_far *= T::ONE + T::TWO * gamma(3);
            t0 = t0.max(t_near);
            t1 = t1.min(t_far);
            if t0 > t1 {
                return None;
            }
        }

        Some((t0, t1))
    }

    /// Calculate the intersection point between the bounds and a given ray
    /// given the inverse of the ray directions pre-computed
    pub fn intersect_inv<U>(&self, ray: Ray<U, T>, inv: Vector3<T>, is_neg: [bool; 3]) -> bool {
        let mut t_min = (self[is_neg[0] as usize].x - ray.origin.x) * inv.x;
        let mut t_max = (self[1 - is_neg[0] as usize].x - ray.origin.x) * inv.x;
        let ty_min = (self[is_neg[1] as usize].y - ray.origin.y) * inv.y;
        let mut ty_max = (self[1 - is_neg[1] as usize].y - ray.origin.y) * inv.y;

        t_max *= T::ONE + T::TWO * gamma(3);
        ty_max *= T::ONE + T::TWO * gamma(3);

        if t_min > ty_max || ty_min > t_max {
            return false;
        }
        t_min = t_min.max(ty_min);
        t_max = t_max.min(ty_max);

        let tz_min = (self[is_neg[2] as usize].z - ray.origin.z) * inv.z;
        let mut tz_max = (self[1 - is_neg[2] as usize].z - ray.origin.z) * inv.z;
        tz_max *= T::ONE + T::TWO * gamma(3);

        if t_min > tz_max || tz_min > t_max {
            return false;
        }
        t_min = t_min.max(tz_min);
        t_max = t_max.min(tz_max);

        t_min < ray.t_max && t_max > T::ZERO
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

struct CornerIter<T: Number> {
    bound: Bounds3<T>,
    iter: Range<usize>,
}

impl<T: Number> Iterator for CornerIter<T> {
    type Item = Point3<T>;

    fn next(&mut self) -> Option<Self::Item> {
        Some(self.bound.corner(self.iter.next()?))
    }
}
