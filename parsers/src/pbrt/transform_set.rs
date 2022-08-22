use std::ops::{Index, IndexMut};

use geometry::{Float, Matrix4x4, Number, Point3f, Transform, Vector3, Vector3f, Point3};

/// Set of transformations for multiple points in time
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct TransformSet<T: Number> {
    transforms: [Transform<T>; 2],
    active_transforms: usize,
}

impl<T: Number> TransformSet<T> {
    /// The number of transformations included
    const MAX_TRANSFORMS: usize = 2;

    /// The starting transform
    const START_TRANSFORM: usize = 1 << 0;

    /// The ending transform
    const END_TRANSFORM: usize = 1 << 1;

    /// All possible transforms
    const ALL_TRANSFORM: usize = (1 << Self::MAX_TRANSFORMS) - 1;

    /// Create a new transform set with all transforms being the identity
    pub fn new() -> Self {
        Self {
            transforms: [Transform::IDENTITY; 2],
            active_transforms: Self::ALL_TRANSFORM,
        }
    }

    /// Create a new transform set that is the inverse of the current one
    pub fn inverse(&self) -> Self {
        Self {
            transforms: self.transforms.map(|t| t.inverse()),
            active_transforms: self.active_transforms,
        }
    }

    /// run a function over all enabled transforms
    fn for_active(&mut self, f: impl Fn(&mut Transform<T>)) {
        for (idx, t) in self.transforms.iter_mut().enumerate() {
            if (self.active_transforms & (1 << idx)) != 0 {
                f(t);
            }
        }
    }

    /// Set all active transforms to the identity matrix
    pub fn identity(&mut self) {
        self.for_active(|t| *t = Transform::IDENTITY);
    }

    /// Translate active transforms by (dx, dy, dz)
    pub fn translate(&mut self, dx: T, dy: T, dz: T) {
        self.for_active(|t| *t *= Transform::translation(Vector3::new(dx, dy, dz)))
    }

    /// Rotate active transforms by angle degrees around axis (ax, ay, az)
    pub fn rotate(&mut self, angle: T, ax: T, ay: T, az: T) {
        self.for_active(|t| *t *= Transform::rotate(Vector3::new(ax, ay, az), angle))
    }

    /// Scale active transforms by (sx, sy, sz)
    pub fn scale(&mut self, sx: T, sy: T, sz: T) {
        self.for_active(|t| *t *= Transform::scale(Vector3::new(sx, sy, sz)))
    }

    /// Set the active transforms to be a projection matrix, looking from pos to look.
    pub fn look_at(&mut self, pos: Point3<T>, look: Point3<T>, up: Vector3<T>) {
        self.for_active(|t| {
            *t = Transform::look_at(pos, look, up).unwrap_or_else(|| {
                panic!("Unable to create look at transform with pos: {pos:?}, look: {look:?}, up: {up:?}")
            })
        })
    }

    /// Post multiply the active transforms by the given matrix
    pub fn concat(&mut self, transform: [T; 16]) {
        let mat = Matrix4x4::from_array(&transform);
        let transform = Transform::from_mat(&mat)
            .unwrap_or_else(|| panic!("Unable to create multiplicative transform {mat:?}"));
        self.for_active(|t| *t *= transform);
    }

    /// Set all active transforms to the given matrix
    pub fn set_transform(&mut self, transform: [T; 16]) {
        let mat = Matrix4x4::from_array(&transform);
        let transform = Transform::from_mat(&mat)
            .unwrap_or_else(|| panic!("Unable to create multiplicative transform {mat:?}"));
        self.for_active(|t| *t = transform);
    }

    /// is this a set of animated transformations
    pub fn is_animated(&self) -> bool {
        for idx in 0..self.transforms.len() - 1 {
            if self.transforms[idx] != self.transforms[idx + 1] {
                return true;
            }
        }
        false
    }

    /// Set all transforms to be active
    pub fn all(&mut self) {
        self.active_transforms = Self::ALL_TRANSFORM;
    }

    /// Set only the start transform to be active
    pub fn start(&mut self) {
        self.active_transforms = Self::START_TRANSFORM;
    }

    /// Set only the end transform to be active
    pub fn end(&mut self) {
        self.active_transforms = Self::END_TRANSFORM;
    }
}

impl<T: Number> Index<usize> for TransformSet<T> {
    type Output = Transform<T>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.transforms[index]
    }
}

impl<T: Number> IndexMut<usize> for TransformSet<T> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.transforms[index]
    }
}
