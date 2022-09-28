use std::{
    fmt::Debug,
    ops::{Mul, MulAssign},
};

use geometry::{
    offset_ray_origin, Applicable, Normal3, Number, Point2, Point3, Ray, Transform, Vector3,
};

use crate::{primitive::Primitive, Shape};

/// interaction at a point on a surface
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default, Hash)]
pub struct Interaction<T: Number> {
    /// location of the interaction
    pub point: Point3<T>,

    /// Time that the interaction occurred
    pub time: T,

    /// Error bound on calculations with this interaction
    pub error: Vector3<T>,

    /// the negative ray direction at the interaction
    pub dir: Vector3<T>,

    /// The direction of the normal at the point this intersection occurred
    pub normal: Option<Normal3<T>>,
}

/// Parametric partial derivatives of a point
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default, Hash)]
pub struct PartialDerivatives<T: Number> {
    // Partial derivatives of a point in the tangent plane to the interaction
    pub dpdu: Vector3<T>,
    pub dpdv: Vector3<T>,

    // Partial derivatives of the normal vector of the interaction
    pub dndu: Normal3<T>,
    pub dndv: Normal3<T>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default, Hash)]
/// Different partial derivatives and normals for particular points in the geometry of a shape.
/// Used in bump mapping and other different surfaces.
pub struct Shading<T: Number> {
    pub normal: Normal3<T>,
    pub derivatives: PartialDerivatives<T>,
}

/// Interaction between a ray and a generic point on a surface
#[derive(Debug, Default)]
pub struct SurfaceInteraction<T: Number> {
    pub interaction: Interaction<T>,
    pub uv: Point2<T>,
    pub derivatives: PartialDerivatives<T>,
    pub shape: Option<Box<dyn Shape<T>>>,
    pub shading: Shading<T>,
    pub primitive: Option<Box<dyn Primitive<T>>>,
}

impl<T: Number> Interaction<T> {
    /// Create a new record of an interaction at a point on a surface
    pub fn new(
        point: Point3<T>,
        normal: Normal3<T>,
        error: Vector3<T>,
        dir: Vector3<T>,
        time: T,
    ) -> Self {
        Self {
            point,
            time,
            error,
            dir,
            normal: Some(normal),
        }
    }

    /// Is the interaction on a surface
    pub fn is_surface(&self) -> bool {
        self.normal.is_some()
    }

    /// Spawn a ray at the intersection point
    pub fn spawn_ray(&self, direction: Vector3<T>) -> Ray<T> {
        let origin = offset_ray_origin(
            self.point,
            self.error,
            self.normal.unwrap_or_default(),
            direction,
        );
        Ray {
            origin,
            direction,
            t_max: T::INFINITY,
            time: self.time,
        }
    }

    /// Spawn a ray ending just before a given point
    pub fn spawn_ray_to(&self, point: Point3<T>) -> Ray<T> {
        let origin = offset_ray_origin(
            self.point,
            self.error,
            self.normal.unwrap_or_default(),
            point - self.point,
        );
        let dir = point - origin;
        Ray {
            origin,
            direction: dir,
            t_max: T::ONE - T::cast(0.0001),
            time: self.time,
        }
    }

    /// Spawn a ray ending just before a given surface interaction
    pub fn spawn_ray_intersect(&self, isect: Interaction<T>) -> Ray<T> {
        let origin = offset_ray_origin(
            self.point,
            self.error,
            self.normal.unwrap_or_default(),
            isect.point - self.point,
        );
        let dir = isect.point - origin;
        Ray {
            origin,
            direction: dir,
            t_max: T::ONE - T::cast(0.0001),
            time: self.time,
        }
    }
}

impl<T: Number> SurfaceInteraction<T> {
    /// Create a new surface interaction
    pub fn new(
        point: Point3<T>,
        error: Vector3<T>,
        uv: Point2<T>,
        dir: Vector3<T>,
        derivatives: PartialDerivatives<T>,
        time: T,
    ) -> Self {
        let normal = Normal3::from_vector(derivatives.dpdu.cross(derivatives.dpdv).normalise());
        let interaction = Interaction::new(point, normal, error, dir, time);

        Self {
            interaction,
            uv,
            derivatives,
            shape: None,
            shading: Shading {
                normal,
                derivatives,
            },
            primitive: None,
        }
    }

    /// Associate a shape with this surface interaction
    pub fn with_shape<S2: Shape<T> + 'static>(self, shape: S2) -> Self {
        let mut this = Self {
            shape: Some(Box::new(shape)),
            ..self
        };

        if let Some(shape) = &this.shape {
            if shape.reverses_orientation() {
                this.interaction.normal = Some(this.shading.normal * -T::ONE);
                this.shading.normal *= -T::ONE;
            }
        }

        this
    }

    /// Associate a primitive with the interaction
    pub fn with_primitive<P2: Primitive<T> + 'static>(self, prim: P2) -> Self {
        SurfaceInteraction {
            primitive: Some(Box::new(prim)),
            ..self
        }
    }

    /// Change the shading parameters associated with an interaction
    pub fn set_shading_geometry(
        &mut self,
        derivatives: PartialDerivatives<T>,
        orientation_is_authoritative: bool,
    ) {
        self.shading.normal =
            Normal3::from_vector(derivatives.dpdu.cross(derivatives.dpdv)).normalise();

        if let Some(shape) = &self.shape {
            if shape.reverses_orientation() {
                self.shading.normal = -self.shading.normal;
            }
        }

        if let Some(n) = self.interaction.normal {
            if orientation_is_authoritative {
                self.interaction.normal = Some(n.face_forward(self.shading.normal.to_vector()))
            } else {
                self.shading.normal = self.shading.normal.face_forward(n.to_vector());
            }
        }

        self.shading.derivatives = derivatives;
    }
}

impl<T: Number> Mul<Transform<T>> for PartialDerivatives<T> {
    type Output = Self;

    fn mul(self, rhs: Transform<T>) -> Self::Output {
        Self {
            dpdu: self.dpdu * rhs,
            dpdv: self.dpdv * rhs,
            dndu: self.dndu * rhs,
            dndv: self.dndv * rhs,
        }
    }
}

impl<T: Number> Mul<Transform<T>> for Interaction<T> {
    type Output = Self;

    fn mul(self, rhs: Transform<T>) -> Self::Output {
        let (point, error) = rhs.apply_err((self.point, self.error));

        Interaction {
            point,
            time: self.time,
            error,
            dir: self.dir * rhs,
            normal: self.normal.map(|n| rhs.apply(n)),
        }
    }
}

impl<T: Number> Mul<Transform<T>> for SurfaceInteraction<T> {
    type Output = SurfaceInteraction<T>;

    fn mul(self, rhs: Transform<T>) -> Self::Output {
        SurfaceInteraction {
            interaction: self.interaction * rhs,
            derivatives: self.derivatives * rhs,
            shading: Shading {
                normal: (self.shading.normal * rhs).normalise(),
                derivatives: self.derivatives * rhs,
            },
            ..self
        }
    }
}

impl<T: Number> MulAssign<Transform<T>> for SurfaceInteraction<T> {
    fn mul_assign(&mut self, rhs: Transform<T>) {
        self.shading = Shading {
            normal: (self.shading.normal * rhs).normalise(),
            derivatives: self.derivatives * rhs,
        };
        self.derivatives = self.derivatives * rhs;
        self.interaction = self.interaction * rhs;
    }
}

impl<T: Number> Applicable<SurfaceInteraction<T>> for Transform<T> {
    fn apply(&self, other: SurfaceInteraction<T>) -> SurfaceInteraction<T> {
        other * *self
    }
}
