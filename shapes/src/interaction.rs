use std::{
    fmt::Debug,
    ops::{Mul, MulAssign},
};

use geometry::{
    offset_ray_origin, Applicable, Normal3, Number, Point2, Point3, Ray, Transform,
    Vector3,
};

/// interaction at a point on a surface
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default, Hash)]
pub struct Interaction<T, F: Number> {
    /// location of the interaction
    pub point: Point3<F>,

    /// Time that the interaction occurred
    pub time: F,

    /// Error bound on calculations with this interaction
    pub error: Vector3<F>,

    /// the negative ray direction at the interaction
    pub dir: Vector3<F>,

    /// The direction of the normal at the point this intersection occurred
    pub normal: Option<Normal3<F>>,

    /// Additional data at this interaction point (supposed to be a medium interface
    /// but cargo does not like recursive dependencies between crates)
    pub medium_interface: T,
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

/// Small portion of the shape interface representing the parts that are used
/// in a surface interaction.  Only exists to fix recursive cargo dependencies
pub trait SurfaceInteractable: Debug {
    /// Does the shape's transform reverse its orientation.
    /// Equal to shape.reverse_orientation ^ shape.transform.swaps_handedness()
    fn reverses_orientation(&self) -> bool;
}

/// Interaction between a ray and a generic point on a surface
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default, Hash)]
pub struct SurfaceInteraction<S: SurfaceInteractable, T, F: Number> {
    pub interaction: Interaction<T, F>,
    pub uv: Point2<F>,
    pub derivatives: PartialDerivatives<F>,
    pub shape: Option<S>,
    pub shading: Shading<F>,
}

impl<T: Default, F: Number> Interaction<T, F> {
    /// Create a new record of an interaction at a point on a surface
    pub fn new(
        point: Point3<F>,
        normal: Normal3<F>,
        error: Vector3<F>,
        dir: Vector3<F>,
        time: F,
    ) -> Self {
        Self {
            point,
            time,
            error,
            dir,
            normal: Some(normal),
            medium_interface: Default::default(),
        }
    }
}

impl<T, F: Number> Interaction<T, F> {
    /// Add a medium interface to this interaction
    pub fn with_medium<U>(self, medium: U) -> Interaction<U, F> {
        Interaction {
            medium_interface: medium,
            point: self.point,
            time: self.time,
            error: self.error,
            dir: self.dir,
            normal: self.normal,
        }
    }

    /// Is the interaction on a surface
    pub fn is_surface(&self) -> bool {
        self.normal.is_some()
    }

    /// Spawn a ray at the intersection point
    pub fn spawn_ray(&self, direction: Vector3<F>) -> Ray<(), F> {
        let origin = offset_ray_origin(
            self.point,
            self.error,
            self.normal.unwrap_or_default(),
            direction,
        );
        Ray {
            origin,
            direction,
            t_max: F::INFINITY,
            time: self.time,
            material: (),
        }
    }

    /// Spawn a ray ending just before a given point
    pub fn spawn_ray_to(&self, point: Point3<F>) -> Ray<(), F> {
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
            t_max: F::ONE - F::cast(0.0001),
            time: self.time,
            material: (),
        }
    }

    /// Spawn a ray ending just before a given surface interaction
    pub fn spawn_ray_intersect<U>(&self, isect: Interaction<U, F>) -> Ray<(), F> {
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
            t_max: F::ONE - F::cast(0.0001),
            time: self.time,
            material: (),
        }
    }
}

/// SurfaceInteractable for &dyn SurfaceInteractable
impl<T: SurfaceInteractable + ?Sized> SurfaceInteractable for &T {
    fn reverses_orientation(&self) -> bool {
        (*self).reverses_orientation()
    }
}

impl SurfaceInteractable for () {
    fn reverses_orientation(&self) -> bool {
        false
    }
}

impl<F: Number> SurfaceInteraction<(), (), F> {
    /// Create a new surface interaction
    pub fn new(
        point: Point3<F>,
        error: Vector3<F>,
        uv: Point2<F>,
        dir: Vector3<F>,
        derivatives: PartialDerivatives<F>,
        time: F,
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
        }
    }
}

impl<S: SurfaceInteractable, T, F: Number> SurfaceInteraction<S, T, F> {
    /// Associate a medium interface with this surface interaction
    pub fn with_medium<U>(self, medium: U) -> SurfaceInteraction<S, U, F> {
        SurfaceInteraction {
            interaction: self.interaction.with_medium(medium),
            uv: self.uv,
            derivatives: self.derivatives,
            shape: self.shape,
            shading: self.shading,
        }
    }

    /// Associate a shape with this surface interaction
    pub fn with_shape<S2: SurfaceInteractable>(self, shape: S2) -> SurfaceInteraction<S2, T, F> {
        let mut this = SurfaceInteraction {
            interaction: self.interaction,
            uv: self.uv,
            derivatives: self.derivatives,
            shape: Some(shape),
            shading: self.shading,
        };

        if let Some(shape) = &this.shape {
            if shape.reverses_orientation() {
                this.interaction.normal = Some(this.shading.normal * -F::ONE);
                this.shading.normal *= -F::ONE;
            }
        }

        this
    }

    fn remove_params(&self) -> SurfaceInteraction<(), (), F> {
        SurfaceInteraction {
            interaction: Interaction {
                point: self.interaction.point,
                time: self.interaction.time,
                error: self.interaction.error,
                dir: self.interaction.dir,
                normal: self.interaction.normal,
                medium_interface: (),
            },
            uv: self.uv,
            derivatives: self.derivatives,
            shape: Some(()),
            shading: self.shading,
        }
    }

    /// Change the shading parameters associated with an interaction
    pub fn set_shading_geometry(
        &mut self,
        derivatives: PartialDerivatives<F>,
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

impl<S: SurfaceInteractable, T, F: Number> Mul<Transform<F>> for SurfaceInteraction<S, T, F> {
    type Output = SurfaceInteraction<S, T, F>;

    fn mul(self, rhs: Transform<F>) -> Self::Output {
        let (point, error) = rhs.apply_err((self.interaction.point, self.interaction.error));

        SurfaceInteraction {
            interaction: Interaction {
                point,
                time: self.interaction.time,
                error,
                dir: self.interaction.dir * rhs,
                normal: self.interaction.normal.map(|n| rhs.apply(n)),
                medium_interface: self.interaction.medium_interface,
            },
            uv: self.uv,
            derivatives: PartialDerivatives {
                dpdu: self.derivatives.dpdu * rhs,
                dpdv: self.derivatives.dpdv * rhs,
                dndu: self.derivatives.dndu * rhs,
                dndv: self.derivatives.dndv * rhs,
            },
            shape: self.shape,
            shading: Shading {
                normal: (self.shading.normal * rhs).normalise(),
                derivatives: PartialDerivatives {
                    dpdu: self.shading.derivatives.dpdu * rhs,
                    dpdv: self.shading.derivatives.dpdv * rhs,
                    dndu: self.shading.derivatives.dndu * rhs,
                    dndv: self.shading.derivatives.dndv * rhs,
                },
            },
        }
    }
}

impl<S: SurfaceInteractable, T, F: Number> MulAssign<Transform<F>> for SurfaceInteraction<S, T, F> {
    fn mul_assign(&mut self, rhs: Transform<F>) {
        let res = self.remove_params() * rhs;
        self.derivatives = res.derivatives;
        self.interaction.dir = res.interaction.dir;
        self.interaction.error = res.interaction.error;
        self.interaction.normal = res.interaction.normal;
        self.interaction.point = res.interaction.point;
        self.interaction.time = res.interaction.time;
        self.shading = res.shading;
        self.uv = res.uv;
    }
}

impl<S: SurfaceInteractable, T, F: Number> Applicable<SurfaceInteraction<S, T, F>>
    for Transform<F>
{
    fn apply(&self, other: SurfaceInteraction<S, T, F>) -> SurfaceInteraction<S, T, F> {
        other * *self
    }
}
