use std::ops::{Mul, MulAssign};

use crate::{transform::Applicable, Float, Normal3f, Point2f, Point3f, Transform, Vector3f};

/// interaction at a point on a surface
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Interaction<T = ()> {
    /// location of the interaction
    pub point: Point3f,

    /// Time that the interaction occured
    pub time: Float,

    /// Error bound on calculations with this interaction
    pub error: Vector3f,

    /// the negative ray direction at the interaction
    pub dir: Vector3f,

    /// The direction of the normal at the point this intersection occurred
    pub normal: Option<Normal3f>,

    /// Additional data at this interaction point (supposed to be a medium interface
    /// but cargo does not like recursive dependencies between crates)
    pub medium_interface: T,
}

/// Parametric partial derivatives of a point
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct PartialDerivatives {
    // Partial derivatives of a point in the tangent plane to the interaction
    pub dpdu: Vector3f,
    pub dpdv: Vector3f,

    // Partial derivatives of the normal vector of the interaction
    pub dndu: Normal3f,
    pub dndv: Normal3f,
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
/// Different partial derivatives and normals for particular points in the geometry of a shape.
/// Used in bump mapping and other different surfaces.
pub struct Shading {
    pub normal: Normal3f,
    pub derivatives: PartialDerivatives,
}

/// Small portion of the shape interface representing the parts that are used
/// in a surface interaction.  Only exists to fix recursive cargo dependencies
pub trait SurfaceInteractable {
    /// Does the shape's transform reverse its orientation.
    /// Equal to shape.reverse_orientation ^ shape.transform.swaps_handedness()
    fn reverses_orientation(&self) -> bool;
}

/// Interaction between a ray and a generic point on a surface
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct SurfaceInteraction<S: SurfaceInteractable, T = ()> {
    pub interaction: Interaction<T>,
    pub uv: Point2f,
    pub derivatives: PartialDerivatives,
    pub shape: Option<S>,
    pub shading: Shading,
}

impl<T: Default> Interaction<T> {
    /// Create a new record of an interaction at a point on a surface
    pub fn new(
        point: Point3f,
        normal: Normal3f,
        error: Vector3f,
        dir: Vector3f,
        time: Float,
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

impl<T> Interaction<T> {
    /// Add a medium interface to this interaction
    pub fn with_medium<U>(self, medium: U) -> Interaction<U> {
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

impl SurfaceInteraction<(), ()> {
    /// Create a new surface interaction
    pub fn new(
        point: Point3f,
        error: Vector3f,
        uv: Point2f,
        dir: Vector3f,
        derivatives: PartialDerivatives,
        time: Float,
    ) -> Self {
        let normal = Normal3f::from_vector(derivatives.dpdu.cross(derivatives.dpdv).normalise());
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

impl<S: SurfaceInteractable, T> SurfaceInteraction<S, T> {
    /// Associate a medium interface with this surface interaction
    pub fn with_medium<U>(self, medium: U) -> SurfaceInteraction<S, U> {
        SurfaceInteraction {
            interaction: self.interaction.with_medium(medium),
            uv: self.uv,
            derivatives: self.derivatives,
            shape: self.shape,
            shading: self.shading,
        }
    }

    /// Associate a shape with this surface interaction
    pub fn with_shape<S2: SurfaceInteractable>(self, shape: S2) -> SurfaceInteraction<S2, T> {
        let mut this = SurfaceInteraction {
            interaction: self.interaction,
            uv: self.uv,
            derivatives: self.derivatives,
            shape: Some(shape),
            shading: self.shading,
        };

        if let Some(shape) = &this.shape {
            if shape.reverses_orientation() {
                this.interaction.normal = Some(this.shading.normal * -1.0);
                this.shading.normal *= -1.0;
            }
        }

        this
    }

    fn remove_params(&self) -> SurfaceInteraction<(), ()> {
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
        derivatives: PartialDerivatives,
        orientation_is_authoritative: bool,
    ) {
        self.shading.normal =
            Normal3f::from_vector(derivatives.dpdu.cross(derivatives.dpdv)).normalise();

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

impl<S: SurfaceInteractable, T> Mul<Transform> for SurfaceInteraction<S, T> {
    type Output = SurfaceInteraction<S, T>;

    fn mul(self, rhs: Transform) -> Self::Output {
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

impl<S: SurfaceInteractable, T> MulAssign<Transform> for SurfaceInteraction<S, T> {
    fn mul_assign(&mut self, rhs: Transform) {
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

impl<S: SurfaceInteractable, T> Applicable<SurfaceInteraction<S, T>> for Transform {
    fn apply(&self, other: SurfaceInteraction<S, T>) -> SurfaceInteraction<S, T> {
        other * *self
    }
}
