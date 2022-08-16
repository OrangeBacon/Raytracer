use std::ops::Mul;

use crate::{Float, Normal3f, Point2f, Point3f, Transform, Vector3f};

/// interaction at a point on a surface
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Interaction<T = ()> {
    pub point: Point3f,
    pub time: Float,
    pub error: Vector3f,
    pub dir: Vector3f,
    pub normal: Option<Normal3f>,
    pub medium_interface: T,
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
    /// Create a new record of an interaction at a point on a surface with a
    /// given medium interface
    pub fn new_with(
        point: Point3f,
        normal: Normal3f,
        error: Vector3f,
        dir: Vector3f,
        time: Float,
        medium_interface: T,
    ) -> Self {
        Self {
            medium_interface,
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
}

pub trait SurfaceInteractable {
    fn reverses_orientation(&self) -> bool;
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct SurfaceInteraction<S: SurfaceInteractable, T = ()> {
    pub interaction: Interaction<T>,
    pub uv: Point2f,
    pub dpdu: Vector3f,
    pub dpdv: Vector3f,
    pub dndu: Normal3f,
    pub dndv: Normal3f,
    pub shape: Option<S>,
    pub shading: Shading,
}

impl<S: SurfaceInteractable, T: Default> SurfaceInteraction<S, T> {
    /// Create a new surface interaction
    pub fn new(
        point: Point3f,
        error: Vector3f,
        uv: Point2f,
        dir: Vector3f,
        dpdu: Vector3f,
        dpdv: Vector3f,
        dndu: Normal3f,
        dndv: Normal3f,
        time: Float,
    ) -> Self {
        Self::new_shape_medium(
            point,
            error,
            uv,
            dir,
            dpdu,
            dpdv,
            dndu,
            dndv,
            time,
            Default::default(),
            None,
        )
    }
}

impl<S: SurfaceInteractable, T: Default> SurfaceInteraction<S, T> {
    /// Create a new surface interaction
    pub fn new_shape(
        point: Point3f,
        error: Vector3f,
        uv: Point2f,
        dir: Vector3f,
        dpdu: Vector3f,
        dpdv: Vector3f,
        dndu: Normal3f,
        dndv: Normal3f,
        time: Float,
        shape: S,
    ) -> Self {
        Self::new_shape_medium(
            point,
            error,
            uv,
            dir,
            dpdu,
            dpdv,
            dndu,
            dndv,
            time,
            Default::default(),
            Some(shape),
        )
    }
}

impl<S: SurfaceInteractable, T> SurfaceInteraction<S, T> {
    /// Create a new surface interaction
    pub fn new_medium(
        point: Point3f,
        error: Vector3f,
        uv: Point2f,
        dir: Vector3f,
        dpdu: Vector3f,
        dpdv: Vector3f,
        dndu: Normal3f,
        dndv: Normal3f,
        time: Float,
        medium_interface: T,
    ) -> Self {
        Self::new_shape_medium(
            point,
            error,
            uv,
            dir,
            dpdu,
            dpdv,
            dndu,
            dndv,
            time,
            medium_interface,
            None,
        )
    }

    /// Create a new surface interaction
    pub fn new_shape_medium(
        point: Point3f,
        error: Vector3f,
        uv: Point2f,
        dir: Vector3f,
        dpdu: Vector3f,
        dpdv: Vector3f,
        dndu: Normal3f,
        dndv: Normal3f,
        time: Float,
        medium_interface: T,
        shape: Option<S>,
    ) -> Self {
        let normal = Normal3f::from_vector(dpdu.cross(dpdv).normalise());
        let interaction = Interaction::new_with(point, normal, error, dir, time, medium_interface);

        let mut this = Self {
            interaction,
            uv,
            dpdu,
            dpdv,
            dndu,
            dndv,
            shape,
            shading: Shading {
                normal,
                dpdu,
                dpdv,
                dndu,
                dndv,
            },
        };

        if let Some(shape) = &this.shape {
            if shape.reverses_orientation() {
                this.interaction.normal = Some(normal * -1.0);
                this.shading.normal *= -1.0;
            }
        }

        this
    }

    /// Change the shading parameters associated with an interaction
    pub fn set_shading_geometry(
        &mut self,
        dpdus: Vector3f,
        dpdvs: Vector3f,
        dndus: Normal3f,
        dndvs: Normal3f,
        orientation_is_authoritative: bool,
    ) {
        self.shading.normal = Normal3f::from_vector(dpdus.cross(dpdvs)).normalise();

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

        self.shading.dpdu = dpdus;
        self.shading.dpdv = dpdvs;
        self.shading.dndu = dndus;
        self.shading.dndv = dndvs;
    }
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Shading {
    pub normal: Normal3f,
    pub dpdu: Vector3f,
    pub dpdv: Vector3f,
    pub dndu: Normal3f,
    pub dndv: Normal3f,
}

impl<S: SurfaceInteractable, T> Mul<SurfaceInteraction<S, T>> for Transform {
    type Output = SurfaceInteraction<S, T>;

    fn mul(self, _rhs: SurfaceInteraction<S, T>) -> Self::Output {
        todo!()
    }
}
