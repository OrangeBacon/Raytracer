use crate::{Float, Normal3f, Point2f, Point3f, Vector3f};

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Interaction {
    pub point: Point3f,
    pub time: Float,
    pub error: Vector3f,
    pub dir: Vector3f,
    pub normal: Normal3f,
    pub medium_interface: (),
}

impl Interaction {
    pub fn new(
        point: Point3f,
        normal: Normal3f,
        error: Vector3f,
        dir: Vector3f,
        time: Float,
        medium_interface: (),
    ) -> Self {
        Self {
            point,
            time,
            error,
            dir,
            normal,
            medium_interface,
        }
    }

    pub fn is_surface(&self) -> bool {
        self.normal != Default::default()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct SurfaceInteraction {
    pub interaction: Interaction,
    pub uv: Point2f,
    pub dpdu: Vector3f,
    pub dpdv: Vector3f,
    pub dndu: Normal3f,
    pub dndv: Normal3f,
    pub shape: (),
    pub shading: Shading,
}

impl SurfaceInteraction {
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
        shape: (),
    ) -> Self {
        let interaction = Interaction::new(
            point,
            Normal3f::from_vector(dpdu.cross(dpdv).normalise()),
            error,
            dir,
            time,
            (),
        );

        Self {
            interaction,
            uv,
            dpdu,
            dpdv,
            dndu,
            dndv,
            shape,
            shading: Shading {
                normal: interaction.normal,
                dpdu,
                dpdv,
                dndu,
                dndv,
            },
        }
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
