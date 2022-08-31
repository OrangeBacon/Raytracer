use std::{
    ops::{Deref, DerefMut},
    sync::Arc,
};

use geometry::{Bounds3, Number, Ray, SurfaceInteractable, SurfaceInteraction, Transform};

use crate::{triangle_mesh::TriangleMesh, Shape, ShapeData};

pub struct Triangle<F: Number, T> {
    data: ShapeData<F>,
    mesh: Arc<TriangleMesh<F, T>>,
    idx: usize,
}

impl<F: Number, T> Triangle<F, T> {
    /// Create a new triangle from the data stored in a triangle mesh
    pub fn new(
        object_to_world: Transform<F>,
        world_to_object: Transform<F>,
        reverse_orientation: bool,
        mesh: Arc<TriangleMesh<F, T>>,
        tri_number: usize,
    ) -> Self {
        Self {
            data: ShapeData::new(object_to_world, world_to_object, reverse_orientation),
            mesh,
            idx: tri_number * 3,
        }
    }

    /// Create all the triangles represented by a triangle mesh
    pub fn from_mesh(
        mesh: Arc<TriangleMesh<F, T>>,
        world_to_object: Transform<F>,
        reverse_orientation: bool,
    ) -> Vec<Arc<dyn Shape<F>>>
    where
        F: 'static,
        T: 'static,
    {
        let tri_count = mesh.positions().len() / 3;

        let mut tris: Vec<Arc<dyn Shape<F>>> = Vec::with_capacity(tri_count);
        for i in 0..tris.capacity() {
            let tri = Triangle::new(
                *mesh.object_to_world(),
                world_to_object,
                reverse_orientation,
                Arc::clone(&mesh),
                i,
            );
            tris.push(Arc::new(tri));
        }

        tris
    }
}

impl<F: Number, T> Shape<F> for Triangle<F, T> {
    fn object_bound(&self) -> Bounds3<F> {
        let p0 = self.mesh.positions()[self.idx];
        let p1 = self.mesh.positions()[self.idx + 1];
        let p2 = self.mesh.positions()[self.idx + 2];

        Bounds3::new(p0 * self.world_to_object, p1 * self.world_to_object)
            .union_point(p2 * self.world_to_object)
    }

    fn world_bound(&self) -> Bounds3<F> {
        let p0 = self.mesh.positions()[self.idx];
        let p1 = self.mesh.positions()[self.idx + 1];
        let p2 = self.mesh.positions()[self.idx + 2];

        Bounds3::new(p0, p1).union_point(p2)
    }

    fn intersect(
        &self,
        ray: Ray<(), F>,
        _test_alpha: bool,
    ) -> Option<(F, SurfaceInteraction<&dyn SurfaceInteractable, (), F>)> {
        todo!()
    }

    fn area(&self) -> F {
        let p0 = self.mesh.positions()[self.idx];
        let p1 = self.mesh.positions()[self.idx + 1];
        let p2 = self.mesh.positions()[self.idx + 2];

        (F::ONE / F::TWO) * (p1 - p0).cross(p2 - p0).length()
    }
}

impl<F: Number, T> SurfaceInteractable for Triangle<F, T> {
    fn reverses_orientation(&self) -> bool {
        self.reverse_orientation ^ self.transform_swaps_handedness
    }
}

impl<F: Number, T> Deref for Triangle<F, T> {
    type Target = ShapeData<F>;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl<F: Number, T> DerefMut for Triangle<F, T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}
