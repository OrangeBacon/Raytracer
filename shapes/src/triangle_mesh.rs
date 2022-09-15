use std::sync::Arc;

use geometry::{Normal3, Number, Point2, Point3, Transform, Vector3};

/// Container for individual triangle vertex data
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Default, Hash)]
pub struct TriangleMesh<F: Number, T> {
    object_to_world: Transform<F>,

    position: Vec<Point3<F>>,
    index: Vec<usize>,
    tangent: Option<Vec<Vector3<F>>>,
    normal: Option<Vec<Normal3<F>>>,
    uv: Option<Vec<Point2<F>>>,
    alpha_mask: Option<Arc<T>>,
}

impl<F: Number> TriangleMesh<F, ()> {
    /// Create a new triangle mesh
    pub fn new(
        object_to_world: Transform<F>,
        mut positions: Vec<Point3<F>>,
        indices: Vec<usize>,
    ) -> Self {
        for v in &mut positions {
            *v *= object_to_world;
        }

        Self {
            object_to_world,
            position: positions,
            index: indices,
            tangent: None,
            normal: None,
            uv: None,
            alpha_mask: None,
        }
    }
}

impl<F: Number, T> TriangleMesh<F, T> {
    /// Add an alpha mask texture to a triangle mesh
    pub fn with_alpha_mask<U>(self, alpha: Arc<U>) -> TriangleMesh<F, U> {
        TriangleMesh {
            object_to_world: self.object_to_world,
            position: self.position,
            index: self.index,
            tangent: self.tangent,
            normal: self.normal,
            uv: self.uv,
            alpha_mask: Some(alpha),
        }
    }

    /// Add a set of tangent data to a triangle mesh.  Asserts that the input
    /// data has the same length as the number of vertices stored in the mesh.
    pub fn with_tangent(mut self, mut tan: Vec<Vector3<F>>) -> Self {
        assert_eq!(tan.len(), self.position.len());
        for v in &mut tan {
            *v *= self.object_to_world;
        }

        self.tangent = Some(tan);
        self
    }

    /// Add a set of normal vector data to a triangle mesh.  Asserts that the input
    /// data has the same length as the number of vertices stored in the mesh.
    pub fn with_normal(mut self, mut norm: Vec<Normal3<F>>) -> Self {
        assert_eq!(norm.len(), self.position.len());
        for v in &mut norm {
            *v *= self.object_to_world;
        }

        self.normal = Some(norm);
        self
    }

    /// Add a set of uv coordinates to a triangle mesh.  Asserts that the input
    /// data has the same length as the number of vertices stored in the mesh.
    pub fn with_uv(mut self, uv: Vec<Point2<F>>) -> Self {
        assert_eq!(uv.len(), self.position.len());
        self.uv = Some(uv);
        self
    }

    /// Get the vertex positions stored in the mesh
    pub fn positions(&self) -> &[Point3<F>] {
        &self.position
    }

    /// Get the vertex indices stored in the mesh
    pub fn indices(&self) -> &[usize] {
        &self.index
    }

    /// Get the alpha mask texture stored with the mesh
    pub fn alpha_mask(&self) -> Option<&Arc<T>> {
        self.alpha_mask.as_ref()
    }

    /// Get the vertex tangent vectors stored in the mesh
    pub fn tangent(&self) -> Option<&[Vector3<F>]> {
        self.tangent.as_deref()
    }

    /// Get the vertex normal vectors stored in the mesh
    pub fn normal(&self) -> Option<&[Normal3<F>]> {
        self.normal.as_deref()
    }

    /// Get the vertex uv locations stored in the mesh
    pub fn uv(&self) -> Option<&[Point2<F>]> {
        self.uv.as_deref()
    }

    pub fn object_to_world(&self) -> &Transform<F> {
        &self.object_to_world
    }
}
