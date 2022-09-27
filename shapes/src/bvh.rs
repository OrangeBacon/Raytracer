use std::{ops::Range, sync::Arc};

use geometry::{Bounds3, Number, Point3, Ray};

use crate::{primitive::Primitive, SurfaceInteraction};

/// The heuristic that should be used when constructing a BVH
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum SplitMethod {
    /// Surface Area Heuristic
    SAH,

    // Hierarchical Linear Volume Bounding Hierarchy
    HLVBH,

    /// Split on the centroid of every shape
    Middle,

    /// Split so that BVH nodes have equal numbers of primitives on each side
    EqualCounts,
}

/// A Bounding Volume H
#[derive(Debug)]
pub struct BVH<T: Number> {
    /// Maximum number of primitives that can be stored in one BVH node, maximum
    /// value for this is 255, anything higher is taken to be equal to 255.
    primitives_in_node: usize,

    /// How the BVH was constructed.
    split_method: SplitMethod,

    /// A list of all primitives in the BVH
    primitives: Vec<Arc<dyn Primitive<T>>>,
}

/// Information about a single primitive in the BVH
struct BVHPrimitiveInfo<T: Number> {
    /// The index into the BVH's primitive array
    primitive_number: usize,

    /// The bounding box in world space of the primitive
    bounds: Bounds3<T>,

    /// The centroid of the bounding box
    centroid: Point3<T>,
}

/// A Single node in the BVH
struct BVHBuildNode<T: Number> {
    /// The bounding box containing all the children of this node
    bounds: Bounds3<T>,

    /// Indices of the child nodes of this node
    children: Option<[Arc<Self>; 2]>,

    /// Axis the node splits along
    split_axis: usize,

    /// Indices into the BVH's primitive array that this node represents
    primitives: Range<usize>,
}

impl<T: Number> BVH<T> {
    /// Construct a BVH from a list of primitives
    pub fn new(
        primitives: Vec<Arc<dyn Primitive<T>>>,
        primitives_in_node: usize,
        split_method: SplitMethod,
    ) -> Self {
        let mut this = Self {
            primitives_in_node: primitives_in_node.min(255),
            split_method,
            primitives,
        };
        this.construct();
        this
    }

    /// Build the BVH
    fn construct(&mut self) {
        let mut primitive_info: Vec<_> = self
            .primitives
            .iter()
            .enumerate()
            .map(|(idx, prim)| BVHPrimitiveInfo::new(idx, prim.world_bound()))
            .collect();

        let mut total_nodes = 0;
        let mut ordered_primitives = Vec::with_capacity(self.primitives.len());

        if self.split_method == SplitMethod::HLVBH {
            self.hlbvh_build()
        } else {
            self.recursive_build(
                &mut primitive_info,
                0..self.primitives.len(),
                &mut total_nodes,
                &mut ordered_primitives,
            );
        }

        self.primitives = ordered_primitives;
    }

    /// Build the BVH using [`SplitMethod::HLBVH`]
    fn hlbvh_build(&mut self) {}

    /// Build the BVH using all split methods that are not HLVBH
    fn recursive_build(
        &mut self,
        primitive_info: &mut [BVHPrimitiveInfo<T>],
        primitives: Range<usize>,
        total_nodes: &mut usize,
        ordered_primitives: &mut Vec<Arc<dyn Primitive<T>>>,
    ) -> Arc<BVHBuildNode<T>> {
        // compute bounds of all included primitives
        let bounds = primitive_info
            .iter()
            .fold(Bounds3::ZERO, |a, b| a.union_box(b.bounds));

        if primitives.clone().len() == 1 {
            return self.create_leaf_node(
                primitive_info,
                primitives,
                total_nodes,
                ordered_primitives,
                bounds,
            );
        }

        let centroid_bounds = primitive_info
            .iter()
            .fold(Bounds3::ZERO, |a, b| a.union_point(b.centroid));
        let dimension = centroid_bounds.maximum_extent();

        // If the bounds have no volume, don't bother splitting and just create
        // a single leaf node for all the primitives
        if centroid_bounds.max[dimension] == centroid_bounds.min[dimension] {
            return self.create_leaf_node(
                primitive_info,
                primitives,
                total_nodes,
                ordered_primitives,
                bounds,
            );
        }

        let mid = match self.split_method {
            SplitMethod::SAH => {
                match self.sah(primitive_info, centroid_bounds, dimension, bounds) {
                    Some(mid) => mid,
                    None => {
                        return self.create_leaf_node(
                            primitive_info,
                            primitives,
                            total_nodes,
                            ordered_primitives,
                            bounds,
                        )
                    }
                }
            }
            SplitMethod::Middle => self.middle(
                &mut primitive_info[primitives.clone()],
                centroid_bounds,
                dimension,
            ),
            SplitMethod::EqualCounts => {
                self.equal_counts(&mut primitive_info[primitives.clone()], dimension)
            }
            SplitMethod::HLVBH => panic!("HLBVH construction is not valid for recursive_build"),
        };

        // recurse splitting the primitives
        let start = primitives.start;
        let end = primitives.end;
        let c0 = self.recursive_build(
            &mut primitive_info[start..mid],
            start..mid,
            total_nodes,
            ordered_primitives,
        );
        let c1 = self.recursive_build(
            &mut primitive_info[mid..end],
            mid..end,
            total_nodes,
            ordered_primitives,
        );
        let node = BVHBuildNode::interior(dimension, c0, c1);

        Arc::new(node)
    }

    /// Create a leaf node containing multiple primitives
    fn create_leaf_node(
        &mut self,
        primitive_info: &[BVHPrimitiveInfo<T>],
        primitives: Range<usize>,
        total_nodes: &mut usize,
        ordered_primitives: &mut Vec<Arc<dyn Primitive<T>>>,
        bounds: Bounds3<T>,
    ) -> Arc<BVHBuildNode<T>> {
        let start = ordered_primitives.len();
        for num in &primitive_info[primitives.clone()] {
            ordered_primitives.push(Arc::clone(&self.primitives[num.primitive_number]));
        }

        let node = BVHBuildNode::leaf(start..primitives.len(), bounds);
        *total_nodes += 1;

        Arc::new(node)
    }

    /// Partition the bvh using the 'middle' heuristic.  Returns the new mid point.
    fn middle(
        &mut self,
        primitive_info: &mut [BVHPrimitiveInfo<T>],
        centroid_bounds: Bounds3<T>,
        dimension: usize,
    ) -> usize {
        let mid = (centroid_bounds.min[dimension] + centroid_bounds.max[dimension]) / T::TWO;
        let mid = partition(primitive_info, |p| p.centroid[dimension] < mid);

        // middle partition failed, so fallback to equal
        if mid == 0 || mid == primitive_info.len() {
            self.equal_counts(primitive_info, dimension)
        } else {
            mid
        }
    }

    /// Partition the array so the number of elements on each side of the midpoint
    /// is equal and so centroids to the left of the mid point are earlier in the
    /// array, and the opposite if the centroid is to the right of the mid point.
    /// Returns the new midpoint of the primitive array
    fn equal_counts(
        &mut self,
        primitive_info: &mut [BVHPrimitiveInfo<T>],
        dimension: usize,
    ) -> usize {
        let len = primitive_info.len() / 2;
        primitive_info.select_nth_unstable_by(len, |a, b| {
            a.centroid[dimension]
                .partial_cmp(&b.centroid[dimension])
                .unwrap_or(std::cmp::Ordering::Less)
        });

        len
    }

    /// Partition the array using the surface area heuristic.
    /// If the split succeeds, returns the index of the point to recurse around,
    /// Otherwise, returns None and a leaf node should be created.
    fn sah(
        &mut self,
        primitive_info: &mut [BVHPrimitiveInfo<T>],
        centroid_bounds: Bounds3<T>,
        dimension: usize,
        bounds: Bounds3<T>,
    ) -> Option<usize> {
        if primitive_info.len() <= 4 {
            return Some(self.equal_counts(primitive_info, dimension));
        }

        #[derive(Clone, Copy)]
        struct BucketInfo<T: Number> {
            count: usize,
            bounds: Bounds3<T>,
        }

        const BUCKET_COUNT: usize = 12;

        let mut buckets = [BucketInfo {
            count: 0,
            bounds: Bounds3::ZERO,
        }; BUCKET_COUNT];

        for prim in &primitive_info[..] {
            let mut b =
                BUCKET_COUNT * (centroid_bounds.offset(prim.centroid)[dimension].i32() as usize);
            if b == BUCKET_COUNT {
                b -= 1;
            }
            buckets[b].count += 1;
            buckets[b].bounds = buckets[b].bounds.union_box(prim.bounds);
        }

        // compute bucket split costs
        let mut cost = [T::ZERO; BUCKET_COUNT];
        for i in 0..BUCKET_COUNT - 1 {
            let mut bound_0 = Bounds3::ZERO;
            let mut bound_1 = Bounds3::ZERO;
            let mut count_0 = 0;
            let mut count_1 = 0;

            for j in 0..=i {
                bound_0 = bound_0.union_box(buckets[j].bounds);
                count_0 += buckets[j].count;
            }

            for j in i + 1..BUCKET_COUNT {
                bound_1 = bound_1.union_box(buckets[j].bounds);
                count_1 = buckets[j].count;
            }

            cost[i] = T::cast(0.125)
                + (T::cast(count_0 as i32) * bound_0.surface_area()
                    + T::cast(count_1 as i32) * bound_1.surface_area())
                    / bounds.surface_area();
        }

        let mut min_cost = cost[0];
        let mut min_cost_idx = 0;
        for (idx, &cost) in cost.iter().enumerate() {
            if cost < min_cost {
                min_cost = cost;
                min_cost_idx = idx;
            }
        }

        if primitive_info.len() > self.primitives_in_node
            || (min_cost.i32() as usize) < primitive_info.len()
        {
            return Some(partition(primitive_info, |prim| {
                let mut b = BUCKET_COUNT
                    * (centroid_bounds.offset(prim.centroid)[dimension].i32() as usize);
                if b == BUCKET_COUNT {
                    b -= 1;
                }
                b <= min_cost_idx
            }));
        } else {
            None
        }
    }
}

impl<T: Number> BVHPrimitiveInfo<T> {
    /// Construct new information about a primitive
    fn new(primitive_number: usize, bounds: Bounds3<T>) -> Self {
        Self {
            primitive_number,
            bounds,
            centroid: bounds.min * T::HALF + bounds.max * T::HALF,
        }
    }
}

impl<T: Number> BVHBuildNode<T> {
    /// Construct a new Leaf node (no further nodes, only primitive children)
    fn leaf(primitives: Range<usize>, bounds: Bounds3<T>) -> Self {
        Self {
            bounds,
            children: None,
            split_axis: 0,
            primitives,
        }
    }

    /// Build an interior node of the BVH (No primitive children, both next nodes)
    fn interior(axis: usize, c0: Arc<Self>, c1: Arc<Self>) -> Self {
        Self {
            bounds: c0.bounds.union_box(c1.bounds),
            children: Some([c0, c1]),
            split_axis: axis,
            primitives: 0..0,
        }
    }
}

impl<T: Number> Primitive<T> for BVH<T> {
    fn world_bound(&self) -> Bounds3<T> {
        todo!()
    }

    fn intersect(&self, ray: Ray<(), T>) -> Option<SurfaceInteraction<(), T>> {
        todo!()
    }

    fn does_intersect(&self, ray: Ray<(), T>) -> bool {
        todo!()
    }

    fn area_light(&self) {
        panic!("<BVH as Primitive>::area_light should never be called");
    }

    fn material(&self) {
        panic!("<BVH as Primitive>::material should never be called");
    }

    fn scattering_functions(&self) {
        panic!("<BVH as Primitive>::scattering_functions should never be called");
    }
}

/// Sort the input array so that all elements that return true are at the start
/// and all elements that return false are at the end.  Returns the index of the
/// first element that returned false.
fn partition<T>(array: &mut [T], f: impl Fn(&T) -> bool) -> usize {
    let len = array.len();
    if len == 0 {
        return 0;
    }

    let (mut left, mut right) = (0, len - 1);
    loop {
        while left < len && f(&array[left]) {
            left += 1;
        }
        while right > 0 && !f(&array[right]) {
            right -= 1;
        }
        if left >= right {
            break left;
        }
        array.swap(left, right);
    }
}
