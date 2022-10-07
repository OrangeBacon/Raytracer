use std::{ops::Range, sync::Arc};

use geometry::{Bounds3, Number, Point3, Ray, Vector3};
use rayon::{prelude::ParallelIterator, slice::ParallelSlice};

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

    /// The nodes in the BVH
    nodes: Vec<LinearBVHNode<T>>,
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

/// An index into the primitive array and the morton encoding of the primitive's centroid
#[derive(Clone, Copy)]
struct MortonPrimitive {
    /// Index into the primitive array
    primitive_index: usize,

    /// The primitive's centroid, morton encoded
    morton_code: u32,
}

/// Information about traversing the BVH after it has been constructed
#[repr(align(32))]
#[derive(Debug, Default)]
struct LinearBVHNode<T: Number> {
    /// The bounds of all children within the BVH
    bounds: Bounds3<T>,

    /// If this is a leaf node, this is the index into the primitives array for
    /// the children owned by this node.  If this is an interior node it has two
    /// children, child one is directly after this node and the offset stored here
    /// refers to the index of child two.
    child_offset: u32,

    /// The number of primitives referred to by this node.  If this number is 0,
    /// then the node is an interior node, not a leaf node.
    primitive_count: u16,

    /// The axis that the BVH partitioned along
    axis: u8,

    /// Padding to ensure that LinearBVHNode will never straddle a cache line when
    /// in a properly aligned array.
    _padding: [u8; 1],
}

/// Assert that LinearBVHNode fits in power of two size for better cache line alignment
const _: () = assert!(std::mem::size_of::<LinearBVHNode<f32>>() == 32);

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
            nodes: vec![],
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

        let root = if self.split_method == SplitMethod::HLVBH {
            self.hlbvh_build(&primitive_info, &mut total_nodes, &mut ordered_primitives)
        } else {
            self.recursive_build(
                &mut primitive_info,
                0..self.primitives.len(),
                &mut total_nodes,
                &mut ordered_primitives,
            )
        };

        self.primitives = ordered_primitives;

        // compute depth first traversal of the BVH
        self.nodes = Vec::with_capacity(total_nodes);
        self.flatten_tree(root);
    }

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

        if primitives.len() == 1 {
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
        *total_nodes += 1;

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

        for prim in primitive_info.iter() {
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
        for (i, cost) in cost.iter_mut().enumerate().take(BUCKET_COUNT - 1) {
            let mut bound_0 = Bounds3::ZERO;
            let mut bound_1 = Bounds3::ZERO;
            let mut count_0 = 0;
            let mut count_1 = 0;

            for bucket in buckets.iter().take(i + 1) {
                bound_0 = bound_0.union_box(bucket.bounds);
                count_0 += bucket.count;
            }

            for bucket in buckets.iter().take(BUCKET_COUNT).skip(i + 1) {
                bound_1 = bound_1.union_box(bucket.bounds);
                count_1 = bucket.count;
            }

            *cost = T::cast(0.125)
                + (T::cast(count_0 as i32) * bound_0.surface_area()
                    + T::cast(count_1 as i32) * bound_1.surface_area())
                    / bounds.surface_area();
        }

        // find index of minimum cost bucket
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
            Some(partition(primitive_info, |prim| {
                let mut b = BUCKET_COUNT
                    * (centroid_bounds.offset(prim.centroid)[dimension].i32() as usize);
                if b == BUCKET_COUNT {
                    b -= 1;
                }
                b <= min_cost_idx
            }))
        } else {
            None // Not worth partitioning, just create leaf node
        }
    }

    /// Build the BVH using [`SplitMethod::HLBVH`]
    fn hlbvh_build(
        &mut self,
        primitive_info: &[BVHPrimitiveInfo<T>],
        total_nodes: &mut usize,
        ordered_primitives: &mut Vec<Arc<dyn Primitive<T>>>,
    ) -> Arc<BVHBuildNode<T>> {
        // Bounds of the centroids of all primitives
        let bounds = primitive_info
            .iter()
            .fold(Bounds3::ZERO, |a, b| a.union_point(b.centroid));

        let mut morton_primitives: Vec<_> = primitive_info
            .par_chunks(512)
            .flat_map_iter(|prim| {
                prim.iter().map(|prim| {
                    let morton_bits = 10;
                    let morton_scale = 1 << morton_bits;
                    let centroid = bounds.offset(prim.centroid);
                    MortonPrimitive {
                        primitive_index: prim.primitive_number,
                        morton_code: morton_encode3(centroid.cast() * morton_scale),
                    }
                })
            })
            .collect();

        radix_sort(&mut morton_primitives);

        // find primitive ranges for each lvbh treelet
        let mut start = 0;
        let mut end = 1;
        let mut treelets_to_build = vec![];
        loop {
            let mask = 0b00111111111111000000000000000000;
            if end == morton_primitives.len()
                || ((morton_primitives[start].morton_code & mask)
                    != (morton_primitives[end].morton_code & mask))
            {
                treelets_to_build.push(start..end);
                start = end;
            }

            if end > morton_primitives.len() {
                break;
            }
            end += 1;
        }

        let mut finished_treelets = Vec::with_capacity(treelets_to_build.len());

        // TODO: Parallelise this
        for treelet in treelets_to_build {
            let first_bit_index = 29 - 12;
            let node = self.emit_lvbh(
                &primitive_info,
                &morton_primitives[treelet],
                total_nodes,
                ordered_primitives,
                first_bit_index,
            );
            finished_treelets.push(node);
        }

        self.build_upper_sah(&mut finished_treelets, total_nodes)
    }

    /// partition nodes within a hlbvh area
    fn emit_lvbh(
        &self,
        primitive_info: &[BVHPrimitiveInfo<T>],
        morton_primitives: &[MortonPrimitive],
        total_nodes: &mut usize,
        ordered_primitives: &mut Vec<Arc<dyn Primitive<T>>>,
        bit_index: i32,
    ) -> Arc<BVHBuildNode<T>> {
        assert!(morton_primitives.len() > 0);

        if bit_index == -1 || morton_primitives.len() < self.primitives_in_node {
            // create leaf node of lbvh treelet
            *total_nodes += 1;
            let mut bounds = Bounds3::default();
            let first_prim_offset = ordered_primitives.len();
            for prim in morton_primitives {
                let index = prim.primitive_index;
                ordered_primitives.push(Arc::clone(&self.primitives[index]));
                bounds = bounds.union_box(primitive_info[index].bounds);
            }
            Arc::new(BVHBuildNode::leaf(
                first_prim_offset..(first_prim_offset + morton_primitives.len()),
                bounds,
            ))
        } else {
            let mask = 1 << bit_index;

            // advance to next sub tree if this bit does not have an lvbh split
            if (morton_primitives[0].morton_code & mask)
                == (morton_primitives[morton_primitives.len() - 1].morton_code & mask)
            {
                return self.emit_lvbh(
                    primitive_info,
                    morton_primitives,
                    total_nodes,
                    ordered_primitives,
                    bit_index - 1,
                );
            }

            // find lvbh split point
            let mut start = 0;
            let mut end = morton_primitives.len() - 1;
            while start + 1 != end {
                assert_ne!(start, end);
                let mid = (start + end) / 2;
                if (morton_primitives[start].morton_code & mask)
                    == (morton_primitives[mid].morton_code & mask)
                {
                    start = mid
                } else {
                    assert_eq!(
                        (morton_primitives[mid].morton_code & mask),
                        (morton_primitives[end].morton_code & mask)
                    );
                    end = mid;
                }
            }

            let offset = end;
            assert!(offset < morton_primitives.len() - 1);
            assert_ne!(
                morton_primitives[offset - 1].morton_code & mask,
                morton_primitives[offset].morton_code & mask
            );

            // create the interior lbvh node
            *total_nodes += 1;
            let c0 = self.emit_lvbh(
                primitive_info,
                &morton_primitives[0..offset],
                total_nodes,
                ordered_primitives,
                bit_index - 1,
            );
            let c1 = self.emit_lvbh(
                primitive_info,
                &morton_primitives[offset..],
                total_nodes,
                ordered_primitives,
                bit_index - 1,
            );
            let axis = bit_index as usize % 3;

            Arc::new(BVHBuildNode::interior(axis, c0, c1))
        }
    }

    /// Create an sah BVH from LBVH treelets
    fn build_upper_sah(
        &mut self,
        treelet_roots: &mut [Arc<BVHBuildNode<T>>],
        total_nodes: &mut usize,
    ) -> Arc<BVHBuildNode<T>> {
        assert!(treelet_roots.len() != 0);

        if treelet_roots.len() == 1 {
            return Arc::clone(&treelet_roots[0]);
        }

        *total_nodes += 1;

        // compute bound of all included nodes
        let mut bounds = Bounds3::default();
        for root in treelet_roots.iter() {
            bounds = bounds.union_box(root.bounds);
        }

        // compute bound of all node centroids
        let mut centroid_bounds = Bounds3::default();
        for root in treelet_roots.iter() {
            let centroid = (root.bounds.min + root.bounds.max) * T::HALF;
            centroid_bounds = centroid_bounds.union_point(centroid);
        }
        let dim = centroid_bounds.maximum_extent();

        // pbrt lists a fix me, for what should happen when/if this assert hits
        // could continue, but the sah split would do nothing
        assert_ne!(centroid_bounds.max[dim], centroid_bounds.min[dim]);

        // initialise sah buckets
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

        for root in treelet_roots.iter() {
            let centroid = (root.bounds.min[dim] + root.bounds.max[dim]) * T::HALF;
            let mut b = BUCKET_COUNT
                * ((centroid - centroid_bounds.min[dim])
                    / (centroid_bounds.max[dim] - centroid_bounds.min[dim]))
                    .i32() as usize;
            if b == BUCKET_COUNT {
                b = BUCKET_COUNT - 1;
            }

            assert!(b < BUCKET_COUNT);
            buckets[b].count += 1;
            buckets[b].bounds = buckets[b].bounds.union_box(root.bounds);
        }

        // compute bucket split costs
        let mut cost = [T::ZERO; BUCKET_COUNT];
        for (i, cost) in cost.iter_mut().enumerate().take(BUCKET_COUNT - 1) {
            let mut bound_0 = Bounds3::ZERO;
            let mut bound_1 = Bounds3::ZERO;
            let mut count_0 = 0;
            let mut count_1 = 0;

            for bucket in buckets.iter().take(i + 1) {
                bound_0 = bound_0.union_box(bucket.bounds);
                count_0 += bucket.count;
            }

            for bucket in buckets.iter().take(BUCKET_COUNT).skip(i + 1) {
                bound_1 = bound_1.union_box(bucket.bounds);
                count_1 = bucket.count;
            }

            *cost = T::cast(0.125)
                + (T::cast(count_0 as i32) * bound_0.surface_area()
                    + T::cast(count_1 as i32) * bound_1.surface_area())
                    / bounds.surface_area();
        }

        // find index of minimum cost bucket
        let mut min_cost = cost[0];
        let mut min_cost_idx = 0;
        for (idx, &cost) in cost.iter().enumerate() {
            if cost < min_cost {
                min_cost = cost;
                min_cost_idx = idx;
            }
        }

        let mid = partition(treelet_roots, |node| {
            let centroid = (node.bounds.min[dim] + node.bounds.max[dim]) * T::HALF;

            let mut b = BUCKET_COUNT
                * ((centroid - centroid_bounds.min[dim])
                    / (centroid_bounds.max[dim] - centroid_bounds.min[dim]))
                    .i32() as usize;
            if b == BUCKET_COUNT {
                b = BUCKET_COUNT - 1;
            }

            assert!(b < BUCKET_COUNT);

            b <= min_cost_idx
        });

        assert!(mid < treelet_roots.len());

        Arc::new(BVHBuildNode::interior(
            dim,
            self.build_upper_sah(&mut treelet_roots[..mid], total_nodes),
            self.build_upper_sah(&mut treelet_roots[mid..], total_nodes),
        ))
    }

    /// Convert a bvh into a more efficient representation
    fn flatten_tree(&mut self, node: Arc<BVHBuildNode<T>>) -> usize {
        let mut linear_node = LinearBVHNode::default();
        linear_node.bounds = node.bounds;
        let idx = self.nodes.len();
        if let Some(children) = &node.children {
            // axis will only ever be 0, 1 or 2, so this will never fail
            linear_node.axis = node.split_axis.try_into().unwrap();
            linear_node.primitive_count = 0;
            self.nodes.push(linear_node);
            self.flatten_tree(Arc::clone(&children[0]));
            let second_idx = self.flatten_tree(Arc::clone(&children[1]));
            self.nodes[idx].child_offset = second_idx
                .try_into()
                .expect("Error: More than u32::MAX bvh nodes attempted to be constructed");
        } else {
            linear_node.child_offset = node
                .primitives
                .start
                .try_into()
                .expect("Error: More than u32::MAX bvh nodes attempted to be constructed");

            // max primitives per node is clamped at maximum 255, so this will never fail
            linear_node.primitive_count = node.primitives.len().try_into().unwrap();
            self.nodes.push(linear_node);
        }

        idx
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
        if let Some(node) = self.nodes.first() {
            node.bounds
        } else {
            Bounds3::default()
        }
    }

    fn intersect(&self, ray: Ray<T>) -> Option<SurfaceInteraction<T>> {
        if self.nodes.is_empty() {
            return None;
        }
        let mut hit = None;
        let inv_dir = Vector3::new(
            T::ONE / ray.direction.x,
            T::ONE / ray.direction.y,
            T::ONE / ray.direction.z,
        );
        let is_dir_negative = [
            inv_dir.x < T::ZERO,
            inv_dir.y < T::ZERO,
            inv_dir.z < T::ZERO,
        ];

        // stack of nodes to visit in the BVH
        let mut nodes_to_visit = [0; 64];
        // index of next free element in the stack
        let mut nodes_to_visit_idx = 0;
        // index of the node currently being checked
        let mut current_node_index = 0;

        loop {
            let node = &self.nodes[current_node_index];

            if node.bounds.intersect_inv(ray, inv_dir, is_dir_negative) {
                if node.primitive_count > 0 {
                    // intersect ray with primitives in leaf node
                    let start = node.child_offset as usize;
                    let range = start..(start + node.primitive_count as usize);
                    for prim in &self.primitives[range] {
                        if let Some(isect) = prim.intersect(ray) {
                            hit = Some(isect);
                        }
                    }
                    if nodes_to_visit_idx == 0 {
                        break;
                    }
                    nodes_to_visit_idx -= 1;
                    current_node_index = nodes_to_visit[nodes_to_visit_idx];
                } else {
                    // put far node on the stack
                    if is_dir_negative[node.axis as usize] {
                        nodes_to_visit[nodes_to_visit_idx] = current_node_index + 1;
                        nodes_to_visit_idx += 1;
                        current_node_index = node.child_offset as _;
                    } else {
                        nodes_to_visit[nodes_to_visit_idx] = node.child_offset as _;
                        nodes_to_visit_idx += 1;
                        current_node_index = current_node_index + 1;
                    }
                }
            } else {
                if nodes_to_visit_idx == 0 {
                    break;
                }
                nodes_to_visit_idx -= 1;
                current_node_index = nodes_to_visit[nodes_to_visit_idx];
            }
        }

        hit
    }

    fn does_intersect(&self, ray: Ray<T>) -> bool {
        if self.nodes.is_empty() {
            return false;
        }
        let inv_dir = Vector3::new(
            T::ONE / ray.direction.x,
            T::ONE / ray.direction.y,
            T::ONE / ray.direction.z,
        );
        let is_dir_negative = [
            inv_dir.x < T::ZERO,
            inv_dir.y < T::ZERO,
            inv_dir.z < T::ZERO,
        ];

        // stack of nodes to visit in the BVH
        let mut nodes_to_visit = [0; 64];
        // index of next free element in the stack
        let mut nodes_to_visit_idx = 0;
        // index of the node currently being checked
        let mut current_node_index = 0;

        loop {
            let node = &self.nodes[current_node_index];

            if node.bounds.intersect_inv(ray, inv_dir, is_dir_negative) {
                if node.primitive_count > 0 {
                    // intersect ray with primitives in leaf node
                    let start = node.child_offset as usize;
                    let range = start..(start + node.primitive_count as usize);
                    for prim in &self.primitives[range] {
                        if prim.does_intersect(ray) {
                            return true;
                        }
                    }
                    if nodes_to_visit_idx == 0 {
                        break;
                    }
                    nodes_to_visit_idx -= 1;
                    current_node_index = nodes_to_visit[nodes_to_visit_idx];
                } else {
                    // put far node on the stack
                    if is_dir_negative[node.axis as usize] {
                        nodes_to_visit[nodes_to_visit_idx] = current_node_index + 1;
                        nodes_to_visit_idx += 1;
                        current_node_index = node.child_offset as _;
                    } else {
                        nodes_to_visit[nodes_to_visit_idx] = node.child_offset as _;
                        nodes_to_visit_idx += 1;
                        current_node_index = current_node_index + 1;
                    }
                }
            } else {
                if nodes_to_visit_idx == 0 {
                    break;
                }
                nodes_to_visit_idx -= 1;
                current_node_index = nodes_to_visit[nodes_to_visit_idx];
            }
        }

        false
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

/// Decompose a number so that all of the bits are shifted 3 places apart.
fn left_shift3(mut x: u32) -> u32 {
    if x == 1 << 10 {
        x -= 1;
    }
    x = (x | (x << 16)) & 0b00000011000000000000000011111111;
    x = (x | (x << 8)) & 0b00000011000000001111000000001111;
    x = (x | (x << 4)) & 0b00000011000011000011000011000011;
    x = (x | (x << 2)) & 0b00001001001001001001001001001001;
    x
}

/// Return the 3D morton encoding of a vector
fn morton_encode3<T: Number>(v: Vector3<T>) -> u32 {
    let vec = v.cast::<i32>();
    let [x, y, z] = vec.to_array().map(|x| x as u32);
    (left_shift3(z) << 2) | (left_shift3(y) << 1) | left_shift3(x)
}

/// Specialised sorting algorithm that is better for sorting the binary numbers
/// in morton primitives.
fn radix_sort(primitives: &mut Vec<MortonPrimitive>) {
    let mut temp = Vec::with_capacity(primitives.len());
    const BITS_PER_PASS: usize = 6;
    const BIT_COUNT: usize = 30;
    const PASS_COUNT: usize = BIT_COUNT / BITS_PER_PASS;

    for pass in 0..PASS_COUNT {
        let low_bit = pass * BITS_PER_PASS;

        let (in_vec, out_vec) = if pass & 1 == 0 {
            (primitives.as_mut_slice(), temp.as_mut_slice())
        } else {
            (temp.as_mut_slice(), primitives.as_mut_slice())
        };

        // count number of zero bits in array
        const BUCKET_COUNT: usize = 1 << BITS_PER_PASS;
        let mut buckets = [0; BUCKET_COUNT];
        const BIT_MASK: usize = (1 << BITS_PER_PASS) - 1;

        for prim in &mut in_vec[..] {
            let bucket = ((prim.morton_code as usize) >> low_bit) & BIT_MASK;
            buckets[bucket] += 1;
        }

        // compute starting index for each bucket
        let mut out_index = [0; BUCKET_COUNT];
        for i in 0..BUCKET_COUNT {
            out_index[i] = out_index[i - 1] + buckets[i - 1];
        }

        // store sorted values in output array
        for prim in in_vec {
            let bucket = ((prim.morton_code as usize) >> low_bit) & BIT_MASK;
            out_vec[out_index[bucket]] = *prim;
            out_index[bucket] += 1;
        }
    }

    // copy final result from temp if needed
    if PASS_COUNT & 1 != 0 {
        std::mem::swap(primitives, &mut temp);
    }
}
