use std::sync::Arc;

use geometry::{log2usize, Bounds3, Number, Ray, Vector3};

use crate::{primitive::Primitive, SurfaceInteraction};

/// Binary Space Partitioning tree, splitting on only one coordinate axis
/// (octrees split on all three at the same time, so are more efficient, but harder
/// to implement)
#[derive(Debug)]
pub struct KDTree<T: Number> {
    intersect_cost: usize,
    traversal_cost: usize,
    max_primitives: usize,

    empty_bonus: T,

    /// All primitives stored in the tree
    primitives: Vec<Arc<dyn Primitive<T>>>,

    /// List of indexes into the primitives array
    primitive_indices: Vec<usize>,

    /// All the nodes making up the tree
    nodes: Vec<KDNode<T>>,

    /// The index of the next un-used node in the nodes vec
    next_free_node: usize,

    /// The bounding box of the tree
    bounds: Bounds3<T>,
}

/// A single node in the KD tree.  The 2 low order bits of flags are the kind of
/// node, [0,1,2] => interior nodes with splits in the [x,y,z] axis and 3 => leaf
/// nodes.
///
/// If it is a leaf node, the remaining 30 bits of flags are the number
/// of primitives.
///
/// If it is an interior node, then split is a float recording how far along the
/// split axis the split is being performed.  For the two children of the node,
/// the first is immediately after this node in the final array and the second
/// one's index is stored in the remaining 30 bits of `flags`.
#[repr(align(8))]
#[derive(Debug, Default)]
struct KDNode<T: Number> {
    split: T,
    flags: i32,
}

#[derive(Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
enum EdgeType {
    #[default]
    Start,
    End,
}

#[derive(Debug, Default)]
struct BoundEdge<T: Number> {
    t: T,
    primitive_number: usize,
    edge_type: EdgeType,
}

const _: () = assert!(std::mem::size_of::<KDNode<f32>>() == 8);

impl<T: Number> KDTree<T> {
    /// Construct a new KD tree
    pub fn new(
        primitives: Vec<Arc<dyn Primitive<T>>>,
        intersect_cost: usize,
        traversal_cost: usize,
        empty_bonus: T,
        max_primitives: usize,
        max_depth: usize,
    ) -> Self {
        let mut tree = Self {
            intersect_cost,
            traversal_cost,
            max_primitives,
            empty_bonus,
            primitives,
            primitive_indices: vec![],
            nodes: vec![],
            next_free_node: 0,
            bounds: Bounds3::default(),
        };

        tree.build(max_depth);

        tree
    }

    /// Build the KD tree (part of the constructor)
    fn build(&mut self, mut max_depth: usize) {
        if max_depth <= 0 {
            max_depth = (8.0 + 1.3 * (log2usize(self.primitives.len()) as f64)).round() as usize;
        }

        let primitive_bounds: Vec<_> = self.primitives.iter().map(|p| p.world_bound()).collect();
        let primitive_numbers: Vec<_> = (0..self.primitives.len()).collect();

        let mut bound_edges = [(); 3].map(|_| Vec::with_capacity(2 * self.primitives.len()));

        self.build_tree(
            self.bounds,
            &primitive_bounds,
            &primitive_numbers,
            max_depth,
            &mut bound_edges,
            0,
        );
    }

    /// Recursive tree construction algorithm
    fn build_tree(
        &mut self,
        node_bounds: Bounds3<T>,
        primitive_bounds: &[Bounds3<T>],
        primitive_numbers: &[usize],
        depth: usize,
        edges: &mut [Vec<BoundEdge<T>>; 3],
        mut bad_refines: usize,
    ) {
        let primitive_count = primitive_numbers.len();

        // get free node from nodes array
        self.next_free_node += 1;

        // initialise leaf node
        if primitive_count <= self.max_primitives || depth == 0 {
            self.nodes
                .push(KDNode::leaf(primitive_numbers, &mut self.primitive_indices))
        }

        // choose split axis for interior node
        let mut best_axis = usize::MAX;
        let mut best_offset = usize::MAX;
        let mut best_cost = T::INFINITY;
        let old_cost = T::cast((self.intersect_cost * primitive_count) as i64);
        let total_surface = node_bounds.surface_area();
        let inv_total_surface = T::ONE / total_surface;
        let d = node_bounds.max - node_bounds.min;

        for axis_index in 0..3 {
            let axis = (node_bounds.maximum_extent() + axis_index) % 3;

            // initialise edges for axis
            for (i, &prim) in primitive_numbers.iter().enumerate() {
                let bounds = primitive_bounds[prim];
                edges[axis][2 * i] = BoundEdge::new(bounds.min[axis], prim, true);
                edges[axis][2 * i + 1] = BoundEdge::new(bounds.max[axis], prim, false);
            }

            edges[axis].sort_by(|a, b| match a.t.partial_cmp(&b.t).unwrap() {
                a @ (std::cmp::Ordering::Less | std::cmp::Ordering::Greater) => a,
                std::cmp::Ordering::Equal => a.edge_type.cmp(&b.edge_type),
            });

            // compute split cost
            let mut count_below = 0;
            let mut count_above = primitive_count;
            for i in 0..2 * primitive_count {
                if edges[axis][i].edge_type == EdgeType::End {
                    count_above -= 1;
                }

                let edge_t = edges[axis][i].t;
                if edge_t > node_bounds.min[axis] && edge_t < node_bounds.max[axis] {
                    // Compute split cost at ith edge

                    let other_axis0 = (axis + 1) % 3;
                    let other_axis1 = (axis + 2) % 3;
                    let below_surface = T::TWO
                        * (d[other_axis0] * d[other_axis1]
                            + (edge_t - node_bounds.min[axis]) * (d[other_axis0] + d[other_axis1]));
                    let above_surface = T::TWO
                        * (d[other_axis0] * d[other_axis1]
                            + (edge_t - node_bounds.max[axis]) * (d[other_axis0] + d[other_axis1]));

                    let p_below = below_surface * inv_total_surface;
                    let p_above = above_surface * inv_total_surface;
                    let empty_bonus = if count_above == 0 || count_below == 0 {
                        self.empty_bonus
                    } else {
                        T::ZERO
                    };

                    let cost = T::cast(self.traversal_cost as i64)
                        + T::cast(self.intersect_cost as i64)
                            * (T::ONE - empty_bonus)
                            * (p_below * T::cast(count_below as i64)
                                + p_above * T::cast(count_above as i64));

                    if cost < best_cost {
                        best_cost = cost;
                        best_axis = axis;
                        best_offset = i;
                    }
                }

                if edges[axis][i].edge_type == EdgeType::Start {
                    count_below += 1;
                }
            }
        }

        if best_cost > old_cost {
            bad_refines += 1;
        }

        if (best_cost > T::cast(4) * old_cost && primitive_count < 16)
            || best_axis == usize::MAX
            || bad_refines == 3
        {
            self.nodes
                .push(KDNode::leaf(primitive_numbers, &mut self.primitive_indices));
            return;
        }

        // classify primitives with respect to split
        // TODO: reuse these vectors rather than re-allocating
        let mut primitives0 = vec![];
        let mut primitives1 = vec![];
        for i in 0..best_offset {
            if edges[best_axis][i].edge_type == EdgeType::Start {
                primitives0.push(edges[best_axis][i].primitive_number);
            }
        }
        for i in best_offset + 1..2 * primitive_count {
            if edges[best_axis][i].edge_type == EdgeType::End {
                primitives1.push(edges[best_axis][i].primitive_number);
            }
        }

        // recursively initialise child nodes
        let t_split = edges[best_axis][best_offset].t;
        let mut bounds0 = node_bounds;
        let mut bounds1 = node_bounds;
        bounds0.max[best_axis] = t_split;
        bounds1.min[best_axis] = t_split;
        let node_num = self.nodes.len();

        // placeholder item
        self.nodes.push(KDNode::default());

        self.build_tree(
            bounds0,
            primitive_bounds,
            &primitives0,
            depth - 1,
            edges,
            bad_refines,
        );
        self.nodes[node_num] = KDNode::interior(best_axis, self.nodes.len(), t_split);
        self.build_tree(
            bounds1,
            primitive_bounds,
            &primitives1,
            depth - 1,
            edges,
            bad_refines,
        )
    }
}

impl<T: Number> KDNode<T> {
    /// Construct a node for a child
    fn leaf(primitive_numbers: &[usize], primitive_indices: &mut Vec<usize>) -> Self {
        let split = match primitive_numbers.len() {
            0 => 0,
            1 => primitive_numbers[0],
            _ => {
                let split = primitive_indices.len();
                for primitive in primitive_numbers {
                    primitive_indices.push(*primitive);
                }
                split
            }
        };

        let split = T::from_bits(T::Bits::cast(split as u64 as i64));

        Self {
            flags: 3 | ((primitive_numbers.len() as u32 as i32) << 2),
            split,
        }
    }

    /// Construct a KD node for an interior node with a split position along the
    /// given axis and an index to the next child
    fn interior(axis: usize, next_child: usize, split_position: T) -> Self {
        Self {
            split: split_position,
            flags: (axis as u32 as i32) | ((next_child as u32 as i32) << 2),
        }
    }

    /// Position of a split along the given axis, returns unspecified float if
    /// this node is not an interior node
    fn split_position(&self) -> T {
        self.split
    }

    /// Number of primitives stored in the node, unspecified if not a leaf node
    fn primitive_count(&self) -> usize {
        (self.flags >> 2) as u32 as usize
    }

    /// The axis that was used to split the primitives
    fn split_axis(&self) -> usize {
        (self.flags & 3) as u32 as usize
    }

    /// Is this node a leaf node
    fn is_leaf(&self) -> bool {
        (self.flags & 3) == 3
    }

    /// Index of the next child, unspecified if not an interior node
    fn next_child(&self) -> usize {
        (self.flags >> 2) as u32 as usize
    }

    /// Offset of a primitive in the node, unspecified if not a leaf node
    fn primitive_offset(&self) -> usize {
        self.split.to_bits().i64() as usize
    }
}

impl<T: Number> BoundEdge<T> {
    fn new(t: T, primitive_number: usize, starting: bool) -> Self {
        Self {
            t,
            primitive_number,
            edge_type: if starting {
                EdgeType::Start
            } else {
                EdgeType::End
            },
        }
    }
}

impl<T: Number> Primitive<T> for KDTree<T> {
    fn world_bound(&self) -> Bounds3<T> {
        self.bounds
    }

    fn intersect(&self, ray: Ray<T>) -> Option<SurfaceInteraction<T>> {
        let (mut t_min, mut t_max) = self.bounds.intersect_p(ray)?;

        #[derive(Clone, Copy, Default)]
        struct KdTodo<T: Number> {
            node_idx: usize,
            t_min: T,
            t_max: T,
        }

        let inv_dir = Vector3::new(
            T::ONE / ray.direction.x,
            T::ONE / ray.direction.y,
            T::ONE / ray.direction.z,
        );

        let mut todo = [KdTodo::default(); 64];
        let mut todo_pos = 0;
        let mut isect = None;

        let mut node_idx = 0;

        loop {
            let node = &self.nodes[node_idx];

            if ray.t_max < t_min {
                break;
            }

            if node.is_leaf() {
                let primitive_count = node.primitive_count();
                if primitive_count == 1 {
                    let p = &self.primitives[node.primitive_offset()];
                    isect = p.intersect(ray);
                } else {
                    for i in 0..primitive_count {
                        let index = self.primitive_indices[node.primitive_offset() + i];
                        let p = &self.primitives[index];
                        isect = p.intersect(ray);
                    }
                }

                // get next node to process
                if todo_pos > 0 {
                    todo_pos -= 1;
                    KdTodo {
                        node_idx,
                        t_min,
                        t_max,
                    } = todo[todo_pos];
                } else {
                    break;
                }
            } else {
                let axis = node.split_axis();
                let plane = (node.split_position() - ray.origin[axis]) * inv_dir[axis];
                let below_first = (ray.origin[axis] < node.split_position())
                    || (ray.origin[axis] == node.split_position()
                        && ray.direction[axis] <= T::ZERO);
                let (first, second) = if below_first {
                    (node_idx + 1, node.next_child())
                } else {
                    (node.next_child(), node_idx + 1)
                };

                if plane > t_max || plane <= T::ZERO {
                    node_idx = first;
                } else if plane < t_min {
                    node_idx = second;
                } else {
                    todo[todo_pos] = KdTodo {
                        node_idx: second,
                        t_min: plane,
                        t_max,
                    };
                    todo_pos += 1;
                    node_idx = first;
                    t_max = plane;
                }
            }
        }

        isect
    }

    fn does_intersect(&self, ray: Ray<T>) -> bool {
        let (mut t_min, mut t_max) = match self.bounds.intersect_p(ray) {
            Some(a) => a,
            None => return false,
        };

        #[derive(Clone, Copy, Default)]
        struct KdTodo<T: Number> {
            node_idx: usize,
            t_min: T,
            t_max: T,
        }

        let inv_dir = Vector3::new(
            T::ONE / ray.direction.x,
            T::ONE / ray.direction.y,
            T::ONE / ray.direction.z,
        );

        let mut todo = [KdTodo::default(); 64];
        let mut todo_pos = 0;

        let mut node_idx = 0;

        loop {
            let node = &self.nodes[node_idx];

            if ray.t_max < t_min {
                break;
            }

            if node.is_leaf() {
                let primitive_count = node.primitive_count();
                if primitive_count == 1 {
                    let p = &self.primitives[node.primitive_offset()];
                    if p.does_intersect(ray) {
                        return true;
                    }
                } else {
                    for i in 0..primitive_count {
                        let index = self.primitive_indices[node.primitive_offset() + i];
                        let p = &self.primitives[index];
                        if p.does_intersect(ray) {
                            return true;
                        }
                    }
                }

                // get next node to process
                if todo_pos > 0 {
                    todo_pos -= 1;
                    KdTodo {
                        node_idx,
                        t_min,
                        t_max,
                    } = todo[todo_pos];
                } else {
                    break;
                }
            } else {
                let axis = node.split_axis();
                let plane = (node.split_position() - ray.origin[axis]) * inv_dir[axis];
                let below_first = (ray.origin[axis] < node.split_position())
                    || (ray.origin[axis] == node.split_position()
                        && ray.direction[axis] <= T::ZERO);
                let (first, second) = if below_first {
                    (node_idx + 1, node.next_child())
                } else {
                    (node.next_child(), node_idx + 1)
                };

                if plane > t_max || plane <= T::ZERO {
                    node_idx = first;
                } else if plane < t_min {
                    node_idx = second;
                } else {
                    todo[todo_pos] = KdTodo {
                        node_idx: second,
                        t_min: plane,
                        t_max,
                    };
                    todo_pos += 1;
                    node_idx = first;
                    t_max = plane;
                }
            }
        }

        true
    }

    fn area_light(&self) {
        panic!("<KDTree as Primitive>::area_light should never be called");
    }

    fn material(&self) {
        panic!("<KDTree as Primitive>::area_light should never be called");
    }

    fn scattering_functions(&self) {
        panic!("<KDTree as Primitive>::area_light should never be called");
    }
}
