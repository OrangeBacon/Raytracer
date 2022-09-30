use std::sync::Arc;

use geometry::{log2usize, Bounds3, Number};

use crate::primitive::Primitive;

/// Binary Space Partitioning tree, splitting on only one coordinate axis
/// (octrees split on all three at the same time, so are more efficient, but harder
/// to implement)
pub struct KDTree<T: Number> {
    intersect_cost: usize,
    traversal_cost: usize,
    max_primitives: usize,

    empty_bonus: T,

    /// All primitives stored in the tree
    primitives: Vec<Arc<dyn Primitive<T>>>,

    /// All the nodes making up the tree
    nodes: Vec<KDNode<T>>,

    /// How many nodes were allocated
    allocated_nodes: usize,

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
struct KDNode<T: Number> {
    split: T,
    flags: i32,
}

enum EdgeType {
    Start,
    End,
}

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
            nodes: vec![],
            allocated_nodes: 0,
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

        let mut primitive_numbers: Vec<_> = (0..self.primitives.len()).collect();

        // self.build_tree(
        //     0,
        //     self.bounds,
        //     &primitive_bounds,
        //     &mut primitive_numbers,
        //     max_depth,
        //     edges,
        //     primitives0,
        //     primitives1,
        // );
    }

    /// Recursive tree construction algorithm
    fn build_tree(
        &mut self,
        node_id: usize,
        node_bounds: Bounds3<T>,
        primitive_bounds: &[Bounds3<T>],
        primitive_numbers: &mut [usize],
        depth: usize,
        edges: [BoundEdge<T>; 3],
        primitives0: Arc<dyn Primitive<T>>,
        primitives1: Arc<dyn Primitive<T>>,
    ) {
    }
}

impl<T: Number> KDNode<T> {
    /// Construct a node for a child containing up to `primitive_count` primitives
    fn leaf(
        primitive_numbers: &[usize],
        primitive_count: usize,
        primitive_indices: &mut Vec<usize>,
    ) -> Self {
        let split = match primitive_count {
            0 => 0,
            1 => primitive_numbers[0],
            _ => {
                let split = primitive_indices.len();
                for primitive in &primitive_numbers[..primitive_count] {
                    primitive_indices.push(*primitive);
                }
                split
            }
        };

        let split = T::from_bits(T::Bits::cast(split as u64 as i64));

        Self {
            flags: 3 | ((primitive_count as u32 as i32) << 2),
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
}
