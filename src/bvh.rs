use std::{cmp::Ordering, sync::Arc};

use rand::Rng;

use crate::{
    aabb::Aabb,
    shapes::{HitRecord, Hittable},
};

#[derive(Debug)]
pub struct BVHNode {
    left: Arc<dyn Hittable>,
    right: Arc<dyn Hittable>,
    aabb: Aabb,
}

impl BVHNode {
    pub fn new(objects: &[Arc<dyn Hittable>], time0: f64, time1: f64) -> Option<Arc<Self>> {
        let axis = rand::thread_rng().gen_range(0..=2);

        let (left, right) = if objects.len() == 1 {
            (Arc::clone(&objects[0]), Arc::clone(&objects[0]))
        } else if objects.len() == 2 {
            if Self::compare_box(&objects[0], &objects[1], axis)?.is_lt() {
                (Arc::clone(&objects[0]), Arc::clone(&objects[1]))
            } else {
                (Arc::clone(&objects[1]), Arc::clone(&objects[0]))
            }
        } else {
            let mut sorted = objects.to_vec();
            sorted.sort_by(|a, b| Self::compare_box(a, b, axis).unwrap());

            let res: (Arc<dyn Hittable>, Arc<dyn Hittable>) = (
                Self::new(&sorted[..sorted.len() / 2], time0, time1)?,
                Self::new(&sorted[sorted.len() / 2..], time0, time1)?,
            );

            res
        };

        let left_box = left.bounding_box(time0, time1)?;
        let right_box = right.bounding_box(time0, time1)?;

        Some(Arc::new(BVHNode {
            left,
            right,
            aabb: left_box.surrounds(right_box),
        }))
    }

    fn compare_box(a: &Arc<dyn Hittable>, b: &Arc<dyn Hittable>, axis: usize) -> Option<Ordering> {
        let box_a = a.bounding_box(0.0, 0.0)?;
        let box_b = b.bounding_box(0.0, 0.0)?;

        box_a.minimum[axis].partial_cmp(&box_b.minimum[axis])
    }
}

impl Hittable for BVHNode {
    fn hit(&self, ray: crate::ray::Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        if !self.aabb.hit(ray, t_min, t_max) {
            return None;
        }

        let left = self.left.hit(ray, t_min, t_max);
        let right = self
            .right
            .hit(ray, t_min, left.as_ref().map(|x| x.t).unwrap_or(t_max));

        right.or(left)
    }

    fn bounding_box(&self, _time0: f64, _time1: f64) -> Option<Aabb> {
        Some(self.aabb)
    }
}
