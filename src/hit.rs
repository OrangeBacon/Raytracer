use std::{f64::consts::PI, fmt::Debug, sync::Arc};

use glam::{dvec2, dvec3, DVec2, DVec3};

use crate::{aabb::Aabb, material::Material, ray::Ray};

pub struct HitRecord {
    pub point: DVec3,
    pub normal: DVec3,
    pub t: f64,
    pub front_face: bool,
    pub material: Arc<dyn Material>,
    pub uv: DVec2,
}

pub trait Hittable: Send + Sync + Debug {
    fn hit(&self, ray: Ray, t_min: f64, t_max: f64) -> Option<HitRecord>;

    fn bounding_box(&self, time0: f64, time1: f64) -> Option<Aabb>;
}

impl HitRecord {
    pub fn new(
        ray: Ray,
        point: DVec3,
        outward_normal: DVec3,
        t: f64,
        material: Arc<dyn Material>,
        uv: DVec2,
    ) -> Self {
        let mut record = HitRecord {
            point,
            material,
            normal: DVec3::ZERO,
            t,
            front_face: false,
            uv,
        };
        record.set_face_normal(ray, outward_normal);
        record
    }

    fn set_face_normal(&mut self, ray: Ray, outward_normal: DVec3) {
        self.front_face = ray.direction.dot(outward_normal) < 0.0;
        self.normal = if self.front_face {
            outward_normal
        } else {
            -outward_normal
        }
    }
}

#[derive(Clone, Debug)]
pub struct Sphere {
    pub centre: DVec3,
    pub radius: f64,
    pub material: Arc<dyn Material>,
}

impl Hittable for Sphere {
    fn hit(&self, ray: Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let oc = ray.origin - self.centre;
        let a = ray.direction.length_squared();
        let half_b = oc.dot(ray.direction);
        let c = oc.length_squared() - self.radius * self.radius;
        let discriminant = half_b * half_b - a * c;
        if discriminant < 0.0 {
            return None;
        }

        let sqrt_discriminant = discriminant.sqrt();

        let mut root = (-half_b - sqrt_discriminant) / a;
        if root < t_min || t_max < root {
            root = (-half_b + sqrt_discriminant) / a;
            if root < t_min || t_max < root {
                return None;
            }
        }

        let point = ray.at(root);
        let outward_normal = (point - self.centre) / self.radius;
        Some(HitRecord::new(
            ray,
            point,
            outward_normal,
            root,
            Arc::clone(&self.material),
            sphere_uv(outward_normal),
        ))
    }

    fn bounding_box(&self, _time0: f64, _time1: f64) -> Option<Aabb> {
        Some(Aabb {
            minimum: self.centre - DVec3::splat(self.radius),
            maximum: self.centre + DVec3::splat(self.radius),
        })
    }
}

fn sphere_uv(point: DVec3) -> DVec2 {
    let theta = (-point.y).acos();
    let phi = (-point.z).atan2(point.x) + PI;

    dvec2(phi / (2.0 * PI), theta / PI)
}

#[cfg(test)]
#[test]
fn sphere_uv_test() {
    use approx::assert_relative_eq;
    use glam::dvec3;

    assert_relative_eq!(sphere_uv(dvec3(1.0, 0.0, 0.0)), dvec2(0.5, 0.5));
    assert_relative_eq!(sphere_uv(dvec3(0.0, 1.0, 0.0)), dvec2(0.5, 1.0));
    assert_relative_eq!(sphere_uv(dvec3(0.0, 0.0, 1.0)), dvec2(0.25, 0.5));
    assert_relative_eq!(sphere_uv(dvec3(-1.0, 0.0, 0.0)), dvec2(0.0, 0.5));
    assert_relative_eq!(sphere_uv(dvec3(0.0, -1.0, 0.0)), dvec2(0.5, 0.0));
    assert_relative_eq!(sphere_uv(dvec3(0.0, 0.0, -1.0)), dvec2(0.75, 0.5));
}

#[derive(Debug)]
pub struct MovingSphere {
    pub centre0: DVec3,
    pub centre1: DVec3,
    pub time0: f64,
    pub time1: f64,
    pub radius: f64,
    pub material: Arc<dyn Material>,
}

impl MovingSphere {
    fn centre(&self, time: f64) -> DVec3 {
        self.centre0
            + ((time - self.time0) / (self.time1 - self.time0)) * (self.centre1 - self.centre0)
    }
}

impl Hittable for MovingSphere {
    fn hit(&self, ray: Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let oc = ray.origin - self.centre(ray.time);
        let a = ray.direction.length_squared();
        let half_b = oc.dot(ray.direction);
        let c = oc.length_squared() - self.radius * self.radius;
        let discriminant = half_b * half_b - a * c;
        if discriminant < 0.0 {
            return None;
        }

        let sqrt_discriminant = discriminant.sqrt();

        let mut root = (-half_b - sqrt_discriminant) / a;
        if root < t_min || t_max < root {
            root = (-half_b + sqrt_discriminant) / a;
            if root < t_min || t_max < root {
                return None;
            }
        }

        let point = ray.at(root);
        let outward_normal = (point - self.centre(ray.time)) / self.radius;
        Some(HitRecord::new(
            ray,
            point,
            outward_normal,
            root,
            Arc::clone(&self.material),
            sphere_uv(outward_normal),
        ))
    }

    fn bounding_box(&self, time0: f64, time1: f64) -> Option<Aabb> {
        let box0 = Aabb {
            minimum: self.centre(time0) - DVec3::splat(self.radius),
            maximum: self.centre(time0) + DVec3::splat(self.radius),
        };

        let box1 = Aabb {
            minimum: self.centre(time1) - DVec3::splat(self.radius),
            maximum: self.centre(time1) + DVec3::splat(self.radius),
        };

        Some(box0.surrounds(box1))
    }
}

#[derive(Debug)]
pub struct HittableVec {
    pub objects: Vec<Arc<dyn Hittable>>,
}

impl Hittable for HittableVec {
    fn hit(&self, ray: Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let mut closest = t_max;
        let mut output = None;

        for object in &self.objects {
            if let Some(record) = object.hit(ray, t_min, closest) {
                closest = record.t;
                output = Some(record);
            }
        }

        output
    }

    fn bounding_box(&self, time0: f64, time1: f64) -> Option<Aabb> {
        let mut result = None;

        for object in &self.objects {
            if let Some(aabb) = object.bounding_box(time0, time1) {
                match result {
                    None => result = Some(aabb),
                    Some(a) => result = Some(a.surrounds(aabb)),
                }
            } else {
                return None;
            }
        }

        result
    }
}

#[derive(Debug)]
pub struct XyRect {
    pub mat: Arc<dyn Material>,
    pub x0: f64,
    pub x1: f64,
    pub y0: f64,
    pub y1: f64,
    pub k: f64,
}

impl Hittable for XyRect {
    fn hit(&self, ray: Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let t = (self.k - ray.origin.z) / ray.direction.z;
        if t < t_min || t > t_max {
            return None;
        }

        let x = ray.origin.x + t * ray.direction.x;
        let y = ray.origin.y + t * ray.direction.y;
        if x < self.x0 || x > self.x1 || y < self.y0 || y > self.y1 {
            return None;
        }

        Some(HitRecord::new(
            ray,
            ray.at(t),
            DVec3::Z,
            t,
            Arc::clone(&self.mat),
            dvec2(
                (x - self.x0) / (self.x1 - self.x0),
                (y - self.y0) / (self.y1 - self.y0),
            ),
        ))
    }

    fn bounding_box(&self, _: f64, _: f64) -> Option<Aabb> {
        Some(Aabb {
            minimum: dvec3(self.x0, self.y0, self.k - 0.0001),
            maximum: dvec3(self.x1, self.y1, self.k + 0.0001),
        })
    }
}

#[derive(Debug)]
pub struct XzRect {
    pub mat: Arc<dyn Material>,
    pub x0: f64,
    pub x1: f64,
    pub z0: f64,
    pub z1: f64,
    pub k: f64,
}

impl Hittable for XzRect {
    fn hit(&self, ray: Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let t = (self.k - ray.origin.y) / ray.direction.y;
        if t < t_min || t > t_max {
            return None;
        }

        let x = ray.origin.x + t * ray.direction.x;
        let z = ray.origin.z + t * ray.direction.z;
        if x < self.x0 || x > self.x1 || z < self.z0 || z > self.z1 {
            return None;
        }

        Some(HitRecord::new(
            ray,
            ray.at(t),
            DVec3::Y,
            t,
            Arc::clone(&self.mat),
            dvec2(
                (x - self.x0) / (self.x1 - self.x0),
                (z - self.z0) / (self.z1 - self.z0),
            ),
        ))
    }

    fn bounding_box(&self, _: f64, _: f64) -> Option<Aabb> {
        Some(Aabb {
            minimum: dvec3(self.x0, self.k - 0.0001, self.z0),
            maximum: dvec3(self.x1, self.k + 0.0001, self.z1),
        })
    }
}

#[derive(Debug)]
pub struct YzRect {
    pub mat: Arc<dyn Material>,
    pub y0: f64,
    pub y1: f64,
    pub z0: f64,
    pub z1: f64,
    pub k: f64,
}

impl Hittable for YzRect {
    fn hit(&self, ray: Ray, t_min: f64, t_max: f64) -> Option<HitRecord> {
        let t = (self.k - ray.origin.x) / ray.direction.x;
        if t < t_min || t > t_max {
            return None;
        }

        let y = ray.origin.y + t * ray.direction.y;
        let z = ray.origin.z + t * ray.direction.z;
        if y < self.y0 || y > self.y1 || z < self.z0 || z > self.z1 {
            return None;
        }

        Some(HitRecord::new(
            ray,
            ray.at(t),
            DVec3::X,
            t,
            Arc::clone(&self.mat),
            dvec2(
                (y - self.y0) / (self.y1 - self.y0),
                (z - self.z0) / (self.z1 - self.z0),
            ),
        ))
    }

    fn bounding_box(&self, _: f64, _: f64) -> Option<Aabb> {
        Some(Aabb {
            minimum: dvec3(self.k - 0.0001, self.y0, self.z0),
            maximum: dvec3(self.k + 0.0001, self.y1, self.z1),
        })
    }
}
