use glam::DVec3;
use rand::Rng;

pub fn rand_sphere_point() -> DVec3 {
    let mut rng = rand::thread_rng();

    let point = DVec3::from_array(rng.gen());

    if point == DVec3::ZERO {
        DVec3::ZERO
    } else {
        point.length_recip() * point
    }
}

pub fn rand_hemisphere_point(normal: DVec3) -> DVec3 {
    let rand = rand_sphere_point();
    if rand.dot(normal) > 0.0 {
        rand
    } else {
        -rand
    }
}

pub fn near_zero(vec: DVec3) -> bool {
    let e = 1.0e-8;
    vec.x.abs() < e && vec.y.abs() < e && vec.x.abs() < e
}