use glam::{dvec3, DVec3};
use rand::Rng;

pub fn rand_sphere_point() -> DVec3 {
    let mut rng = rand::thread_rng();

    loop {
        let p = dvec3(
            rng.gen_range(-1.0..=1.0),
            rng.gen_range(-1.0..=1.0),
            rng.gen_range(-1.0..=1.0),
        );
        if p.length_squared() < 1.0 {
            break p;
        }
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

pub fn rand_unit_disk() -> DVec3 {
    let mut rng = rand::thread_rng();

    loop {
        let p = dvec3(rng.gen_range(-1.0..=1.0), rng.gen_range(-1.0..=1.0), 0.0);
        if p.length_squared() < 1.0 {
            break p;
        }
    }
}

pub fn near_zero(vec: DVec3) -> bool {
    let e = 1.0e-8;
    vec.x.abs() < e && vec.y.abs() < e && vec.x.abs() < e
}
