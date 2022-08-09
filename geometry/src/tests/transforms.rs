use crate::{tests::rng::Rng, AnimatedTransform, Bounds3f, Point3f, Transform, Vector3f};

#[test]
#[ignore]
fn random_transform_interpolation() {
    let mut rng = Rng::new(1);

    for _ in 0..200 {
        // get random transform
        let t0 = random_transform(&mut rng);
        let t1 = random_transform(&mut rng);
        let animation = AnimatedTransform::new(t0, 0.0, t1, 1.0).unwrap();

        for _ in 0..5 {
            // get random bounding box, moving through the random transform
            let mut r = || -10.0 + 20.0 * rng.float();
            let bounds = Bounds3f::new(Point3f::new(r(), r(), r()), Point3f::new(r(), r(), r()));
            let motion = animation.motion_bounds(bounds);

            let mut t = 0.0;
            while t < 1.0 {
                // Interpolate at random times and get the location of the bounds at that time
                let tr = animation.interpolate(t);
                let mut tb = bounds * tr;

                // allow for floating point error
                tb.min += 1.0e-4 * tb.diagonal();
                tb.max -= 1.0e-4 * tb.diagonal();

                // transformed bounds should be within motion bounds
                assert!(tb.min.x >= motion.min.x);
                assert!(tb.max.x <= motion.max.x);

                assert!(tb.min.y >= motion.min.y);
                assert!(tb.max.y <= motion.max.y);

                assert!(tb.min.z >= motion.min.z);
                assert!(tb.max.z <= motion.max.z);

                t += 1.0e-3 * rng.float();
            }
        }
    }
}

fn random_transform(rng: &mut Rng) -> Transform {
    let mut r = || -10.0 + 20.0 * rng.float();
    let mut transform = Transform::IDENTITY;

    for _ in 0..10 {
        match r() as i32 {
            -10..=-4 => {
                // transform = transform * Transform::scale(Vector3f::new(r(), r(), r()).abs())
            }
            // -3..=3 => transform = transform * Transform::translation(Vector3f::new(r(), r(), r())),
            _ => {
                transform = transform * Transform::rotate(Vector3f::new(r(), r(), r()), r() * 20.0);
            }
        }
    }

    transform
}

#[test]
fn rotate() {
    let mut rng = Rng::new(3);
    let mut r = || -10.0 + 20.0 * rng.float();

    for _ in 0..200 {
        let t = Transform::rotate(Vector3f::new(r(), r(), r()), 90.0);
        let v0 = Vector3f::new(r(), r(), r());
        let v1 = v0 * t;
        let v2 = v1 * t;
        let v3 = v2 * t;
        let v4 = v3 * t;
        let vm1 = v0 * t.inverse();
        let vm2 = vm1 * t.inverse();

        for (a, b) in v0.to_array().into_iter().zip(v4.to_array()) {
            assert!((a - b).abs() <= 0.0001);
        }
        for (a, b) in v2.to_array().into_iter().zip(vm2.to_array()) {
            assert!((a - b).abs() <= 0.0001);
        }
    }
}
