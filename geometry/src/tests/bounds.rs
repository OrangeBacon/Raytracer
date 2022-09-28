use crate::{Bounds2f, Bounds2i, Bounds3f, ConstZero, Point2f, Point2i, Point3f};

#[test]
fn bound_iter() {
    let bound = Bounds2i::new(Point2i::new(0, 1), Point2i::new(2, 3));
    let expected = [
        Point2i::new(0, 1),
        Point2i::new(1, 1),
        Point2i::new(0, 2),
        Point2i::new(1, 2),
    ];
    let mut offset = 0;
    for point in bound {
        assert_eq!(point, expected[offset]);
        offset += 1;
    }
}

#[test]
fn bound_iter_degenerate() {
    let bound = Bounds2i::new(Point2i::new(0, 0), Point2i::new(0, 10));
    for point in bound {
        panic!("point {:?} generated in unreachable loop", point);
    }

    let bound = Bounds2i::new(Point2i::new(0, 0), Point2i::new(4, 0));
    for point in bound {
        panic!("point {:?} generated in unreachable loop", point);
    }

    let bound = Bounds2i::ZERO;
    for point in bound {
        panic!("point {:?} generated in unreachable loop", point);
    }
}

#[test]
fn point_distance() {
    let b = Bounds3f::new(Point3f::ZERO, Point3f::new(1.0, 1.0, 1.0));

    // Points inside the bounding box or on faces
    assert_eq!(0.0, b.distance(Point3f::new(0.5, 0.5, 0.5)));
    assert_eq!(0.0, b.distance(Point3f::new(0.0, 1.0, 1.0)));
    assert_eq!(0.0, b.distance(Point3f::new(0.25, 0.8, 1.0)));
    assert_eq!(0.0, b.distance(Point3f::new(0.0, 0.25, 0.8)));
    assert_eq!(0.0, b.distance(Point3f::new(0.7, 0.0, 0.8)));

    // Aligned with the plane of one of the faces
    assert_eq!(5.0, b.distance(Point3f::new(6.0, 1.0, 1.0)));
    assert_eq!(10.0, b.distance(Point3f::new(0.0, -10.0, 1.0)));

    // 2 of the dimensions inside the box's extent
    assert_eq!(2.0, b.distance(Point3f::new(0.5, 0.5, 3.0)));
    assert_eq!(3.0, b.distance(Point3f::new(0.5, 0.5, -3.0)));
    assert_eq!(2.0, b.distance(Point3f::new(0.5, 3.0, 0.5)));
    assert_eq!(3.0, b.distance(Point3f::new(0.5, -3.0, 0.5)));
    assert_eq!(2.0, b.distance(Point3f::new(3.0, 0.5, 0.5)));
    assert_eq!(3.0, b.distance(Point3f::new(-3.0, 0.5, 0.5)));

    // General points
    assert_eq!(
        3.0 * 3.0 + 7.0 * 7.0 + 10.0 * 10.0,
        b.distance_squared(Point3f::new(4.0, 8.0, -10.0))
    );
    assert_eq!(
        6.0 * 6.0 + 10.0 * 10.0 + 7.0 * 7.0,
        b.distance_squared(Point3f::new(-6.0, -10.0, 8.0))
    );

    // A few with a more irregular box, just to be sure
    let b = Bounds3f::new(Point3f::new(-1.0, -3.0, 5.0), Point3f::new(2.0, -2.0, 18.0));
    assert_eq!(0.0, b.distance(Point3f::new(-0.99, -2.0, 5.0)));
    assert_eq!(
        2.0 * 2.0 + 6.0 * 6.0 + 4.0 * 4.0,
        b.distance_squared(Point3f::new(-3.0, -9.0, 22.0))
    );
}

#[test]
fn union2() {
    let a = Bounds2f::new(Point2f::new(-10.0, -10.0), Point2f::new(0.0, 20.0));
    let b = Bounds2f::ZERO;
    let c = a.union_box(b);

    assert_eq!(a, c);
    assert_eq!(b, b.union_box(b));

    let d = Bounds2f::at_point(Point2f::new(-15.0, 10.0));
    let e = a.union_box(d);
    assert_eq!(
        e,
        Bounds2f::new(Point2f::new(-15.0, -10.0), Point2f::new(0.0, 20.0))
    );
}

#[test]
fn union3() {
    let a = Bounds3f::new(
        Point3f::new(-10.0, -10.0, 5.0),
        Point3f::new(0.0, 20.0, 10.0),
    );
    let b = Bounds3f::ZERO;
    let c = a.union_box(b);

    assert_eq!(a, c);
    assert_eq!(b, b.union_box(b));

    let d = Bounds3f::at_point(Point3f::new(-15.0, 10.0, 30.0));
    let e = a.union_box(d);
    assert_eq!(
        e,
        Bounds3f::new(
            Point3f::new(-15.0, -10.0, 5.0),
            Point3f::new(0.0, 20.0, 30.0)
        )
    );
}
