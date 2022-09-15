use crate::{Matrix4x4, Float, Rng};

#[test]
fn mat_inv() {
    let mut rng = Rng::new(4); // chosen by totally fair dice roll

    for _ in 0..1000 {
        let data = [(); 16].map(|_| rng.uniform_float::<Float>());
        let mat = Matrix4x4::from_array(&data);

        if let Some(inv) = mat.inverse() {
            for (row, ident) in (mat * inv).data.into_iter().zip(Matrix4x4::<Float>::IDENTITY.data) {
                for (element, ident) in row.into_iter().zip(ident) {
                    assert!(element.abs() <= ident + 0.02)
                }
            }
        } else {
            assert!(mat.determinant().abs() <= 0.02);
        }
    }
}
