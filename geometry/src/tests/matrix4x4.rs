use crate::{Float, Matrix4x4, Rng};

#[test]
fn mat_inv() {
    let mut rng = Rng::new(4); // chosen by totally fair dice roll

    for i in 0..1000 {
        let data = [(); 16].map(|_| rng.uniform_float::<Float>());
        let mat = Matrix4x4::from_array(&data);

        if let Some(inv) = mat.inverse() {
            let result = mat * inv;
            for (row, ident) in result
                .data
                .into_iter()
                .zip(Matrix4x4::<Float>::IDENTITY.data)
            {
                for (element, ident) in row.into_iter().zip(ident) {
                    assert!((element - ident).abs() <= 5e-3, "{}\n{:?}", i, result);
                }
            }
        } else {
            assert!(mat.determinant().abs() <= 1e-4);
        }
    }
}
