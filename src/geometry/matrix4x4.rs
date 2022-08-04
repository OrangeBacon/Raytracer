use std::ops::{Index, Mul, MulAssign};

use crate::geometry::Float;

/// 4 by 4 floating point matrix
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct Matrix4x4 {
    pub data: [[Float; 4]; 4],
}

impl Default for Matrix4x4 {
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl Matrix4x4 {
    pub const IDENTITY: Self = Self {
        data: [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ],
    };

    /// Create a new matrix, [row1, row2, ..]
    pub fn new(data: &[[Float; 4]; 4]) -> Self {
        Self { data: *data }
    }

    /// Construct a matrix [row1, row2, ...]
    pub fn from_array(data: &[Float; 16]) -> Self {
        Self {
            data: [
                [data[0], data[1], data[2], data[3]],
                [data[4], data[5], data[6], data[7]],
                [data[8], data[9], data[10], data[11]],
                [data[12], data[13], data[14], data[15]],
            ],
        }
    }

    /// Construct a matrix [row1, row2, ...] from a slice of elements
    /// panics if the slice does not contain at least 16 elements
    pub fn from_slice(data: &[Float]) -> Self {
        Self::from_array(data.try_into().unwrap())
    }

    /// Get the transpose of the matrix
    #[rustfmt::skip]
    pub fn transpose(&self) -> Self {
        Self::new(&[
            [self.data[0][0], self.data[1][0], self.data[2][0], self.data[3][0]],
            [self.data[0][1], self.data[1][1], self.data[2][1], self.data[3][1]],
            [self.data[0][2], self.data[1][2], self.data[2][2], self.data[3][2]],
            [self.data[0][3], self.data[1][3], self.data[2][3], self.data[3][3]],
        ])
    }

    /// Calculate the inverse of the matrix, mat ^ -1
    /// returns None if the inverse does not exist
    /// implementation ported from
    /// https://github.com/mmp/pbrt-v3/blob/aaa552a4b9cbf9dccb71450f47b268e0ed6370e2/src/core/transform.cpp#L82
    /// (numerically stable Gauss–Jordan elimination)
    pub fn inverse(&self) -> Option<Matrix4x4> {
        let mut indxc = [0usize; 4];
        let mut indxr = [0usize; 4];
        let mut ipiv = [0usize; 4];
        let mut minv = self.data;

        for i in 0..4 {
            let mut irow = 0;
            let mut icol = 0;
            let mut big = 0.0;

            // choose pivot
            for j in 0..4 {
                if ipiv[j] != 1 {
                    for k in 0..4 {
                        if ipiv[k] == 0 {
                            if minv[j][k].abs() >= big {
                                big = minv[j][k].abs();
                                irow = j;
                                icol = k;
                            }
                        } else if ipiv[k] > 1 {
                            return None;
                        }
                    }
                }
            }

            ipiv[icol] += 1;

            // swap rows for pivot
            if irow != icol {
                for k in 0..4 {
                    let tmp = minv[irow][k];
                    minv[irow][k] = minv[icol][k];
                    minv[icol][k] = tmp;
                }
            }
            indxr[i] = irow;
            indxc[i] = icol;

            if minv[icol][icol] == 0.0 {
                return None;
            }

            // set m[icol][icol] to one by scaling row icol
            let pivinv = 1.0 / minv[icol][icol];
            minv[icol][icol] = 1.0;
            for j in 0..4 {
                minv[icol][j] *= pivinv;
            }

            // subtract this row from others to zero out their columns
            for j in 0..4 {
                if j != icol {
                    let save = minv[j][icol];
                    minv[j][icol] = 0.0;

                    for k in 0..4 {
                        minv[j][k] -= minv[icol][k] * save;
                    }
                }
            }
        }

        for j in (0..=3).rev() {
            if indxr[j] != indxc[j] {
                for k in 0..4 {
                    minv[k].swap(indxr[j], indxc[j]);
                }
            }
        }

        Some(Matrix4x4 { data: minv })
    }
}

#[cfg(test)]
#[test]
fn mat_inv() {
    let a = format!(
        "{:.2?}",
        Matrix4x4::from_array(&[
            1.0, 4.0, 5.0, -1.0, -2.0, 3.0, -1.0, 0.0, 2.0, 1.0, 1.0, 0.0, 3.0, -1.0, 2.0, 1.0
        ])
        .inverse()
        .unwrap()
    );
    let b = format!(
        "{:.2?}",
        Matrix4x4::from_array(&[
            -0.1, -0.1, 0.6, -0.1, 0.0, 0.25, 0.25, 0.0, 0.2, -0.05, -0.45, 0.2, -0.1, 0.65, -0.65,
            0.9
        ])
    );
    assert_eq!(a, b);

    let a = format!(
        "{:.2?}",
        Matrix4x4::from_array(&[
            1.0, 1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0, 1.0
        ])
        .inverse()
        .unwrap()
    );
    let b = format!(
        "{:.2?}",
        Matrix4x4::from_array(&[
            0.25, 0.25, 0.25, -0.25, 0.25, 0.25, -0.25, 0.25, 0.25, -0.25, 0.25, 0.25, -0.25, 0.25,
            0.25, 0.25,
        ])
    );
    assert_eq!(a, b)
}

impl Mul for Matrix4x4 {
    type Output = Matrix4x4;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut data: [[Float; 4]; 4] = Default::default();

        for i in 0..4 {
            for j in 0..4 {
                data[i][j] = self.data[i][0] * rhs.data[0][j]
                    + self.data[i][1] * rhs.data[1][j]
                    + self.data[i][2] * rhs.data[2][j]
                    + self.data[i][3] * rhs.data[3][j]
            }
        }

        Matrix4x4 { data }
    }
}

impl MulAssign for Matrix4x4 {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl Index<usize> for Matrix4x4 {
    type Output = [Float; 4];

    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}
