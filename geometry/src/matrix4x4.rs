use std::ops::{Index, IndexMut, Mul, MulAssign};

use crate::Float;

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
                    for (k, &ipiv) in ipiv.iter().enumerate() {
                        if ipiv == 0 {
                            if minv[j][k].abs() >= big {
                                big = minv[j][k].abs();
                                irow = j;
                                icol = k;
                            }
                        } else if ipiv > 1 {
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
                for k in &mut minv {
                    k.swap(indxr[j], indxc[j]);
                }
            }
        }

        Some(Matrix4x4 { data: minv })
    }

    /// Calculate the determinant of this 4x4 matrix
    pub fn determinant(&self) -> Float {
        fn det2(mat: [&[Float]; 2]) -> Float {
            mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]
        }
        fn det3(mat: [[Float; 3]; 3]) -> Float {
            mat[0][0] * det2([&mat[1][1..=2], &mat[2][1..=2]])
                - mat[0][1] * det2([&[mat[1][0], mat[1][2]], &[mat[2][0], mat[2][2]]])
                + mat[0][2] * det2([&mat[1][0..=1], &mat[2][0..=1]])
        }

        self[0][0] * det3(self.cofactor([1, 2, 3])) - self[0][1] * det3(self.cofactor([0, 2, 3]))
            + self[0][2] * det3(self.cofactor([0, 1, 3]))
            - self[0][3] * det3(self.cofactor([0, 1, 2]))
    }

    /// Get the 3*3 cofactor of the matrix, rows 1,2,3 (excluding 0)
    /// with the supplied columns
    fn cofactor(&self, cols: [usize; 3]) -> [[Float; 3]; 3] {
        let mut res = [[0.0; 3]; 3];

        for row in 1..4 {
            for (&col, res_col) in cols.iter().zip(0..3) {
                res[row - 1][res_col] = self[row][col]
            }
        }

        res
    }
}

impl Mul for Matrix4x4 {
    type Output = Matrix4x4;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut data: [[Float; 4]; 4] = Default::default();

        for (i, row) in self.data.iter().enumerate() {
            for j in 0..4 {
                data[i][j] = row[0] * rhs.data[0][j]
                    + row[1] * rhs.data[1][j]
                    + row[2] * rhs.data[2][j]
                    + row[3] * rhs.data[3][j]
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

impl IndexMut<usize> for Matrix4x4 {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index]
    }
}
