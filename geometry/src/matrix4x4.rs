use std::ops::{Index, IndexMut, Mul, MulAssign};

use crate::Number;

/// 4 by 4 floating point matrix
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Matrix4x4<T: Number> {
    pub data: [[T; 4]; 4],
}

impl<T: Number> Default for Matrix4x4<T> {
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl<T: Number> Matrix4x4<T> {
    pub const IDENTITY: Self = Self {
        data: [
            [T::ONE, T::ZERO, T::ZERO, T::ZERO],
            [T::ZERO, T::ONE, T::ZERO, T::ZERO],
            [T::ZERO, T::ZERO, T::ONE, T::ZERO],
            [T::ZERO, T::ZERO, T::ZERO, T::ONE],
        ],
    };

    /// Create a new matrix, [row1, row2, ..]
    pub fn new(data: &[[T; 4]; 4]) -> Self {
        Self { data: *data }
    }

    /// Construct a matrix [row1, row2, ...]
    pub fn from_array(data: &[T; 16]) -> Self {
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
    pub fn from_slice(data: &[T]) -> Self {
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
    pub fn inverse(&self) -> Option<Self> {
        let mut indxc = [0usize; 4];
        let mut indxr = [0usize; 4];
        let mut ipiv = [0usize; 4];
        let mut minv = self.data;

        for i in 0..4 {
            let mut irow = 0;
            let mut icol = 0;
            let mut big = T::ZERO;

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

            if minv[icol][icol] == T::ZERO {
                return None;
            }

            // set m[icol][icol] to one by scaling row icol
            let pivinv = T::ONE / minv[icol][icol];
            minv[icol][icol] = T::ONE;
            for j in 0..4 {
                minv[icol][j] *= pivinv;
            }

            // subtract this row from others to zero out their columns
            for j in 0..4 {
                if j != icol {
                    let save = minv[j][icol];
                    minv[j][icol] = T::ZERO;

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
    pub fn determinant(&self) -> T {
        fn det2<T: Number>(mat: [&[T]; 2]) -> T {
            mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]
        }
        fn det3<T: Number>(mat: [[T; 3]; 3]) -> T {
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
    fn cofactor(&self, cols: [usize; 3]) -> [[T; 3]; 3] {
        let mut res = [[T::ZERO; 3]; 3];

        for row in 1..4 {
            for (&col, res_col) in cols.iter().zip(0..3) {
                res[row - 1][res_col] = self[row][col]
            }
        }

        res
    }
}

impl<T: Number> Mul for Matrix4x4<T> {
    type Output = Matrix4x4<T>;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut data: [[T; 4]; 4] = Default::default();

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

impl<T: Number> MulAssign for Matrix4x4<T> {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

impl<T: Number> Index<usize> for Matrix4x4<T> {
    type Output = [T; 4];

    fn index(&self, index: usize) -> &Self::Output {
        &self.data[index]
    }
}

impl<T: Number> IndexMut<usize> for Matrix4x4<T> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[index]
    }
}
