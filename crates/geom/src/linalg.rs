//! Dense linear algebra for small systems: a row-major [`Matrix`], its singular value
//! decomposition by one-sided Jacobi rotations (Hestenes), and the Moore–Penrose pseudo-inverse
//! built from it, which solves a system in the least-squares sense with the least-norm solution
//! (exact where the system is consistent, least squares where it is overdetermined).
//!
//! One-sided Jacobi is the textbook choice for small dense matrices: it orthogonalises the
//! columns pairwise, so the singular values come out to full relative accuracy and no
//! bidiagonalisation is needed (Golub and Van Loan, Matrix Computations, 4th ed., 8.6.3; Demmel
//! and Veselić 1992).

/// A dense matrix, its entries row by row.
#[derive(Clone, Debug, PartialEq)]
pub struct Matrix {
    rows: usize,
    cols: usize,
    data: Vec<f64>,
}

/// Most sweeps the Jacobi rotations take before they are taken as not converging: each sweep
/// squares the off-diagonal mass once the columns are nearly orthogonal, so a well-posed small
/// matrix converges in well under ten.
const MAX_SWEEPS: usize = 60;

impl Matrix {
    /// The `rows` x `cols` zero matrix.
    pub fn zeros(rows: usize, cols: usize) -> Matrix {
        Matrix { rows, cols, data: vec![0.0; rows.saturating_mul(cols)] }
    }

    /// The matrix with these rows, `None` if they are not all as long as the first or an entry
    /// is not finite.
    pub fn from_rows<R: AsRef<[f64]>>(rows: &[R]) -> Option<Matrix> {
        let cols = rows.first().map_or(0, |r| r.as_ref().len());
        let mut data = Vec::with_capacity(rows.len().saturating_mul(cols));
        for r in rows {
            let r = r.as_ref();
            if r.len() != cols || r.iter().any(|x| !x.is_finite()) {
                return None;
            }
            data.extend_from_slice(r);
        }
        Some(Matrix { rows: rows.len(), cols, data })
    }

    pub fn rows(&self) -> usize {
        self.rows
    }

    pub fn cols(&self) -> usize {
        self.cols
    }

    /// The entry at row `i`, column `j` (0 outside the matrix).
    pub fn get(&self, i: usize, j: usize) -> f64 {
        if i < self.rows && j < self.cols { self.data.get(i * self.cols + j).copied().unwrap_or(0.0) } else { 0.0 }
    }

    fn set(&mut self, i: usize, j: usize, x: f64) {
        if i < self.rows
            && j < self.cols
            && let Some(slot) = self.data.get_mut(i * self.cols + j)
        {
            *slot = x;
        }
    }

    pub fn transpose(&self) -> Matrix {
        let mut t = Matrix::zeros(self.cols, self.rows);
        for i in 0..self.rows {
            for j in 0..self.cols {
                t.set(j, i, self.get(i, j));
            }
        }
        t
    }

    /// The product `self x`, `None` if `x` is not as long as a row.
    pub fn mul_vec(&self, x: &[f64]) -> Option<Vec<f64>> {
        (x.len() == self.cols).then(|| (0..self.rows).map(|i| x.iter().enumerate().map(|(j, xj)| self.get(i, j) * xj).sum()).collect())
    }

    /// The singular values, largest first, with `U` (`rows` x `k`) and `V` (`cols` x `k`),
    /// `k = min(rows, cols)`, so that `self = U diag(sigma) V^T`; `None` if the rotations do not
    /// converge or the matrix is empty.
    pub fn svd(&self) -> Option<(Matrix, Vec<f64>, Matrix)> {
        if self.rows == 0 || self.cols == 0 {
            return None;
        }
        if self.rows < self.cols {
            let (u, s, v) = self.transpose().svd()?;
            return Some((v, s, u));
        }
        let (m, n) = (self.rows, self.cols);
        let mut u = self.clone();
        let mut v = Matrix::zeros(n, n);
        for i in 0..n {
            v.set(i, i, 1.0);
        }
        let column = |a: &Matrix, p: usize, q: usize| -> (f64, f64, f64) {
            (0..m).fold((0.0, 0.0, 0.0), |(app, aqq, apq), i| {
                let (x, y) = (a.get(i, p), a.get(i, q));
                (app + x * x, aqq + y * y, apq + x * y)
            })
        };
        let mut converged = false;
        for _ in 0..MAX_SWEEPS {
            let mut rotated = false;
            for p in 0..n {
                for q in p + 1..n {
                    let (app, aqq, apq) = column(&u, p, q);
                    if apq.abs() <= f64::EPSILON * (app * aqq).sqrt() || apq == 0.0 {
                        continue;
                    }
                    rotated = true;
                    // The rotation making columns p and q orthogonal.
                    let zeta = (aqq - app) / (2.0 * apq);
                    let t = zeta.signum() / (zeta.abs() + (1.0 + zeta * zeta).sqrt());
                    let t = if zeta == 0.0 { 1.0 } else { t };
                    let c = 1.0 / (1.0 + t * t).sqrt();
                    let s = c * t;
                    for a in [&mut u, &mut v] {
                        for i in 0..a.rows {
                            let (x, y) = (a.get(i, p), a.get(i, q));
                            a.set(i, p, c * x - s * y);
                            a.set(i, q, s * x + c * y);
                        }
                    }
                }
            }
            if !rotated {
                converged = true;
                break;
            }
        }
        if !converged {
            return None;
        }
        let mut sigma: Vec<(f64, usize)> = (0..n).map(|j| ((0..m).map(|i| u.get(i, j).powi(2)).sum::<f64>().sqrt(), j)).collect();
        sigma.sort_by(|a, b| b.0.total_cmp(&a.0));
        let (mut uo, mut vo) = (Matrix::zeros(m, n), Matrix::zeros(n, n));
        for (k, (s, j)) in sigma.iter().enumerate() {
            for i in 0..m {
                uo.set(i, k, if *s > 0.0 { u.get(i, *j) / s } else { 0.0 });
            }
            for i in 0..n {
                vo.set(i, k, v.get(i, *j));
            }
        }
        let values: Vec<f64> = sigma.into_iter().map(|(s, _)| s).collect();
        values.iter().all(|s| s.is_finite()).then_some((uo, values, vo))
    }

    /// The Moore–Penrose pseudo-inverse (`cols` x `rows`): `V diag(1 / sigma) U^T` over the
    /// singular values above `max(rows, cols) eps sigma_max` (the usual rank cut-off), so
    /// `pseudo_inverse() b` is the least-norm least-squares solution of `self x = b`. `None` if
    /// the decomposition fails.
    pub fn pseudo_inverse(&self) -> Option<Matrix> {
        let (u, sigma, v) = self.svd()?;
        let largest = sigma.first().copied().unwrap_or(0.0);
        let cut = self.rows.max(self.cols) as f64 * f64::EPSILON * largest;
        let mut p = Matrix::zeros(self.cols, self.rows);
        for (k, s) in sigma.iter().enumerate() {
            if *s <= cut || *s == 0.0 {
                continue;
            }
            for i in 0..self.cols {
                for j in 0..self.rows {
                    p.set(i, j, p.get(i, j) + v.get(i, k) * u.get(j, k) / s);
                }
            }
        }
        Some(p)
    }

    /// The rank: the singular values above the pseudo-inverse's cut-off.
    pub fn rank(&self) -> usize {
        let Some((_, sigma, _)) = self.svd() else { return 0 };
        let largest = sigma.first().copied().unwrap_or(0.0);
        let cut = self.rows.max(self.cols) as f64 * f64::EPSILON * largest;
        sigma.iter().filter(|s| **s > cut && **s > 0.0).count()
    }
}

#[cfg(test)]
#[path = "linalg_tests.rs"]
mod tests;
