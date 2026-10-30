//! An ellipsoid centred at the origin, `x^T A x = 1` with `A` symmetric positive definite (six
//! unknowns, its axes in any direction), fitted to conditions by the Moore–Penrose
//! pseudo-inverse ([`crate::linalg::Matrix::pseudo_inverse`]): exact where the conditions fix
//! it, least squares where there are more than six independent ones.
//!
//! A G2 fillet corner's setback comes from it (`solidcraft_kernel`'s `corner_g2`): through the
//! points where neighbouring rounds' contact lines cross on each face, crossing each face square
//! to it there.

use crate::Vec3;
use crate::linalg::Matrix;

/// One condition on the ellipsoid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Condition {
    /// The ellipsoid passes through the point: `p^T A p = 1`.
    On(Vec3),
    /// At `at`, the ellipsoid's normal `A at` is square to `normal`: it crosses the plane through
    /// the centre with that normal at right angles there, as an ellipsoid whose axes lie in or
    /// across the plane does everywhere.
    Square { at: Vec3, normal: Vec3 },
}

/// Why no ellipsoid fits.
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
pub enum EllipsoidError {
    /// The conditions are not finite, or name no point away from the centre.
    #[error("the ellipsoid's conditions are not finite or name no point away from its centre")]
    Degenerate,
    /// The fitted quadric is not an ellipsoid (`A` has an eigenvalue at or below zero: a
    /// hyperboloid or a cylinder).
    #[error("the fitted quadric is not an ellipsoid (its matrix is not positive definite)")]
    NotPositiveDefinite,
}

/// `x^T A x = 1`, `A` symmetric positive definite.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CentredEllipsoid {
    a: [[f64; 3]; 3],
}

/// The unknowns' order in a condition's row: `a11, a22, a33, a12, a13, a23`.
const PAIRS: [(usize, usize); 6] = [(0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2)];

fn comp(v: Vec3, i: usize) -> f64 {
    match i {
        0 => v.x,
        1 => v.y,
        _ => v.z,
    }
}

/// The row of `u^T A w` in the unknowns (both orders of an off-diagonal pair; for `u = w`,
/// `p^T A p` with its off-diagonal terms twice).
fn bilinear_row(u: Vec3, w: Vec3) -> [f64; 6] {
    PAIRS.map(|(i, j)| if i == j { comp(u, i) * comp(w, i) } else { comp(u, i) * comp(w, j) + comp(u, j) * comp(w, i) })
}

impl CentredEllipsoid {
    /// The ellipsoid with matrix `a`, refused unless symmetric, finite and positive definite
    /// (every pivot of its Cholesky factorisation above `1e-12` of its largest diagonal entry).
    pub fn new(a: [[f64; 3]; 3]) -> Result<CentredEllipsoid, EllipsoidError> {
        if a.iter().flatten().any(|x| !x.is_finite()) {
            return Err(EllipsoidError::Degenerate);
        }
        let scale = (0..3).map(|i| a[i][i].abs()).fold(0.0f64, f64::max);
        let symmetric = (0..3).all(|i| (0..3).all(|j| (a[i][j] - a[j][i]).abs() <= 1e-12 * scale.max(f64::MIN_POSITIVE)));
        if !symmetric {
            return Err(EllipsoidError::Degenerate);
        }
        // Cholesky: A = L L^T exists with positive pivots iff A is positive definite.
        let mut l = [[0.0f64; 3]; 3];
        for j in 0..3 {
            let pivot = a[j][j] - (0..j).map(|k| l[j][k] * l[j][k]).sum::<f64>();
            if !(pivot > 1e-12 * scale) {
                return Err(EllipsoidError::NotPositiveDefinite);
            }
            l[j][j] = pivot.sqrt();
            for i in j + 1..3 {
                l[i][j] = (a[i][j] - (0..j).map(|k| l[i][k] * l[j][k]).sum::<f64>()) / l[j][j];
            }
        }
        Ok(CentredEllipsoid { a })
    }

    /// The ellipsoid the conditions fix, by the pseudo-inverse of their system: exact where they
    /// are six independent ones, the least-squares fit (each row made dimensionless by the
    /// points' root-mean-square distance from the centre) where there are more, the least-norm
    /// one of those that meet them where there are fewer. Refused where the fitted quadric is
    /// not an ellipsoid.
    pub fn fit(conditions: &[Condition]) -> Result<CentredEllipsoid, EllipsoidError> {
        let points: Vec<Vec3> = conditions
            .iter()
            .map(|c| match c {
                Condition::On(p) | Condition::Square { at: p, .. } => *p,
            })
            .collect();
        if points.is_empty() || points.iter().any(|p| !p.is_finite()) {
            return Err(EllipsoidError::Degenerate);
        }
        let s = (points.iter().map(|p| p.len2()).sum::<f64>() / points.len() as f64).sqrt();
        if !(s > 0.0 && s.is_finite()) {
            return Err(EllipsoidError::Degenerate);
        }
        // Unknowns B = A s^2 over points q = p / s: a point's row is q^T B q = 1, a square
        // condition's n^T B q = 0 (n unit), both dimensionless.
        let mut rows: Vec<[f64; 6]> = Vec::with_capacity(conditions.len());
        let mut rhs: Vec<f64> = Vec::with_capacity(conditions.len());
        for c in conditions {
            match c {
                Condition::On(p) => {
                    let q = *p * (1.0 / s);
                    rows.push(bilinear_row(q, q));
                    rhs.push(1.0);
                }
                Condition::Square { at, normal } => {
                    let n = normal.normalized().ok_or(EllipsoidError::Degenerate)?;
                    rows.push(bilinear_row(n, *at * (1.0 / s)));
                    rhs.push(0.0);
                }
            }
        }
        let m = Matrix::from_rows(&rows).ok_or(EllipsoidError::Degenerate)?;
        let b = m.pseudo_inverse().and_then(|p| p.mul_vec(&rhs)).ok_or(EllipsoidError::Degenerate)?;
        let k = 1.0 / (s * s);
        let x = |i: usize| b.get(i).copied().unwrap_or(0.0) * k;
        CentredEllipsoid::new([[x(0), x(3), x(4)], [x(3), x(1), x(5)], [x(4), x(5), x(2)]])
    }

    /// The matrix `A`.
    pub fn matrix(&self) -> [[f64; 3]; 3] {
        self.a
    }

    /// `u^T A w`.
    pub fn form(&self, u: Vec3, w: Vec3) -> f64 {
        (0..3).map(|i| (0..3).map(|j| comp(u, i) * self.a[i][j] * comp(w, j)).sum::<f64>()).sum()
    }

    /// The distance from the centre to the ellipsoid along `dir`, `None` for a zero direction.
    pub fn intercept(&self, dir: Vec3) -> Option<f64> {
        let u = dir.normalized()?;
        let q = self.form(u, u);
        (q > 0.0).then(|| 1.0 / q.sqrt()).filter(|t| t.is_finite())
    }

    /// Where the line `origin + t dir` (`dir` unit) leaves the ellipsoid scaled by `scale` about
    /// its centre (`x^T A x = scale^2`): the larger `t`, `None` where the line misses it.
    pub fn line_exit(&self, origin: Vec3, dir: Vec3, scale: f64) -> Option<f64> {
        let u = dir.normalized()?;
        let (qa, qb, qc) = (self.form(u, u), self.form(origin, u), self.form(origin, origin) - scale * scale);
        let disc = qb * qb - qa * qc;
        (qa > 0.0 && disc >= 0.0).then(|| (-qb + disc.sqrt()) / qa).filter(|t| t.is_finite())
    }

    /// The largest absolute residual of the conditions, each as [`CentredEllipsoid::fit`] weighs it.
    pub fn residual(&self, conditions: &[Condition]) -> f64 {
        let points: Vec<Vec3> = conditions
            .iter()
            .map(|c| match c {
                Condition::On(p) | Condition::Square { at: p, .. } => *p,
            })
            .collect();
        let s = (points.iter().map(|p| p.len2()).sum::<f64>() / points.len().max(1) as f64).sqrt();
        conditions
            .iter()
            .map(|c| match c {
                Condition::On(p) => (self.form(*p, *p) - 1.0).abs(),
                Condition::Square { at, normal } => normal.normalized().map_or(f64::INFINITY, |n| (self.form(n, *at) * s).abs()),
            })
            .fold(0.0, f64::max)
    }
}

#[cfg(test)]
#[path = "ellipsoid_tests.rs"]
mod tests;
