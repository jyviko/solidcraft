use super::*;

fn close(a: &Matrix, b: &Matrix, tol: f64) -> bool {
    a.rows() == b.rows() && a.cols() == b.cols() && (0..a.rows()).all(|i| (0..a.cols()).all(|j| (a.get(i, j) - b.get(i, j)).abs() <= tol))
}

fn product(a: &Matrix, b: &Matrix) -> Matrix {
    let mut p = Matrix::zeros(a.rows(), b.cols());
    for i in 0..a.rows() {
        for j in 0..b.cols() {
            p.set(i, j, (0..a.cols()).map(|k| a.get(i, k) * b.get(k, j)).sum());
        }
    }
    p
}

/// The four Penrose conditions: A A+ A = A, A+ A A+ = A+, A A+ and A+ A symmetric.
fn penrose(a: &Matrix) {
    let p = a.pseudo_inverse().unwrap();
    assert_eq!((p.rows(), p.cols()), (a.cols(), a.rows()));
    assert!(close(&product(&product(a, &p), a), a, 1e-12));
    assert!(close(&product(&product(&p, a), &p), &p, 1e-12));
    let (ap, pa) = (product(a, &p), product(&p, a));
    assert!(close(&ap, &ap.transpose(), 1e-12));
    assert!(close(&pa, &pa.transpose(), 1e-12));
}

/// A square invertible matrix: the pseudo-inverse is the inverse.
#[test]
fn the_pseudo_inverse_of_an_invertible_matrix_is_its_inverse() {
    let a = Matrix::from_rows(&[[4.0, 1.0, 0.0], [1.0, 3.0, 1.0], [0.0, 1.0, 2.0]]).unwrap();
    let p = a.pseudo_inverse().unwrap();
    let mut id = Matrix::zeros(3, 3);
    for i in 0..3 {
        id.set(i, i, 1.0);
    }
    assert!(close(&product(&a, &p), &id, 1e-13));
    penrose(&a);
}

/// More rows than columns: the least-squares solution, whose residual is square to the columns.
#[test]
fn an_overdetermined_system_takes_its_least_squares_solution() {
    // A line through (0, 1), (1, 2.9), (2, 5.2), (3, 6.8).
    let a = Matrix::from_rows(&[[1.0, 0.0], [1.0, 1.0], [1.0, 2.0], [1.0, 3.0]]).unwrap();
    let b = [1.0, 2.9, 5.2, 6.8];
    let x = a.pseudo_inverse().unwrap().mul_vec(&b).unwrap();
    let r: Vec<f64> = a.mul_vec(&x).unwrap().iter().zip(b).map(|(ax, b)| ax - b).collect();
    let at_r = a.transpose().mul_vec(&r).unwrap();
    assert!(at_r.iter().all(|v| v.abs() < 1e-12), "{at_r:?}");
    assert!((x[0] - 1.02).abs() < 1e-12 && (x[1] - 1.97).abs() < 1e-12, "{x:?}");
    penrose(&a);
}

/// Fewer rows than columns, and a rank-deficient matrix: the least-norm solution.
#[test]
fn an_underdetermined_or_singular_system_takes_the_least_norm_solution() {
    let a = Matrix::from_rows(&[[1.0, 1.0, 0.0]]).unwrap();
    let x = a.pseudo_inverse().unwrap().mul_vec(&[2.0]).unwrap();
    assert!((x[0] - 1.0).abs() < 1e-14 && (x[1] - 1.0).abs() < 1e-14 && x[2].abs() < 1e-14, "{x:?}");
    let singular = Matrix::from_rows(&[[1.0, 2.0], [2.0, 4.0], [3.0, 6.0]]).unwrap();
    assert_eq!(singular.rank(), 1);
    penrose(&singular);
    penrose(&a);
    assert!(Matrix::from_rows(&[vec![1.0, f64::NAN]]).is_none());
    assert!(Matrix::from_rows(&[vec![1.0, 2.0], vec![1.0]]).is_none());
    assert!(Matrix::zeros(0, 0).pseudo_inverse().is_none());
}
