use super::*;

/// A corner's conditions as the G2 corner sets them: on each face between edges `u` and `w`
/// (unit, away from the vertex), the point where the contact lines `d` from each edge cross,
/// `(d w_rail u + d u_rail w) / sin`, on the ellipsoid and crossing the face square to it.
fn corner(edges: &[Vec3], rail: &[f64]) -> Vec<Condition> {
    let mut c = Vec::new();
    for i in 0..edges.len() {
        let j = (i + 1) % edges.len();
        let (u, w) = (edges[i].normalized().unwrap(), edges[j].normalized().unwrap());
        let n = u.cross(w);
        let sin = n.len();
        let p = (u * rail[j] + w * rail[i]) * (1.0 / sin);
        c.push(Condition::On(p));
        c.push(Condition::Square { at: p, normal: n });
    }
    c
}

/// A box corner, three edges square to each other, rounds of radius r: six independent
/// conditions, met exactly by the sphere through the crossings `(r, r, 0)`..., intercept `sqrt 2 r`
/// on each edge. A skewed three-edge corner with unequal rails is met exactly too.
#[test]
fn a_three_edge_corner_fixes_the_ellipsoid_exactly() {
    for r in [3.0, 5.0] {
        let c = corner(&[Vec3::X, Vec3::Y, Vec3::Z], &[r, r, r]);
        let e = CentredEllipsoid::fit(&c).unwrap();
        assert!(e.residual(&c) < 1e-12, "{}", e.residual(&c));
        for (i, row) in e.matrix().iter().enumerate() {
            for (j, a) in row.iter().enumerate() {
                let want = if i == j { 1.0 / (2.0 * r * r) } else { 0.0 };
                assert!((a - want).abs() < 1e-14, "{:?}", e.matrix());
            }
        }
        for u in [Vec3::X, Vec3::Y, Vec3::Z] {
            assert!((e.intercept(u).unwrap() - 2f64.sqrt() * r).abs() < 1e-12);
        }
        // Scaled 1.25, the sphere leaves the contact line (t, r, 0) at t = r sqrt(2 1.25^2 - 1).
        let t = e.line_exit(Vec3::new(0.0, r, 0.0), Vec3::X, 1.25).unwrap();
        assert!((t - r * (2.0 * 1.5625f64 - 1.0).sqrt()).abs() < 1e-12, "{t}");
    }
    let skew = [Vec3::new(1.0, 0.1, 0.0), Vec3::new(0.2, 1.0, 0.1), Vec3::new(0.1, -0.2, 1.0)];
    let c = corner(&skew, &[2.0, 3.0, 4.0]);
    let e = CentredEllipsoid::fit(&c).unwrap();
    assert!(e.residual(&c) < 1e-12, "{}", e.residual(&c));
}

/// The fit's own weighted residuals (as rows of its system).
fn sum_of_squares(a: [[f64; 3]; 3], c: &[Condition]) -> f64 {
    let s = (c
        .iter()
        .map(|x| match x {
            Condition::On(p) | Condition::Square { at: p, .. } => p.len2(),
        })
        .sum::<f64>()
        / c.len() as f64)
        .sqrt();
    let form = |u: Vec3, w: Vec3| (0..3).map(|i| (0..3).map(|j| comp(u, i) * a[i][j] * comp(w, j)).sum::<f64>()).sum::<f64>();
    c.iter()
        .map(|x| match x {
            Condition::On(p) => (form(*p, *p) - 1.0).powi(2),
            Condition::Square { at, normal } => (form(normal.normalized().unwrap(), *at) * s).powi(2),
        })
        .sum()
}

/// A four-edge vertex (an irregular pyramid's apex): eight conditions on six unknowns, not all
/// met; the fit is the least-squares one: no change of the six unknowns lowers the sum of the
/// squared residuals.
#[test]
fn a_four_edge_corner_takes_the_least_squares_ellipsoid() {
    let edges = [Vec3::new(1.0, 0.2, -0.8), Vec3::new(-0.1, 1.0, -0.6), Vec3::new(-1.0, 0.1, -1.0), Vec3::new(0.2, -1.0, -0.7)];
    let c = corner(&edges, &[3.0, 3.0, 3.0, 3.0]);
    let e = CentredEllipsoid::fit(&c).unwrap();
    let a = e.matrix();
    let best = sum_of_squares(a, &c);
    assert!(e.residual(&c) > 1e-6, "the four-edge conditions are met exactly: {}", e.residual(&c));
    let h = 1e-6 * a[0][0].abs();
    for (i, j) in PAIRS {
        for step in [h, -h] {
            let mut b = a;
            b[i][j] += step;
            if i != j {
                b[j][i] += step;
            }
            assert!(sum_of_squares(b, &c) >= best * (1.0 - 1e-9), "({i}, {j}) {step}: {} < {best}", sum_of_squares(b, &c));
        }
    }
    for u in edges {
        assert!(e.intercept(u).is_some());
    }
}

/// Points that only a hyperboloid passes through: refused, as a matrix that is not positive
/// definite is.
#[test]
fn a_quadric_that_is_no_ellipsoid_is_refused() {
    let c = [
        Condition::On(Vec3::X),
        Condition::On(Vec3::Y),
        Condition::On(Vec3::Z),
        Condition::On(Vec3::new(0.4, 0.4, 0.0)),
        Condition::On(Vec3::new(0.4, 0.0, 0.4)),
        Condition::On(Vec3::new(0.0, 0.4, 0.4)),
    ];
    assert_eq!(CentredEllipsoid::fit(&c), Err(EllipsoidError::NotPositiveDefinite));
    assert_eq!(CentredEllipsoid::new([[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]]), Err(EllipsoidError::NotPositiveDefinite));
    assert_eq!(CentredEllipsoid::new([[1.0, 0.5, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]), Err(EllipsoidError::Degenerate));
    assert_eq!(CentredEllipsoid::fit(&[]), Err(EllipsoidError::Degenerate));
    assert_eq!(CentredEllipsoid::fit(&[Condition::On(Vec3::new(f64::NAN, 0.0, 0.0))]), Err(EllipsoidError::Degenerate));
}
