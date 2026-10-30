use super::*;

/// A quintic Bézier's derivative of order `d` at `t`.
fn bez(p: &[Vec3], t: f64, d: usize) -> Vec3 {
    let mut q = p.to_vec();
    for _ in 0..d {
        let n = q.len() - 1;
        q = (0..n).map(|i| (q[i + 1] - q[i]) * n as f64).collect();
    }
    let n = q.len() - 1;
    (0..=n).fold(Vec3::ZERO, |s, i| {
        let b = (1..=i).fold(1.0, |c, k| c * (n + 1 - k) as f64 / k as f64) * t.powi(i as i32) * (1.0 - t).powi((n - i) as i32);
        s + q[i] * b
    })
}

/// Patch `p`'s derivative `∂ᵘᵃ∂ᵛᵇS` at `(u, v)`, from its 6 × 6 net.
fn deriv(p: &NurbsSurface, u: f64, v: f64, a: usize, b: usize) -> Vec3 {
    let rows: Vec<Vec3> = (0..6)
        .map(|i| {
            let col: Vec<Vec3> = (0..6).map(|j| p.point(i, j).project().unwrap()).collect();
            bez(&col, v, b)
        })
        .collect();
    bez(&rows, u, a)
}

/// `[S, S_u, S_v, S_uu, S_uv, S_vv]` at `(u, v)`.
fn jet(p: &NurbsSurface, u: f64, v: f64) -> [Vec3; 6] {
    [deriv(p, u, v, 0, 0), deriv(p, u, v, 1, 0), deriv(p, u, v, 0, 1), deriv(p, u, v, 2, 0), deriv(p, u, v, 1, 1), deriv(p, u, v, 0, 2)]
}

/// The normal curvature of a patch with derivatives `s` along the tangent `dir`, about `normal`.
fn second_form(s: &[Vec3; 6], normal: Vec3, dir: Vec3) -> f64 {
    let [_, su, sv, suu, suv, svv] = *s;
    let (e, f, g) = (su.dot(su), su.dot(sv), sv.dot(sv));
    let det = e * g - f * f;
    let (du, dv) = (dir.dot(su), dir.dot(sv));
    let (a, b) = ((g * du - f * dv) / det, (e * dv - f * du) / det);
    (normal.dot(suu) * a * a + 2.0 * normal.dot(suv) * a * b + normal.dot(svv) * b * b) / (dir.len2())
}

fn unit_normal(s: &[Vec3; 6]) -> Vec3 {
    s[1].cross(s[2]).normalized().unwrap()
}

/// De Casteljau halves of a quintic.
fn halves(p: [Vec3; 6]) -> ([Vec3; 6], [Vec3; 6]) {
    let mut left = [Vec3::ZERO; 6];
    let mut right = [Vec3::ZERO; 6];
    let mut q = p.to_vec();
    for i in 0..6 {
        left[i] = q[0];
        right[5 - i] = q[q.len() - 1];
        q = q.windows(2).map(|w| (w[0] + w[1]) * 0.5).collect();
    }
    (left, right)
}

/// The hole a set-back corner of three rounds of radius `r` leaves at the corner (0, 0, 0) of
/// the box x, y, z ≤ 0, set back `s`: alternately a face's curve (on z = 0, x = 0, y = 0) between
/// two rounds' contact lines, leaving them along them with no curvature (its inner control points
/// on them), and a round's section, a quintic leaving each plane tangent to it with no curvature
/// across (the eased section, its inner control points at thirds towards the edge), the round
/// beside it its section moved along the edge. Each side split at its middle.
fn box_hole(r: f64, s: f64) -> Vec<QuadEdge> {
    let (x, y, z) = (Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0), Vec3::new(0.0, 0.0, 1.0));
    let h = (s - r) / 3.0;
    let face = |p: Vec3, tp: Vec3, q: Vec3, tq: Vec3| [p, p + tp * h, p + tp * (2.0 * h), q + tq * (2.0 * h), q + tq * h, q];
    let section = |p: Vec3, tp: Vec3, q: Vec3, tq: Vec3| [p, p + tp * (r / 3.0), p + tp * (2.0 * r / 3.0), q + tq * (2.0 * r / 3.0), q + tq * (r / 3.0), q];
    let plane = |n: Vec3| Beside::Plane { normal: n, offset: 0.0 };
    let rail = |e: Vec3| Beside::Rail { rail: e };
    let sides: [([Vec3; 6], Beside); 6] = [
        (face(Vec3::new(-s, -r, 0.0), x, Vec3::new(-r, -s, 0.0), y), plane(z)),
        (section(Vec3::new(-r, -s, 0.0), x, Vec3::new(0.0, -s, -r), z), rail(y)),
        (face(Vec3::new(0.0, -s, -r), y, Vec3::new(0.0, -r, -s), z), plane(x)),
        (section(Vec3::new(0.0, -r, -s), y, Vec3::new(-r, 0.0, -s), x), rail(z)),
        (face(Vec3::new(-r, 0.0, -s), z, Vec3::new(-s, 0.0, -r), x), plane(y)),
        (section(Vec3::new(-s, 0.0, -r), z, Vec3::new(-s, -r, 0.0), y), rail(x)),
    ];
    let knots = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
    let mut edges = Vec::new();
    for (pts, beside) in sides {
        let (a, b) = halves(pts);
        for half in [a, b] {
            edges.push(QuadEdge { degree: 5, knots: knots.clone(), points: half.to_vec(), weights: Vec::new(), beside });
        }
    }
    edges
}

const SPEC: QuadSpec = QuadSpec { tol: 1e-9, angle_tol: 1e-6, curvature_tol: 1e-6 };

/// What the fill achieves: the worst angle and curvature jump across every seam (the twelve hole
/// edges against what is beside them, the six spokes between their patches), ends included, and
/// how far any control point lies past the box's planes.
fn achieved(fill: &QuadFill, edges: &[QuadEdge]) -> (f64, f64, f64) {
    let (mut ang, mut jump) = (0.0f64, 0.0f64);
    let ts: Vec<f64> = (0..=40).map(|i| i as f64 / 40.0).collect();
    for k in 0..6 {
        let (a, b) = (&fill.patches[k], &fill.patches[(k + 1) % 6]);
        for &t in &ts {
            let (ja, jb) = (jet(a, 1.0, t), jet(b, t, 1.0));
            assert!(ja[0].dist(jb[0]) < 1e-12, "spoke {k} at {t}: the patches part");
            let (na, nb) = (unit_normal(&ja), unit_normal(&jb));
            ang = ang.max(na.cross(nb).len());
            let across = na.cross(ja[2]).normalized().unwrap();
            jump = jump.max((second_form(&ja, na, across) - second_form(&jb, na, across)).abs());
        }
    }
    for (j, e) in edges.iter().enumerate() {
        // Edge 2k is patch k's v = 0; edge 2k + 1 patch k + 1's u = 0, run back.
        let (patch, bottom) = if j % 2 == 0 { (j / 2, true) } else { ((j / 2 + 1) % 6, false) };
        for &t in &ts {
            let s = if bottom { jet(&fill.patches[patch], t, 0.0) } else { jet(&fill.patches[patch], 0.0, t) };
            let tangent = if bottom { s[1] } else { s[2] };
            let n = unit_normal(&s);
            let across = n.cross(tangent).normalized().unwrap();
            let (theirs, mine) = match e.beside {
                Beside::Plane { normal, .. } => {
                    ang = ang.max(n.cross(normal).len());
                    (0.0, second_form(&s, n, across))
                }
                Beside::Rail { rail } => {
                    // The round, the section moved along the rail: normal square to both; its
                    // normal curvature along X = α e + β E' is β² n·E'' / |X|².
                    ang = ang.max(n.dot(rail).abs());
                    // The edge's own parameter: the patch's along v = 0, run back up u = 0.
                    let te = if bottom { t } else { 1.0 - t };
                    let (e1, e2) = (bez(&e.points, te, 1), bez(&e.points, te, 2));
                    let (ee, ef, ff) = (1.0, rail.dot(e1), e1.dot(e1));
                    let det = ee * ff - ef * ef;
                    let beta = (ee * across.dot(e1) - ef * across.dot(rail)) / det;
                    (beta * beta * n.dot(e2), second_form(&s, n, across))
                }
            };
            jump = jump.max((mine - theirs).abs());
        }
    }
    let past = fill.patches.iter().flat_map(|p| (0..6).flat_map(move |i| (0..6).map(move |j| p.point(i, j).project().unwrap()))).fold(f64::NEG_INFINITY, |m, q| m.max(q.x).max(q.y).max(q.z));
    (ang, jump, past)
}

/// The patches' smallest principal curvature about the outward normal (+(1, 1, 1) at the
/// corner), over a grid on each: never negative on a convex corner.
fn least_curvature(fill: &QuadFill) -> f64 {
    let out = Vec3::new(1.0, 1.0, 1.0);
    let mut least = f64::INFINITY;
    for p in &fill.patches {
        for i in 0..=20 {
            for j in 0..=20 {
                let s = jet(p, i as f64 / 20.0, j as f64 / 20.0);
                let mut n = unit_normal(&s);
                if n.dot(out) < 0.0 {
                    n = -n;
                }
                let [_, su, sv, suu, suv, svv] = s;
                let (e, f, g) = (su.dot(su), su.dot(sv), sv.dot(sv));
                let (l, m, nn) = (-n.dot(suu), -n.dot(suv), -n.dot(svv));
                let det = e * g - f * f;
                let k = (l * nn - m * m) / det;
                let h = (e * nn - 2.0 * f * m + g * l) / (2.0 * det);
                least = least.min(h - (h * h - k).max(0.0).sqrt());
            }
        }
    }
    least
}

/// The corner of three r5 rounds set back 8.84 (the ellipsoid's default at a box corner): six
/// patches, each holding its hole edges as its boundary rows exactly, G1 and G2 across every
/// seam to rounding, ends included (the hole's corners, the sides' middles and the centre), every
/// control point inside the box and the patches convex.
#[test]
fn a_box_corners_hole_fills_g2_exactly() {
    let edges = box_hole(5.0, 8.838_834_764_831_846);
    let fill = quad_fill(&edges, &SPEC).unwrap();
    assert_eq!(fill.patches.len(), 6);
    for (j, e) in edges.iter().enumerate() {
        let (patch, bottom) = if j % 2 == 0 { (j / 2, true) } else { ((j / 2 + 1) % 6, false) };
        for i in 0..6 {
            let q = if bottom { fill.patches[patch].point(i, 0) } else { fill.patches[patch].point(0, 5 - i) };
            assert!(q.project().unwrap().dist(e.points[i]) < 1e-12, "edge {j} point {i}");
        }
    }
    let (ang, jump, past) = achieved(&fill, &edges);
    let least = least_curvature(&fill);
    eprintln!("box corner: angle {ang:e}, curvature jump {jump:e}, past the box {past:e}, least curvature {least:e}, centre {:?}", fill.centre);
    assert!(ang < 1e-12 && jump < 1e-10 && past < 1e-12 && least > -1e-9, "{ang} {jump} {past} {least}");
    // On the corner's diagonal, by its symmetry.
    let c = fill.centre;
    assert!((c.x - c.y).abs() < 1e-12 && (c.y - c.z).abs() < 1e-12 && c.x < 0.0, "{c:?}");
}

/// Set back 15 at r5, three times the radius: the hole is large against the rounds; the patches
/// still stay inside the box and convex, exactly G2 round.
#[test]
fn a_long_setback_stays_inside_the_box() {
    for s in [11.0, 15.0, 20.0] {
        let edges = box_hole(5.0, s);
        let fill = quad_fill(&edges, &SPEC).unwrap();
        let (ang, jump, past) = achieved(&fill, &edges);
        let least = least_curvature(&fill);
        eprintln!("set back {s}: angle {ang:e}, jump {jump:e}, past the box {past:e}, least curvature {least:e}");
        assert!(ang < 1e-12 && jump < 1e-10 && past < 1e-12 && least > -1e-9, "{s}: {ang} {jump} {past} {least}");
    }
}

/// Hostile and unsupported holes are refused, never filled wrongly: a hole of another size,
/// edges that part, a rational edge, two halves of a side meeting at a kink, a face curve
/// leaving off its round's rail, sides that do not alternate, non-finite points.
#[test]
fn a_hole_the_construction_does_not_take_is_refused() {
    let good = box_hole(5.0, 9.0);
    assert!(quad_fill(&good[..10], &SPEC).is_err());
    let mut apart = good.clone();
    apart[3].points[0].x += 1e-3;
    assert!(quad_fill(&apart, &SPEC).is_err());
    let mut rational = good.clone();
    rational[0].weights = vec![1.0, 2.0, 1.0, 1.0, 1.0, 1.0];
    assert!(matches!(quad_fill(&rational, &SPEC), Err(SplineError::Unsupported(_))));
    let mut kink = good.clone();
    kink[1].points[1] = kink[1].points[0] + (kink[1].points[1] - kink[1].points[0]) * 1.5;
    assert!(quad_fill(&kink, &SPEC).is_err());
    let mut off_rail = good.clone();
    if let Beside::Rail { rail } = off_rail[2].beside {
        let turned = Beside::Rail { rail: (rail + Vec3::new(0.0, 0.0, 0.01)).normalized().unwrap() };
        off_rail[2].beside = turned;
        off_rail[3].beside = turned;
    }
    assert!(matches!(quad_fill(&off_rail, &SPEC), Err(SplineError::Unsupported(_))));
    let mut same = good.clone();
    let plane = same[0].beside;
    same[2].beside = plane;
    same[3].beside = plane;
    assert!(quad_fill(&same, &SPEC).is_err());
    let mut nan = good.clone();
    nan[5].points[2].y = f64::NAN;
    assert!(quad_fill(&nan, &SPEC).is_err());
    let bad_spec = QuadSpec { tol: f64::NAN, ..SPEC };
    assert!(quad_fill(&good, &bad_spec).is_err());
}
