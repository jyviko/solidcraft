use solidcraft_geom::Vec3;
use solidcraft_geom::brep as b;

use crate::*;

/// The three r5 rounds at the corner (60, 40, 20) of a 60 x 40 x 20 box (ladder 07/009 and
/// 07/015).
fn box_corner() -> (Body, [Vec3; 3]) {
    let b = box_solid(Vec3::new(0.0, 0.0, 0.0), Vec3::new(60.0, 40.0, 20.0)).unwrap();
    (b, [Vec3::new(30.0, 40.0, 20.0), Vec3::new(60.0, 20.0, 20.0), Vec3::new(60.0, 40.0, 10.0)])
}

/// The corner's seams as the ladder grades them ([`corner_continuity`]): every tangent seam
/// within `within` of `near`, its curvature jump held `clear` mm from the vertices where seams
/// turn (the hole's corners, where a round's own G1 seam with a face ends; the sides' middles
/// and the centre, where the patches' spokes meet). Each seam is printed.
fn graded(name: &str, body: &Body, near: Vec3, within: f64, clear: f64) -> CornerContinuity {
    let c = corner_continuity(body, near, within, clear).unwrap();
    eprintln!(
        "{name}: {} seams, gap {:.1e} mm, angle {:.4}°, curvature jump {:.3e} per mm, {:.3e} clear of the turns ({} stations)",
        c.seams.len(),
        c.gap_mm,
        c.normal_angle_deg,
        c.curvature_jump_per_mm,
        c.curvature_jump_clear_per_mm,
        c.stations_clear
    );
    for s in &c.seams {
        eprintln!(
            "  seam {:>3} {:.2} mm: gap {:.1e}, angle {:.4}°, jump {:.3e}, clear {}",
            s.edge,
            s.length,
            s.gap_mm.value,
            s.normal_angle_deg.value,
            s.curvature_jump_per_mm.value,
            s.curvature_jump_clear_per_mm.map(|w| format!("{:.3e}", w.value)).unwrap_or_else(|| "-".into())
        );
    }
    c
}

/// The G2 corner of 07/015, set back the default 7.5 (1.5 r): one closed valid solid of the box's
/// six faces, the three rounds and the six patches about the centre (15 faces, as SMLib's own
/// setback corner, 07/013), within the box; graded as the ladder grades it, every seam of the
/// corner tangent and, 1.2 mm clear of the turns, curvature continuous to `CURVATURE_SHARE / r`
/// where the rolling ball's corner jumps 1 / r; the patches' edges are the hole's own, so the
/// gaps are rounding.
#[test]
fn the_box_corner_g2_patch_is_a_closed_solid_curvature_continuous_to_its_neighbours() {
    let (body, edges) = box_corner();
    let ball = fillet(&body, &edges, 5.0).unwrap();
    let g2 = fillet_g2_corners(&body, &edges, 5.0, &[]).unwrap();
    assert!(g2.validity().is_empty(), "{:?}", g2.validity());
    let v = crate::validity::validity(&g2);
    assert!(v.ok(), "{:?}", v);
    let m = measure(&g2).unwrap();
    eprintln!(
        "G2 corner: faces {} edges {} vertices {}, volume {:.6} (ball {:.6})",
        m.faces,
        m.edges,
        m.vertices,
        m.volume,
        measure(&ball).unwrap().volume
    );
    assert_eq!((m.faces, m.edges, m.vertices), (15, 36, 23));
    assert!((m.bbox.min - Vec3::ZERO).len() < 1e-6 && (m.bbox.max - Vec3::new(60.0, 40.0, 20.0)).len() < 1e-6, "{:?}", m.bbox);
    let corner = Vec3::new(60.0, 40.0, 20.0);
    let c = graded("G2 corner", &g2, corner, 15.0, 1.2);
    let b = graded("ball corner", &ball, corner, 15.0, 1.2);
    assert_eq!(c.seams.len(), 18);
    assert!(c.gap_mm < 1e-9 && c.normal_angle_deg <= crate::fill::PATCH_ANGLE_TOL_DEG, "{c:?}");
    assert!(c.stations_clear > 0 && c.curvature_jump_clear_per_mm <= CURVATURE_SHARE / 5.0, "{}", c.curvature_jump_clear_per_mm);
    assert!(b.curvature_jump_clear_per_mm > 0.19);
}

/// Set back 8, as SMLib's own setback test: outside the corner cube [52, 60] x [32, 40] x
/// [12, 20] the body is the sphere corner's, and inside it the convex body holds the hull of its
/// points on the cube's faces (299.0201, by a hull of them sampled) and lies within the three
/// rounds' prisms (`fillet_setback_tests`).
#[test]
fn set_back_8_the_g2_corner_lies_between_the_hull_and_the_prisms() {
    let pi = std::f64::consts::PI;
    let (body, edges) = box_corner();
    let sphere = measure(&fillet(&body, &edges, 5.0).unwrap()).unwrap().volume;
    let g2 = fillet_g2_corners(&body, &edges, 5.0, &[(Vec3::new(60.0, 40.0, 20.0), 8.0)]).unwrap();
    let v = measure(&g2).unwrap().volume;
    let outside = sphere - (512.0 - 225.0 * (1.0 - pi / 4.0) - 125.0 * (1.0 - pi / 6.0));
    let prisms = 512.0 - 225.0 * (1.0 - pi / 4.0) - 125.0 * (2f64.sqrt() - 1.0);
    eprintln!("G2 corner set back 8: volume {v:.6} in [{:.6}, {:.6}]", outside + 299.0201, outside + prisms);
    assert!(v > outside + 299.0201 && v < outside + prisms);
    let c = graded("G2 corner set back 8", &g2, Vec3::new(60.0, 40.0, 20.0), 15.0, 1.2);
    assert!(c.normal_angle_deg <= crate::fill::PATCH_ANGLE_TOL_DEG && c.curvature_jump_clear_per_mm <= CURVATURE_SHARE / 5.0, "{c:?}");
}

/// How far past the box (60, 40, 20) face `face` of `rec` reaches within its trim: the highest
/// `max(x - 60, y - 40, z - 20)` on it, bracketed `(lo, hi)` from its exact surface alone, with no
/// mesh. The surface (a B-spline, polynomial or rational with positive weights) is cut into its
/// Bézier patches and those split in four (de Casteljau on the homogeneous controls, exact up to
/// rounding) while they can still hold a higher point: a patch lies in the hull of its
/// (projected) control points, so their highest value bounds it from above
/// (`hi`), and its corners are points of the surface (`lo`, where inside the trim). A patch
/// clear of the trim (the face's pcurves, sampled densely in its (u, v)) is dropped; one whose
/// hull cannot beat `lo` by 1e-12 is done; one split 18 times counts its hull in `hi`. Returns
/// the bracket and the (u, v) where `lo` is reached.
fn past_the_box(rec: &b::Brep, face: usize) -> (f64, f64, [f64; 2]) {
    use solidcraft_geom::bspline::{HPoint, segment};
    let f = &rec.faces[face];
    let surface = f.nurbs.as_ref().unwrap_or(&f.surface);
    let (u_knots, v_knots, controls): (&Vec<f64>, &Vec<f64>, Vec<Vec<HPoint>>) = match surface {
        b::Surface::Spline { u_knots, v_knots, points } => {
            (u_knots, v_knots, points.iter().map(|r| r.iter().map(|c| HPoint::weighted(Vec3::from(*c), 1.0)).collect()).collect())
        }
        b::Surface::RationalSpline { u_knots, v_knots, homogeneous } => {
            (u_knots, v_knots, homogeneous.iter().map(|r| r.iter().map(|h| HPoint { x: h[0], y: h[1], z: h[2], w: h[3] }).collect()).collect())
        }
        _ => panic!("face {face} is no B-spline: {:?}", std::mem::discriminant(surface)),
    };
    assert!(controls.iter().flatten().all(|h| h.w > 0.0), "face {face}: a weight is not positive");
    let point = |h: &HPoint| h.project().unwrap();
    let excess = |p: Vec3| (p.x - 60.0).max(p.y - 40.0).max(p.z - 20.0);
    // The trim: every loop's pcurves as one polygon set, inside by the even-odd rule.
    let mut rings: Vec<Vec<[f64; 2]>> = Vec::new();
    for l in &f.loops {
        let b::Loop::Edges { edges } = l else { continue };
        let mut ring = Vec::new();
        for oe in edges {
            let [t0, t1] = oe.pcurve.interval().unwrap();
            let mut seg: Vec<[f64; 2]> = (0..=400).map(|k| oe.pcurve.eval(t0 + (t1 - t0) * f64::from(k) / 400.0).unwrap()).collect();
            if !oe.forward {
                seg.reverse();
            }
            ring.extend(seg);
        }
        rings.push(ring);
    }
    let inside = |w: [f64; 2]| {
        let mut odd = false;
        for r in &rings {
            for k in 0..r.len() {
                let (p, q) = (r[k], r[(k + 1) % r.len()]);
                if (p[1] > w[1]) != (q[1] > w[1]) && w[0] < p[0] + (w[1] - p[1]) * (q[0] - p[0]) / (q[1] - p[1]) {
                    odd = !odd;
                }
            }
        }
        odd
    };
    let touches = |u: [f64; 2], v: [f64; 2]| {
        rings.iter().flatten().any(|w| (u[0]..=u[1]).contains(&w[0]) && (v[0]..=v[1]).contains(&w[1]))
            || [[u[0], v[0]], [u[1], v[0]], [u[0], v[1]], [u[1], v[1]]].into_iter().any(inside)
    };
    // The Bézier patches: each knot span in u, then in v (control point [i][j], i along u).
    let (nu, nv) = (controls.len(), controls[0].len());
    let (p, q) = (u_knots.len() - nu - 1, v_knots.len() - nv - 1);
    let spans = |k: &[f64], n: usize, d: usize| {
        let mut x: Vec<f64> = k[d..=n].to_vec();
        x.dedup();
        x.windows(2).map(|w| [w[0], w[1]]).collect::<Vec<_>>()
    };
    type Patch = (Vec<Vec<HPoint>>, [f64; 2], [f64; 2], u32);
    let mut todo: Vec<Patch> = Vec::new();
    for su in spans(u_knots, nu, p) {
        // Each v column cut to the u span: (p + 1) x nv.
        let cols: Vec<Vec<HPoint>> = (0..nv)
            .map(|j| {
                let c: Vec<HPoint> = (0..nu).map(|i| controls[i][j]).collect();
                segment(p, u_knots, &c, su[0], su[1]).unwrap().1
            })
            .collect();
        for sv in spans(v_knots, nv, q) {
            let net: Vec<Vec<HPoint>> =
                (0..=p).map(|i| segment(q, v_knots, &cols.iter().map(|c| c[i]).collect::<Vec<_>>(), sv[0], sv[1]).unwrap().1).collect();
            todo.push((net, su, sv, 0));
        }
    }
    let half = |c: &[HPoint]| {
        let (mut l, mut r, mut pts) = (Vec::new(), Vec::new(), c.to_vec());
        while let (Some(&a), Some(&z)) = (pts.first(), pts.last()) {
            l.push(a);
            r.push(z);
            pts = pts.windows(2).map(|w| (w[0] + w[1]) * 0.5).collect();
        }
        r.reverse();
        (l, r)
    };
    let (mut lo, mut hi, mut at) = (f64::NEG_INFINITY, f64::NEG_INFINITY, [f64::NAN; 2]);
    while let Some((net, u, v, depth)) = todo.pop() {
        if !touches(u, v) {
            continue;
        }
        for (i, j, uv) in [(0, 0, [u[0], v[0]]), (p, 0, [u[1], v[0]]), (0, q, [u[0], v[1]]), (p, q, [u[1], v[1]])] {
            let e = excess(point(&net[i][j]));
            if e > lo && inside(uv) {
                (lo, at) = (e, uv);
            }
        }
        let hull = net.iter().flatten().map(|c| excess(point(c))).fold(f64::NEG_INFINITY, f64::max);
        if hull <= lo + 1e-12 {
            continue;
        }
        if depth == 18 {
            hi = hi.max(hull);
            continue;
        }
        // Split in u (each column of the net is a v row: split the rows across i), then in v.
        let (mut a, mut z) = (vec![Vec::new(); p + 1], vec![Vec::new(); p + 1]);
        for j in 0..=q {
            let (l, r) = half(&net.iter().map(|row| row[j]).collect::<Vec<_>>());
            for i in 0..=p {
                a[i].push(l[i]);
                z[i].push(r[i]);
            }
        }
        let (um, vm) = (0.5 * (u[0] + u[1]), 0.5 * (v[0] + v[1]));
        for (n, uu) in [(a, [u[0], um]), (z, [um, u[1]])] {
            let (l, r): (Vec<_>, Vec<_>) = n.iter().map(|row| half(row)).unzip();
            todo.push((l, uu, [v[0], vm], depth + 1));
            todo.push((r, uu, [vm, v[1]], depth + 1));
        }
    }
    (lo, hi.max(lo), at)
}

/// Long setbacks, where SMLib's own patch bulges past the faces (beyond about 2 r, by 0.25 at
/// 3 r): the G2 patches meet the faces' own planes along their curves and stay inside the box,
/// tangent across every seam and curvature continuous clear of the turns by the corner's fade.
/// The body is convex (a box with its edges and corner rounded), and a convex face lies on the
/// inner side of every tangent plane: each patch meets the box's planes tangentially along its
/// seams there, so it reaches past none of them, up to rounding (1e-9 mm). Measured on each
/// patch's exact surface within its trim ([`past_the_box`]), not on a mesh of it.
#[test]
fn a_long_setback_keeps_the_g2_patch_inside_the_box() {
    let (body, edges) = box_corner();
    let corner = Vec3::new(60.0, 40.0, 20.0);
    let mut wrong = Vec::new();
    for s in [11.0, 12.0, 13.0, 14.0, 15.0] {
        let g2 = fillet_g2_corners(&body, &edges, 5.0, &[(corner, s)]).unwrap();
        let rec = brep_record(&g2).unwrap();
        let patches: Vec<usize> =
            (0..rec.faces.len()).filter(|f| matches!(rec.faces[*f].surface, b::Surface::Spline { .. } | b::Surface::RationalSpline { .. })).collect();
        assert_eq!(patches.len(), 6, "{s}: six patches");
        let (mut lo, mut hi, mut at) = (f64::NEG_INFINITY, f64::NEG_INFINITY, [f64::NAN; 2]);
        for f in &patches {
            let (l, h, a) = past_the_box(&rec, *f);
            if l > lo {
                (lo, at) = (l, a);
            }
            hi = hi.max(h);
        }
        let c = graded(&format!("set back {s}"), &g2, corner, 1.1 * s + 1.0, 0.15 * s);
        eprintln!("set back {s}: {} patches, past the box by {lo:.4e} at (u, v) {at:?} (at most {hi:.4e})", patches.len());
        if !(hi <= 1e-9 && c.normal_angle_deg <= crate::fill::PATCH_ANGLE_TOL_DEG && c.curvature_jump_clear_per_mm <= CURVATURE_SHARE / 5.0) {
            wrong.push(format!("set back {s}: past the box by {lo:.4e}..{hi:.4e} at (u, v) {at:?}, {c:?}"));
        }
    }
    assert!(wrong.is_empty(), "{wrong:#?}");
}

/// Every edge of a box rounded: eight corners, each set back and filled, one closed solid.
#[test]
fn a_box_with_every_edge_rounded_has_eight_g2_corners() {
    let body = box_solid(Vec3::ZERO, Vec3::new(30.0, 20.0, 10.0)).unwrap();
    let mut edges = Vec::new();
    for (x, y) in [(0.0, 0.0), (30.0, 0.0), (0.0, 20.0), (30.0, 20.0)] {
        edges.push(Vec3::new(x, y, 5.0));
    }
    for (y, z) in [(0.0, 0.0), (20.0, 0.0), (0.0, 10.0), (20.0, 10.0)] {
        edges.push(Vec3::new(15.0, y, z));
    }
    for (x, z) in [(0.0, 0.0), (30.0, 0.0), (0.0, 10.0), (30.0, 10.0)] {
        edges.push(Vec3::new(x, 10.0, z));
    }
    let g2 = fillet_g2_corners(&body, &edges, 2.0, &[]).unwrap();
    assert!(g2.validity().is_empty(), "{:?}", g2.validity());
    let v = crate::validity::validity(&g2);
    let m = measure(&g2).unwrap();
    eprintln!("eight G2 corners: {} faces, volume {:.4}, valid {}", m.faces, m.volume, v.ok());
    assert!(v.ok(), "{v:?}");
    assert_eq!(m.faces, 6 + 12 + 8 * 6);
}

/// What the G2 corner cannot do is refused: a setback naming no corner, one that does not reach
/// past where the rounds meet; two rounded edges at a corner are no G2 corner, the plain fillet.
#[test]
fn what_the_g2_corner_cannot_build_is_refused() {
    let (body, edges) = box_corner();
    let off = fillet_g2_corners(&body, &edges, 5.0, &[(Vec3::new(30.0, 20.0, 10.0), 8.0)]).unwrap_err().to_string();
    assert!(off.contains("names no G2 corner"), "{off}");
    let short = fillet_g2_corners(&body, &edges, 5.0, &[(Vec3::new(60.0, 40.0, 20.0), 4.0)]);
    eprintln!("set back 4 at r5: {short:?}");
    assert!(short.is_err());
    let two = fillet_g2_corners(&body, &edges[..2], 5.0, &[]).unwrap();
    let plain = fillet(&body, &edges[..2], 5.0).unwrap();
    assert_eq!(measure(&two).unwrap().faces, measure(&plain).unwrap().faces);
}

/// The 60 x 40 x 20 box with a 20 x 20 x 10 notch cut from its top corner (60, 40, 20): its
/// inner corner (40, 20, 10) where three concave edges meet, and the notch's outer corner
/// (60, 20, 10) where the concave edge along x meets two convex ones.
fn notched() -> Body {
    let a = box_solid(Vec3::ZERO, Vec3::new(60.0, 40.0, 20.0)).unwrap();
    let b = box_solid(Vec3::new(40.0, 20.0, 10.0), Vec3::new(70.0, 50.0, 30.0)).unwrap();
    boolean(&a, &b, BoolOp::Cut).unwrap().unwrap()
}

/// The notch's inner corner, three concave r5 rounds (ladder 07/019): set back 7.5 and filled,
/// one closed valid solid with the ball corner's faces but its sphere (13), and the six
/// patches; tangent across every seam and curvature continuous clear of the turns, where the
/// ball's jump 1 / r.
#[test]
fn a_concave_corner_takes_the_g2_patch_too() {
    let body = notched();
    let edges = [Vec3::new(50.0, 20.0, 10.0), Vec3::new(40.0, 30.0, 10.0), Vec3::new(40.0, 20.0, 15.0)];
    let ball = fillet(&body, &edges, 5.0).unwrap();
    let g2 = fillet_g2_corners(&body, &edges, 5.0, &[]).unwrap();
    assert!(g2.validity().is_empty(), "{:?}", g2.validity());
    let v = crate::validity::validity(&g2);
    assert!(v.ok(), "{v:?}");
    let (m, mb) = (measure(&g2).unwrap(), measure(&ball).unwrap());
    let near = Vec3::new(40.0, 20.0, 10.0);
    let c = graded("concave G2 corner", &g2, near, 11.0, 1.2);
    let b = graded("concave ball corner", &ball, near, 11.0, 1.2);
    eprintln!("concave G2 corner: {} faces (ball {}), volume {:.4} (ball {:.4})", m.faces, mb.faces, m.volume, mb.volume);
    assert_eq!(m.faces, mb.faces - 1 + 6);
    assert!(c.normal_angle_deg <= crate::fill::PATCH_ANGLE_TOL_DEG && c.curvature_jump_clear_per_mm <= CURVATURE_SHARE / 5.0, "{c:?}");
    assert!(b.curvature_jump_clear_per_mm > 0.19);
}

/// The notch's outer corner, two convex rounds and a concave one: SMLib sets back only a corner
/// of one convexity, so the G2 corner is refused as not supported (the rolling ball builds).
#[test]
fn a_mixed_convexity_corner_is_refused() {
    let body = notched();
    let edges = [Vec3::new(50.0, 20.0, 10.0), Vec3::new(60.0, 30.0, 10.0), Vec3::new(60.0, 20.0, 15.0)];
    assert!(fillet(&body, &edges, 5.0).is_ok());
    let err = fillet_g2_corners(&body, &edges, 5.0, &[]).unwrap_err().to_string();
    eprintln!("mixed corner: {err}");
    // The refusal names the corner and why: its rounds are not all of one convexity.
    assert!(err.contains("not supported yet") && err.contains("(60.000, 20.000, 10.000)") && err.contains("all convex or all concave"), "{err}");
}

/// Curvature-continuous rounds (`fillet_g2`, equivalent radius 5) of the three edges at the box
/// corner: SMLib's own corner patch between them is only tangent to the rounds, its curvature
/// jumping across their seams (it showed as a dome on the corner), so the corner is the G2 patch
/// instead, filled between the curvature rounds' setback sections: one closed valid solid, the
/// patch curvature-continuous to every round and face round it clear of the hole's corners.
#[test]
fn curvature_rounds_meet_at_a_g2_corner() {
    let (body, edges) = box_corner();
    let smlib_corner = fillet_g2(&body, &edges, 5.0).unwrap();
    let g2 = fillet_g2_with_g2_corners(&body, &edges, 5.0, &[]).unwrap();
    assert!(g2.validity().is_empty(), "{:?}", g2.validity());
    let v = crate::validity::validity(&g2);
    assert!(v.ok(), "{:?}", v);
    let corner = Vec3::new(60.0, 40.0, 20.0);
    let c = graded("curvature rounds, G2 corner", &g2, corner, 15.0, 1.2);
    let (clear, n) = (c.curvature_jump_clear_per_mm, c.stations_clear);
    let smlib_clear = graded("curvature rounds, SMLib's corner", &smlib_corner, corner, 15.0, 1.2).curvature_jump_clear_per_mm;
    assert!(n > 0 && clear <= CURVATURE_SHARE / 5.0, "{clear} over {n} stations");
    assert!(smlib_clear > 10.0 * CURVATURE_SHARE / 5.0, "SMLib's own corner patch: {smlib_clear}");
    // Both differ only within the setback cube (side 1.5 r = 7.5) at the corner.
    let (a, b) = (measure(&g2).unwrap().volume, measure(&smlib_corner).unwrap().volume);
    assert!((a - b).abs() < 7.5f64.powi(3), "{a} vs {b}");
}
