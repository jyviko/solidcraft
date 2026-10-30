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

/// The corner's patches: the faces round their centre, the one vertex six faces meet at.
fn patches_of(rec: &b::Brep) -> Vec<usize> {
    let mut at_vertex: Vec<Vec<usize>> = vec![Vec::new(); rec.vertices.len()];
    for e in &rec.edges {
        for v in e.vertices {
            at_vertex[v as usize].extend(e.uses.iter().map(|u| u.face as usize));
        }
    }
    at_vertex
        .iter()
        .map(|f| {
            let mut f = f.clone();
            f.sort_unstable();
            f.dedup();
            f
        })
        .find(|f| f.len() == 6)
        .expect("the patches' centre")
}

/// The worst curvature jump across the seams between the corner's patches (no plane or round of
/// the box) and the rounds (cylinders, and the eased stretches before the setbacks), over each
/// seam's whole length, its ends included: the exact surfaces' normal sections square to the
/// seam ([`seam_diagnostics`]), no mesh.
fn round_seams(body: &Body) -> f64 {
    let rec = brep_record(body).unwrap();
    let plane = |f: usize| matches!(rec.faces[f].surface, b::Surface::Plane { .. } | b::Surface::PlaneThrough { .. });
    let centre = patches_of(&rec);
    let patch = |f: usize| centre.contains(&f);
    let edges: Vec<usize> = (0..rec.edges.len())
        .filter(|e| {
            let faces: Vec<usize> = rec.edges[*e].uses.iter().map(|u| u.face as usize).collect();
            faces.iter().any(|f| patch(*f)) && faces.iter().any(|f| !patch(*f) && !plane(*f))
        })
        .collect();
    assert!(!edges.is_empty(), "no seam between the corner and a round");
    let report =
        seam_diagnostics(body, &edges, 65, ContinuityTolerances { gap_mm: 1e-3, angle_deg: 1.0, curvature_per_mm: 1e-6, ..Default::default() })
            .unwrap();
    report.seams.iter().flat_map(|s| s.stations.iter()).map(|st| st.curvature_jump_per_mm).fold(0.0, f64::max)
}

/// How far the explicit patches may miss G2 anywhere on a seam of the corner, its ends included
/// (the hole's corners, the sides' middles, the centre): rounding, and the hole data's own (SMLib
/// splits a section at its parameter's middle to 1e-9): an angle in degrees (the grading reads it
/// as the arc cosine of the normals' dot, which cannot resolve less than about 1.5e-6°) and a
/// curvature jump per mm.
const EXACT_ANGLE_DEG: f64 = 1e-5;
const EXACT_JUMP_PER_MM: f64 = 1e-9;

/// The G2 corner of 07/015 at its default setback (1.25 x the ellipsoid's intercept sqrt 2 r,
/// 8.84 at r5): one closed valid solid of the box's six faces, the three rounds and the six
/// patches about the centre (15 faces, as SMLib's own setback corner, 07/013), within the box;
/// graded as the ladder grades it, every seam of the corner G1 and G2 along its whole length,
/// its ends included, where the rolling ball's corner jumps 1 / r; the patches' edges are the
/// hole's own, so the gaps are rounding.
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
    assert!(c.gap_mm < 1e-9 && c.normal_angle_deg <= EXACT_ANGLE_DEG && c.curvature_jump_per_mm <= EXACT_JUMP_PER_MM, "{c:?}");
    assert!(b.curvature_jump_clear_per_mm > 0.19);
    // G2 to each round along the whole seam, to its ends on the faces the rounds end on; and on
    // a 20 mm box with r3 rounds the same.
    let rounds = round_seams(&g2);
    eprintln!("G2 corner: worst jump to the rounds along their whole seams {rounds:.4e} per mm");
    assert!(rounds <= EXACT_JUMP_PER_MM, "{rounds}");
    let small = box_solid(Vec3::ZERO, Vec3::new(20.0, 20.0, 20.0)).unwrap();
    let edges = [Vec3::new(10.0, 20.0, 20.0), Vec3::new(20.0, 10.0, 20.0), Vec3::new(20.0, 20.0, 10.0)];
    let g2 = fillet_g2_corners(&small, &edges, 3.0, &[]).unwrap();
    let rounds = round_seams(&g2);
    eprintln!("G2 corner, 20 mm box at r3: worst jump to the rounds along their whole seams {rounds:.4e} per mm");
    assert!(rounds <= EXACT_JUMP_PER_MM, "{rounds}");
}

/// Set back 8, as SMLib's own setback test. The rounds keep their circular sections up to where
/// their ease starts, `a = 8 + EASE_FACTOR r` from the vertex (12 at r5), so outside the cube
/// `[60 - a, 60] x [40 - a, 40] x [20 - a, 20]` the body is the ball corner's exactly: each cut
/// by that cube takes the same volume, to what the two cuts and their measures hold (SMLib's
/// exact volumes, to 1e-8 of each; and each cut's faces on the cube's three faces, at most 3 a²,
/// within its edges' tolerance of the exact planes). Inside the cube the eased rounds are not convex along their edges (the
/// curvature-continuous section's middle lies 0.240 r from the edge, the circle's 0.293 r, and
/// the ease moves it by a smoothstep, not a concave law), so no tangent plane bounds the body
/// there. What holds is what the construction removes: the rounds take material only within the
/// r x r strip along each edge, and the patches only inside the setback cube of side 8, so
/// inside the cube the body keeps at least `a^3 - 8^3 - 3 (a - 8) r^2`; and it lies in the box.
#[test]
fn set_back_8_the_g2_corner_is_the_ball_corner_outside_its_eased_cube() {
    let (body, edges) = box_corner();
    let (r, setback) = (5.0, 8.0);
    let a = setback + crate::corner_g2::EASE_FACTOR * r;
    let corner = Vec3::new(60.0, 40.0, 20.0);
    let ball = fillet(&body, &edges, r).unwrap();
    let g2 = fillet_g2_corners(&body, &edges, r, &[(corner, setback)]).unwrap();
    let cube = box_solid(corner - Vec3::new(a, a, a), corner).unwrap();
    let cut = |b: &Body| boolean(b, &cube, BoolOp::Cut).unwrap().unwrap();
    let (g2_cut, ball_cut) = (cut(&g2), cut(&ball));
    let (g2_out, ball_out) = (measure(&g2_cut).unwrap().volume, measure(&ball_cut).unwrap().volume);
    let edge_tol = |b: &Body| brep_record(b).unwrap().edges.iter().map(|e| e.tolerance).fold(0.0, f64::max);
    let held = 1e-8 * (g2_out + ball_out) + 3.0 * a * a * (edge_tol(&g2_cut) + edge_tol(&ball_cut));
    let m = measure(&g2).unwrap();
    let inside = m.volume - g2_out;
    let least = a.powi(3) - setback.powi(3) - 3.0 * (a - setback) * r * r;
    eprintln!(
        "G2 corner set back 8: outside the {a} mm cube {g2_out:.9} (ball corner {ball_out:.9}, to {held:.3e}), inside {inside:.6} (at least {least:.6})"
    );
    assert!((g2_out - ball_out).abs() <= held, "{g2_out} vs {ball_out} (to {held})");
    assert!(inside >= least && inside <= a.powi(3), "{inside} not in [{least}, {}]", a.powi(3));
    assert!((m.bbox.min - Vec3::ZERO).len() < 1e-6 && (m.bbox.max - corner).len() < 1e-6, "{:?}", m.bbox);
    let c = graded("G2 corner set back 8", &g2, corner, 15.0, 1.2);
    assert!(c.normal_angle_deg <= EXACT_ANGLE_DEG && c.curvature_jump_per_mm <= EXACT_JUMP_PER_MM, "{c:?}");
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
/// G1 and G2 across every seam along its whole length.
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
        let patches = patches_of(&rec);
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
        let rounds = round_seams(&g2);
        if !(hi <= 1e-9 && rounds <= EXACT_JUMP_PER_MM && c.normal_angle_deg <= EXACT_ANGLE_DEG && c.curvature_jump_per_mm <= EXACT_JUMP_PER_MM)
        {
            wrong.push(format!("set back {s}: past the box by {lo:.4e}..{hi:.4e} at (u, v) {at:?}, to the rounds {rounds:.3e}, {c:?}"));
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

/// The default setback of a G2 corner, by its ellipsoid through the points where its rounds'
/// contact lines cross on each face (crossing each face square to it there): at a box corner the
/// sphere through (r, r, 0), (r, 0, r), (0, r, r), its intercept sqrt 2 r on each edge the least
/// setback, the setback 1.25 times that, and the rounds' ease before it 0.8 r. The same at the 60 x 40
/// x 20 box's r5 corner, the 20 mm box's r3 corner, the 40 x 30 x 20 box's r3 corner and each of
/// the eight r2 corners of a 30 x 20 x 10 box. A given setback below the intercept is refused
/// naming it.
#[test]
fn the_corner_ellipsoid_sets_the_g2_setback_and_its_ease() {
    let k = crate::corner_g2::SETBACK_SCALE;
    let check = |name: &str, body: &Body, corner: Vec3, r: f64| {
        let s = fillet_g2_corner_setback(body, corner, r).unwrap().unwrap();
        eprintln!("{name}: least {:.4}, setback {:.4}, ease {:.4} (r {r})", s.least, s.setback, s.ease);
        let least = 2f64.sqrt() * r;
        assert!((s.least - least).abs() < 1e-9 * r, "{name}: {s:?}");
        assert!((s.setback - k * least).abs() < 1e-9 * r, "{name}: {s:?}");
        assert!((s.ease - crate::corner_g2::EASE_FACTOR * r).abs() < 1e-12 * r, "{name}: {s:?}");
    };
    let (body, _) = box_corner();
    check("60 x 40 x 20 box, r5", &body, Vec3::new(60.0, 40.0, 20.0), 5.0);
    let small = box_solid(Vec3::ZERO, Vec3::new(20.0, 20.0, 20.0)).unwrap();
    check("20 mm box, r3", &small, Vec3::new(20.0, 20.0, 20.0), 3.0);
    let rig = box_solid(Vec3::ZERO, Vec3::new(40.0, 30.0, 20.0)).unwrap();
    check("40 x 30 x 20 box, r3", &rig, Vec3::new(40.0, 30.0, 20.0), 3.0);
    let eight = box_solid(Vec3::ZERO, Vec3::new(30.0, 20.0, 10.0)).unwrap();
    for x in [0.0, 30.0] {
        for y in [0.0, 20.0] {
            for z in [0.0, 10.0] {
                check(&format!("30 x 20 x 10 box corner ({x}, {y}, {z}), r2"), &eight, Vec3::new(x, y, z), 2.0);
            }
        }
    }
    let edges = [Vec3::new(20.0, 30.0, 20.0), Vec3::new(40.0, 15.0, 20.0), Vec3::new(40.0, 30.0, 10.0)];
    let short = fillet_g2_with_g2_corners(&rig, &edges, 3.0, &[(Vec3::new(40.0, 30.0, 20.0), 3.5)]).unwrap_err().to_string();
    eprintln!("set back 3.5 at r3: {short}");
    assert!(short.contains("the setback must be at least 4.2426 mm at the corner (40.000, 30.000, 20.000)"), "{short}");
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

/// The 60 x 40 x 20 box with a 20 x 20 x `depth` notch cut from its top corner (60, 40, 20): its
/// inner corner (40, 20, 20 - depth) where three concave edges meet, and the notch's outer
/// corner (60, 20, 20 - depth) where the concave edge along x meets two convex ones.
fn notched(depth: f64) -> Body {
    let a = box_solid(Vec3::ZERO, Vec3::new(60.0, 40.0, 20.0)).unwrap();
    let b = box_solid(Vec3::new(40.0, 20.0, 20.0 - depth), Vec3::new(70.0, 50.0, 30.0)).unwrap();
    boolean(&a, &b, BoolOp::Cut).unwrap().unwrap()
}

/// The inner corner of a notch 15 deep, three concave r5 rounds (ladder 07/019), at its default
/// setback (8.84, with the ease 4 just under the vertical edge's 15): one closed valid solid with
/// the ball corner's faces but its sphere (12) and the six patches; G1 and G2 across every seam
/// along its whole length, where the ball jumps 1 / r. A notch 12 deep leaves its vertical edge
/// too short for 8.84 and the ease, so the scale shrinks to where they fit, (12 - 4) / sqrt 2 r =
/// 1.1314; one 11 deep is too short even at the bare intercept, 7.07 and the ease: refused,
/// naming the length that works.
#[test]
fn a_concave_corner_takes_the_g2_patch_too() {
    let body = notched(15.0);
    let edges = [Vec3::new(50.0, 20.0, 5.0), Vec3::new(40.0, 30.0, 5.0), Vec3::new(40.0, 20.0, 12.0)];
    let near = Vec3::new(40.0, 20.0, 5.0);
    let s = fillet_g2_corner_setback(&body, near, 5.0).unwrap().unwrap();
    assert!((s.scale - crate::corner_g2::SETBACK_SCALE).abs() < 1e-12, "{s:?}");
    let ball = fillet(&body, &edges, 5.0).unwrap();
    let g2 = fillet_g2_corners(&body, &edges, 5.0, &[]).unwrap();
    assert!(g2.validity().is_empty(), "{:?}", g2.validity());
    let v = crate::validity::validity(&g2);
    assert!(v.ok(), "{v:?}");
    let (m, mb) = (measure(&g2).unwrap(), measure(&ball).unwrap());
    let c = graded("concave G2 corner", &g2, near, 11.0, 1.2);
    let b = graded("concave ball corner", &ball, near, 11.0, 1.2);
    eprintln!("concave G2 corner: {} faces (ball {}), volume {:.4} (ball {:.4})", m.faces, mb.faces, m.volume, mb.volume);
    assert_eq!(m.faces, mb.faces - 1 + 6);
    assert!(c.normal_angle_deg <= EXACT_ANGLE_DEG && c.curvature_jump_per_mm <= EXACT_JUMP_PER_MM, "{c:?}");
    assert!(b.curvature_jump_clear_per_mm > 0.19);
    let body = notched(12.0);
    let edges = [Vec3::new(50.0, 20.0, 8.0), Vec3::new(40.0, 30.0, 8.0), Vec3::new(40.0, 20.0, 12.0)];
    let s = fillet_g2_corner_setback(&body, Vec3::new(40.0, 20.0, 8.0), 5.0).unwrap().unwrap();
    eprintln!("12 mm notch: scale {:.6}, setback {:.6}, ease {:.6} (together {:.9})", s.scale, s.setback, s.ease, s.setback + s.ease);
    assert!((s.scale - 8.0 / (2f64.sqrt() * 5.0)).abs() < 1e-9 && s.setback + s.ease < 12.0, "{s:?}");
    let g2 = fillet_g2_corners(&body, &edges, 5.0, &[]).unwrap();
    assert!(crate::validity::validity(&g2).ok());
    let body = notched(11.0);
    let edges = [Vec3::new(50.0, 20.0, 9.0), Vec3::new(40.0, 30.0, 9.0), Vec3::new(40.0, 20.0, 15.0)];
    let err = fillet_g2_corners(&body, &edges, 5.0, &[]).unwrap_err().to_string();
    eprintln!("11 mm notch: {err}");
    assert!(err.contains("no room for its setback") && err.contains("an edge longer than 11.0711 mm works"), "{err}");
}

/// The notch's outer corner, two convex rounds and a concave one: SMLib sets back only a corner
/// of one convexity, so the G2 corner is refused as not supported (the rolling ball builds). A
/// notch 15 deep, set back 7.5 (with the ease, 11.5, it fits the 15 mm edge).
#[test]
fn a_mixed_convexity_corner_is_refused() {
    let body = notched(15.0);
    let edges = [Vec3::new(50.0, 20.0, 5.0), Vec3::new(60.0, 30.0, 5.0), Vec3::new(60.0, 20.0, 12.0)];
    assert!(fillet(&body, &edges, 5.0).is_ok());
    let err = fillet_g2_corners(&body, &edges, 5.0, &[(Vec3::new(60.0, 20.0, 5.0), 7.5)]).unwrap_err().to_string();
    eprintln!("mixed corner: {err}");
    // The refusal names the corner and why: its rounds are not all of one convexity.
    assert!(err.contains("not supported yet") && err.contains("(60.000, 20.000, 5.000)") && err.contains("all convex or all concave"), "{err}");
}

/// Curvature-continuous rounds (`fillet_g2`, equivalent radius 5) of the three edges at the box
/// corner: SMLib's own corner patch between them is only tangent to the rounds, its curvature
/// jumping across their seams (it showed as a dome on the corner), so the corner is the G2 patch
/// instead, filled between the curvature rounds' setback sections: one closed valid solid, the
/// patches curvature-continuous to every round and face round them along the whole of each seam.
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
    assert!(n > 0 && clear <= EXACT_JUMP_PER_MM && c.curvature_jump_per_mm <= EXACT_JUMP_PER_MM && c.normal_angle_deg <= EXACT_ANGLE_DEG, "{c:?}");
    assert!(smlib_clear > 1e-3, "SMLib's own corner patch: {smlib_clear}");
    // Both differ only within the cube of the setback and its ease at the corner.
    let s = fillet_g2_corner_setback(&body, corner, 5.0).unwrap().unwrap();
    let (a, b) = (measure(&g2).unwrap().volume, measure(&smlib_corner).unwrap().volume);
    assert!((a - b).abs() < (s.setback + s.ease).powi(3), "{a} vs {b}");
}


/// The rounds' ease into their curvature-continuous sections before the setback turns their
/// normal along the edge, most at the section's middle (it moves towards the edge by 0.075 r
/// there): measured on the exact surface of the round along the box's y edge, along the
/// bisector of its section, at setback 15 and at the default, the normal turns less than
/// `atan(0.141 / EASE_FACTOR)` (10.0° at 0.8; the quarter-radius ease turned it 29.5°, the dark
/// streak where each round met the corner).
#[test]
fn the_rounds_ease_turns_their_normal_by_about_ten_degrees() {
    use solidcraft_geom::brep::eval::SurfaceEval;
    let (body, edges) = box_corner();
    let corner = Vec3::new(60.0, 40.0, 20.0);
    let r = 5.0;
    let bound = (0.1414 / crate::corner_g2::EASE_FACTOR).atan().to_degrees() + 0.05;
    for given in [Some(15.0), None] {
        let g2 = fillet_g2_corners(&body, &edges, r, &given.map(|d| vec![(corner, d)]).unwrap_or_default()).unwrap();
        let setback = given.unwrap_or_else(|| fillet_g2_corner_setback(&body, corner, r).unwrap().unwrap().setback);
        let ease = crate::corner_g2::ease(r);
        let rec = brep_record(&g2).unwrap();
        // The round along y (x = 60, z = 20) at the middle of the bisector, clear of the ease.
        let probe = Vec3::new(60.0 - 1.46, 40.0 - setback - ease - 2.0, 20.0 - 1.46);
        let (face, s) = rec
            .faces
            .iter()
            .filter(|f| !matches!(f.surface, b::Surface::Plane { .. } | b::Surface::PlaneThrough { .. }))
            .filter_map(|f| SurfaceEval::new(f.parameters()).ok())
            .enumerate()
            .find(|(_, s)| s.project(probe, None, 1.0).and_then(|uv| s.point(uv)).is_some_and(|q| q.dist(probe) < 0.2))
            .unwrap();
        let mut worst: f64 = 0.0;
        let bisector = Vec3::new(1.0, 0.0, 1.0).normalized().unwrap();
        for k in 0..=400 {
            let y = 40.0 - setback - ease - 0.5 + (ease + 0.5) * k as f64 / 400.0;
            // The section's middle: the surface point on the bisector through (60, y, 20).
            let uv = s.project(Vec3::new(60.0, y, 20.0) - bisector * 1.6, None, 2.0).unwrap();
            let (_, su, sv) = s.derivs(uv[0], uv[1]).unwrap();
            let mut n = su.cross(sv).normalized().unwrap();
            if n.dot(bisector) < 0.0 {
                n = -n;
            }
            worst = worst.max(n.dot(bisector).clamp(-1.0, 1.0).acos().to_degrees());
        }
        eprintln!("ease {ease} at setback {setback:.4} (face {face}): the normal turns at most {worst:.3}° (bound {bound:.3}°)");
        assert!(worst < bound && worst > 5.0, "{worst}");
    }
}
