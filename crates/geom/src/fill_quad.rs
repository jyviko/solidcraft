//! Filling the hole of a set-back corner where three rounds meet with six biquintic patches
//! about a centre point, each meeting the faces beside the hole and its two neighbouring patches
//! with G2 exactly: built, not fitted.
//!
//! The hole has six sides, alternately a round's section (a straight edge's round, which near its
//! section is the section moved along the edge, its *rail*) and a planar face's curve between
//! two rounds' contact lines, each split at its middle `M_k` into two edges that are single
//! quintic spans (a set-back corner's hole from SMLib is so: the eased setback, whose sections
//! meet the planes with no curvature across and whose face curves meet the contact lines
//! curvature-continuously). Patch `k` has the corners `V_k` (where side `k` starts), `M_k`, the
//! centre `C` and `M_{k−1}`: on the unit square, edge `2k` (`V_k → M_k`) along `v = 0`, spoke `k`
//! (`M_k → C`) up `u = 1`, spoke `k − 1` along `v = 1` and edge `2k − 1` run back up `u = 0`. Each
//! hole edge is its patch's boundary row exactly, and each spoke one quintic shared by the two
//! patches beside it.
//!
//! A biquintic net is four 3 × 3 corner blocks, and the block at a corner is the patch's
//! derivatives `∂ᵘⁱ∂ᵛʲS`, `i, j ≤ 2` there (Hermite data). Along each side of a patch, its
//! position, cross derivative and second cross derivative are quintics fixed by the blocks at the
//! side's two ends; so the patches are G2 to everything round them when, along every side, those
//! three curves meet the side's conditions, and the blocks at every corner are the ones both sides
//! there ask for (the corner compatibility of Hermite data: the twist and its derivatives).
//!
//! - **Beside a planar face** a patch is G2 exactly when its position, cross derivative and second
//!   cross derivative lie in the plane: its rows 0–2 are in it.
//! - **Beside a round** (the section moved along the rail `e`, to the second order) a patch is G2
//!   exactly when its cross derivatives are `a(t)·e` and `c(t)·e` (the patch is the round
//!   reparameterised to the second order across the section: `S(t, v) = R(t, a v + c v²/2)`, the
//!   round `R(t, w) = E(t) + w e`; its normal curvature along `e` is the round's, zero, and its
//!   twist term the round's, zero). `a` and `c` are scalar quintics; at the hole's corner `a e` is
//!   the face curve's tangent there and `c e` its second derivative (both along the contact line).
//! - **Across a spoke**, between patches `A = S_k` (`u` across) and `B = S_{k+1}` (`v` across),
//!   the regular reparameterisation: `B_v = −A_u + β A_t`, `B_vv = A_uu − 2β A_ut + β² A_tt + δ A_t`
//!   with `β = 2cos(2π/6)·t = t` and `δ = t` (`t` along the spoke from `M_k`): G1 and G2 by a
//!   reparameterisation (`B = A ∘ φ` to the second order), each side a quintic. At `M_k` it ties
//!   the two halves of the side's cross derivatives together: their slopes towards `M_k` add up to
//!   the spoke's tangent there (a round's `a` has a kink at `M_k`).
//! - **The centre** is a G2 enclosure: the patches' blocks there follow one another round `C`
//!   by the spokes' relations; the first patch's block is free but for the closure, and its
//!   orders 0–2 (the tangent plane, one quadratic form and the regular star of spoke tangents,
//!   opposite spokes in one line) close round the centre by themselves. Orders 3–4 close but for
//!   one vector condition on the patch's own third derivatives and one on the sides' middles (an
//!   alternating sum of the cross derivatives' slopes there, zero at a symmetric corner); the
//!   first is met by the centre's block, the second by the planar sides' second cross
//!   derivatives' slopes at `M`.
//!
//! What the conditions leave free is chosen by explicit rules, never fitted: the centre on the
//! curves between opposite middles, the spoke tangents the regular star of the centre's chords to
//! the middles, the quadratic form from the spokes' heights over the tangent plane (each the cubic
//! through its ends), each spoke's speed at its middle its chord, its second derivative there the
//! quartic's through its other end data, each round side's `a` and `c` the cubics through their
//! end data, and the centre's third derivatives those that make the spokes' cross derivatives
//! quartics on average round the centre. Nothing is solved but four small exact systems (3 × 3 and
//! one 9 × 9), and every condition holds to the input's own rounding.
//!
//! Input is hostile: sizes are fixed and every value checked; nothing here panics.

use crate::Vec3;
use crate::bspline::{HPoint, NurbsSurface, SplineError, SplineResult, invalid};

/// Sides of the hole: three rounds' sections and three faces' curves, alternately.
pub const QUAD_SIDES: usize = 6;
const N: usize = QUAD_SIDES;
/// The patches' degree each way.
pub const QUAD_DEGREE: usize = 5;
/// The regular reparameterisation across a spoke: `β(t) = BETA·t` (`2cos(2π/6)`) and
/// `δ(t) = DELTA·t`, which closes the centre's orders 0–2 by themselves.
const BETA: f64 = 1.0;
const DELTA: f64 = 1.0;

/// What a hole edge has beside it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Beside {
    /// The planar face `normal·p = offset` (`normal` a unit vector).
    Plane { normal: Vec3, offset: f64 },
    /// A round that, beside the edge (one of its sections), is the section moved along `rail` to
    /// the second order: its cross derivative along `rail` everywhere on the edge, with no normal
    /// curvature along it and no twist. `rail` is a unit vector into the hole.
    Rail { rail: Vec3 },
}

/// One edge of the hole, from the loop's vertex before it to the one after, on `[0, 1]`: one
/// polynomial span of degree at most [`QUAD_DEGREE`] (`weights` empty or all one).
#[derive(Clone, Debug, PartialEq)]
pub struct QuadEdge {
    pub degree: usize,
    pub knots: Vec<f64>,
    pub points: Vec<Vec3>,
    pub weights: Vec<f64>,
    pub beside: Beside,
}

/// The tolerances the hole's own data must meet the construction's conditions within.
#[derive(Clone, Copy, Debug)]
pub struct QuadSpec {
    /// How far apart two edges' ends may be where they meet, and an edge from its plane.
    pub tol: f64,
    /// How far (radians) a direction may be from the one the conditions ask for: a face curve's
    /// tangent from its round's rail at the hole's corner, a section's from its plane there, the
    /// two halves of a side at its middle.
    pub angle_tol: f64,
    /// How far (1/length) a curvature may be from the one the conditions ask for: a section's
    /// across its plane at the hole's corner (the eased setback's zero), a face curve's at the
    /// contact line.
    pub curvature_tol: f64,
}

/// A spoke: from its side's middle `M_k` (t = 0) to the centre (t = 1).
#[derive(Clone, Debug, PartialEq)]
pub struct Spoke {
    pub degree: usize,
    pub knots: Vec<f64>,
    pub points: Vec<Vec3>,
}

/// The fill: patch `k` on the unit square as the module documentation lays it out.
#[derive(Clone, Debug, PartialEq)]
pub struct QuadFill {
    pub patches: Vec<NurbsSurface>,
    /// Spoke `k` from `M_k` to the centre.
    pub spokes: Vec<Spoke>,
    pub centre: Vec3,
}

/// A curve's value and first two derivatives at a point.
type Jet = [Vec3; 3];
/// A patch's derivatives at a corner: `b[i][j] = ∂ᵘⁱ∂ᵛʲS`.
type Block = [[Vec3; 3]; 3];
/// A scalar's value and first two derivatives at a point.
type SJet = [f64; 3];
const ZERO_BLOCK: Block = [[Vec3::ZERO; 3]; 3];

/// Patch `k`'s corner blocks: at `V_k` (0, 0), `M_k` (1, 0), `M_{k−1}` (0, 1) and the centre.
#[derive(Clone, Copy)]
struct Blocks {
    v: Block,
    m: Block,
    p: Block,
    c: Block,
}

/// Fill the hole `edges` (twelve of them, each starting where the one before ends; edge `2k`
/// starts at a corner of the hole and meets edge `2k + 1` smoothly at its side's middle; the
/// sides alternately beside a round and a planar face) with six biquintic patches G2 to the
/// faces beside the edges and to each other (module docs). Holes of other shapes, rational or
/// multi-span edges, and data that misses the construction's conditions beyond `spec` are
/// refused.
pub fn quad_fill(edges: &[QuadEdge], spec: &QuadSpec) -> SplineResult<QuadFill> {
    let positive = |x: f64| x.is_finite() && x > 0.0;
    if !(positive(spec.tol) && positive(spec.angle_tol) && positive(spec.curvature_tol)) {
        return invalid("the fill's tolerances must be positive");
    }
    if edges.len() != 2 * N {
        return Err(SplineError::Unsupported(format!(
            "a hole of {} edges (the explicit corner takes six sides, three rounds' and three faces', two edges a side)",
            edges.len()
        )));
    }
    let hole = Hole::new(edges, spec)?;
    hole.fill()
}

/// A half side, from its hole corner `V` (t = 0) to its side's middle `M` (t = 1): a quintic
/// Bézier and what is beside it.
#[derive(Clone, Copy)]
struct Half {
    bez: [Vec3; 6],
    beside: Beside,
}

impl Half {
    /// The value and first two derivatives at the corner (`end` false) or the middle.
    fn jet(&self, end: bool) -> Jet {
        let b = &self.bez;
        if end {
            [b[5], (b[5] - b[4]) * 5.0, (b[5] - b[4] * 2.0 + b[3]) * 20.0]
        } else {
            [b[0], (b[1] - b[0]) * 5.0, (b[2] - b[1] * 2.0 + b[0]) * 20.0]
        }
    }
}

/// The hole as the construction reads it: per side `k`, its first half `a[k]` (edge `2k`,
/// `V_k → M_k`, patch `k`'s bottom) and its second `b[k]` (edge `2k + 1` run back,
/// `V_{k+1} → M_k`, patch `k + 1`'s left side).
struct Hole {
    a: [Half; N],
    b: [Half; N],
    spec: QuadSpec,
}

/// An edge as one quintic Bézier, raised from a lower degree when it has one.
fn quintic(e: &QuadEdge) -> SplineResult<[Vec3; 6]> {
    let p = e.degree;
    if !(1..=QUAD_DEGREE).contains(&p) {
        return Err(SplineError::Unsupported(format!("a hole edge of degree {p} (the explicit corner takes one span of degree 5 at most)")));
    }
    let one_span = e.knots.len() == 2 * (p + 1) && e.knots.iter().take(p + 1).all(|k| *k == 0.0) && e.knots.iter().skip(p + 1).all(|k| *k == 1.0);
    if !one_span || e.points.len() != p + 1 {
        return Err(SplineError::Unsupported("a hole edge of more than one span, or not on [0, 1] (the explicit corner takes one quintic span)".into()));
    }
    if !(e.weights.is_empty() || e.weights.iter().all(|w| *w == 1.0)) {
        return Err(SplineError::Unsupported(
            "a rational hole edge (the explicit corner takes polynomial edges: an eased setback's sections and face curves)".into(),
        ));
    }
    if !e.points.iter().all(|q| q.is_finite()) {
        return invalid("a hole edge's control points are not finite");
    }
    // Degree elevation of a Bézier: Q_i = (i/(p+1)) P_{i−1} + (1 − i/(p+1)) P_i.
    let mut pts = e.points.clone();
    while pts.len() < QUAD_DEGREE + 1 {
        let q = pts.len();
        let raised: Vec<Vec3> = (0..=q)
            .map(|i| {
                let f = i as f64 / q as f64;
                let before = if i == 0 { Vec3::ZERO } else { pts.get(i - 1).copied().unwrap_or(Vec3::ZERO) };
                let here = pts.get(i).copied().unwrap_or(Vec3::ZERO);
                before * f + here * (1.0 - f)
            })
            .collect();
        pts = raised;
    }
    let mut out = [Vec3::ZERO; 6];
    for (o, q) in out.iter_mut().zip(&pts) {
        *o = *q;
    }
    Ok(out)
}

/// The angle (radians) between two vectors, `π` where either vanishes.
fn angle(a: Vec3, b: Vec3) -> f64 {
    let d = a.len() * b.len();
    if !(d > 0.0) {
        return std::f64::consts::PI;
    }
    (a.cross(b).len() / d).atan2(a.dot(b) / d)
}

/// The patch net's control points from a corner block: `P_a = Σ_{i ≤ a} c[a][i] ∂ⁱS` along one
/// direction at the near end (`c0`) and at the far end (`c1`, counted from it).
const C0: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [1.0, 0.2, 0.0], [1.0, 0.4, 0.05]];
const C1: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [1.0, -0.2, 0.0], [1.0, -0.4, 0.05]];

/// The quintic's third and fourth derivatives at `t = 1` from its Hermite data at both ends.
fn third_at_one(f0: Jet, f1: Jet) -> Vec3 {
    f0[0] * -60.0 + f0[1] * -24.0 + f0[2] * -3.0 + f1[0] * 60.0 + f1[1] * -36.0 + f1[2] * 9.0
}
fn fourth_at_one(f0: Jet, f1: Jet) -> Vec3 {
    f0[0] * -360.0 + f0[1] * -168.0 + f0[2] * -24.0 + f1[0] * 360.0 + f1[1] * -192.0 + f1[2] * 36.0
}

/// `f''(1)` of the quartic with `f, f', f''` at 0 (`f0`) and `f, f'` at 1: the quintic's whose
/// fifth derivative vanishes.
fn quartic_second_at_one(f0: Jet, f1: [Vec3; 2]) -> Vec3 {
    f0[0] * 12.0 + f0[1] * 6.0 + f0[2] - f1[0] * 12.0 + f1[1] * 6.0
}

/// The cubic with `p(0) = v0` and `p, p', p''` at 1: its first two derivatives at 0.
fn cubic_start(v0: f64, at1: SJet) -> [f64; 2] {
    let [v1, d1, dd1] = at1;
    let c = v1 - d1 + dd1 / 2.0 - v0;
    [d1 - dd1 + 3.0 * c, dd1 - 6.0 * c]
}

/// Solve the small dense system `a x = b` (row-major `n × n`) by Gaussian elimination with
/// partial pivoting; `None` when a pivot is below `1e-12` of the matrix's largest entry.
fn dense_solve(mut a: Vec<f64>, mut b: Vec<f64>, n: usize) -> Option<Vec<f64>> {
    if a.len() != n * n || b.len() != n {
        return None;
    }
    let scale = a.iter().fold(0.0f64, |m, x| m.max(x.abs()));
    if !(scale > 0.0 && scale.is_finite()) {
        return None;
    }
    for c in 0..n {
        let pivot = (c..n).max_by(|i, j| a.get(i * n + c).map(|x| x.abs()).unwrap_or(0.0).total_cmp(&a.get(j * n + c).map(|x| x.abs()).unwrap_or(0.0)))?;
        if !(a.get(pivot * n + c).copied().unwrap_or(0.0).abs() > 1e-12 * scale) {
            return None;
        }
        if pivot != c {
            for k in 0..n {
                a.swap(pivot * n + k, c * n + k);
            }
            b.swap(pivot, c);
        }
        let d = a.get(c * n + c).copied()?;
        for r in 0..n {
            if r == c {
                continue;
            }
            let f = a.get(r * n + c).copied()? / d;
            if f == 0.0 {
                continue;
            }
            for k in c..n {
                let v = a.get(c * n + k).copied()?;
                if let Some(x) = a.get_mut(r * n + k) {
                    *x -= f * v;
                }
            }
            let bc = b.get(c).copied()?;
            if let Some(x) = b.get_mut(r) {
                *x -= f * bc;
            }
        }
    }
    let x: Vec<f64> = (0..n).map(|i| b.get(i).copied().unwrap_or(f64::NAN) / a.get(i * n + i).copied().unwrap_or(f64::NAN)).collect();
    x.iter().all(|v| v.is_finite()).then_some(x)
}

impl Hole {
    fn new(edges: &[QuadEdge], spec: &QuadSpec) -> SplineResult<Hole> {
        let mut nets = [[Vec3::ZERO; 6]; 2 * N];
        for (slot, e) in nets.iter_mut().zip(edges) {
            *slot = quintic(e)?;
        }
        let beside = |j: usize| edges.get(j).map(|e| e.beside).ok_or_else(|| SplineError::Invalid("edge".into()));
        let mut a = [Half { bez: [Vec3::ZERO; 6], beside: Beside::Rail { rail: Vec3::X } }; N];
        let mut b = a;
        for k in 0..N {
            let (first, second) = (nets[2 * k], nets[2 * k + 1]);
            let mut back = second;
            back.reverse();
            a[k] = Half { bez: first, beside: beside(2 * k)? };
            b[k] = Half { bez: back, beside: beside(2 * k + 1)? };
        }
        // Each edge starts where the one before ends.
        for j in 0..2 * N {
            let (prev, here) = (nets[(j + 2 * N - 1) % (2 * N)], nets[j]);
            let gap = prev[5].dist(here[0]);
            if !(gap <= spec.tol) {
                return invalid(format!("hole edge {j} does not start where edge {} ends ({gap:e} apart)", (j + 2 * N - 1) % (2 * N)));
            }
        }
        let hole = Hole { a, b, spec: *spec };
        hole.check()?;
        Ok(hole)
    }

    /// The hole's data against the construction's conditions (module docs), within the spec.
    fn check(&self) -> SplineResult<()> {
        let s = &self.spec;
        let unsupported = |m: String| Err(SplineError::Unsupported(m));
        for k in 0..N {
            let (a, b) = (self.a[k], self.b[k]);
            // One side, one thing beside it.
            match (a.beside, b.beside) {
                (Beside::Plane { normal: n1, offset: c1 }, Beside::Plane { normal: n2, offset: c2 }) => {
                    if !(angle(n1, n2) <= s.angle_tol && (c1 - c2).abs() <= s.tol && (n1.len() - 1.0).abs() <= 1e-9) {
                        return unsupported(format!("side {k}'s two halves lie on different planes"));
                    }
                    for (h, which) in [(a, "first"), (b, "second")] {
                        let off = h.bez.iter().map(|q| (n1.dot(*q) - c1).abs()).fold(0.0f64, f64::max);
                        if !(off <= s.tol) {
                            return unsupported(format!("side {k}'s {which} half leaves its plane by {off:e}"));
                        }
                    }
                }
                (Beside::Rail { rail: r1 }, Beside::Rail { rail: r2 }) => {
                    if !(angle(r1, r2) <= s.angle_tol && (r1.len() - 1.0).abs() <= 1e-9) {
                        return unsupported(format!("side {k}'s two halves have rounds of different rails"));
                    }
                }
                _ => return unsupported(format!("side {k}'s two halves have different faces beside them")),
            }
            let next = self.a[(k + 1) % N].beside;
            if matches!((a.beside, next), (Beside::Plane { .. }, Beside::Plane { .. }) | (Beside::Rail { .. }, Beside::Rail { .. })) {
                return unsupported(
                    "a hole whose sides do not alternate between a round's section and a planar face's curve (the explicit corner of three rounds)"
                        .into(),
                );
            }
            // The two halves meet smoothly at the middle, at the same speed (the side split at
            // its parameter's middle).
            let (ja, jb) = (a.jet(true), b.jet(true));
            if !(ja[0].dist(jb[0]) <= s.tol) {
                return invalid(format!("side {k}'s halves do not meet at its middle"));
            }
            let speed = ja[1].len();
            if !((ja[1] + jb[1]).len() <= s.angle_tol * speed) {
                return unsupported(format!(
                    "side {k}'s halves meet at its middle at different speeds or directions ({:e} apart): the explicit corner takes a side split at its middle",
                    (ja[1] + jb[1]).len()
                ));
            }
            if !((ja[2] - jb[2]).len() <= s.curvature_tol * speed * speed) {
                return unsupported(format!("side {k}'s halves meet at its middle with different curvatures"));
            }
        }
        // At each hole corner a round's section meets a face's curve: the curve leaves along the
        // rail with no curvature, and the section leaves tangent to the plane with no curvature
        // across it.
        for k in 0..N {
            let (bottom, left) = (self.a[k], self.b[(k + N - 1) % N]);
            let (round, face) = match (bottom.beside, left.beside) {
                (Beside::Rail { rail }, Beside::Plane { normal, .. }) => ((bottom, rail), (left, normal)),
                (Beside::Plane { normal, .. }, Beside::Rail { rail }) => ((left, rail), (bottom, normal)),
                _ => return unsupported("a hole corner between two faces of one kind".into()),
            };
            let (r, e) = (round.0.jet(false), round.1);
            let (f, n) = (face.0.jet(false), face.1);
            let (rs, fs) = (r[1].len(), f[1].len());
            if !(rs > 0.0 && fs > 0.0) {
                return invalid(format!("a hole edge stops at the hole's corner {k}"));
            }
            if !(angle(f[1], e) <= s.angle_tol) {
                return unsupported(format!(
                    "at the hole's corner {k} the face's curve leaves {:.4}° off its round's rail (the explicit corner takes a face curve tangent to the contact line)",
                    angle(f[1], e).to_degrees()
                ));
            }
            let bend = (f[2] - e * f[2].dot(e)).len() / (fs * fs);
            if !(bend <= s.curvature_tol) {
                return unsupported(format!(
                    "at the hole's corner {k} the face's curve turns ({bend:e} per length): the explicit corner takes one that leaves the contact line curvature-continuously (an eased setback)"
                ));
            }
            if !(r[1].dot(n).abs() <= s.angle_tol * rs) {
                return unsupported(format!("at the hole's corner {k} the round's section leaves its plane"));
            }
            let across = r[2].dot(n).abs() / (rs * rs);
            if !(across <= s.curvature_tol) {
                return unsupported(format!(
                    "at the hole's corner {k} the round's section curves across its plane ({across:e} per length): the explicit corner takes an eased section, flat across"
                ));
            }
        }
        Ok(())
    }

    /// The direction a side's cross derivative takes at its middle, into the hole: a round's rail;
    /// on a plane, square to the curve in the plane, to the side the hole's corner shows.
    fn leave(&self, k: usize) -> SplineResult<Vec3> {
        let a = self.a[k];
        match a.beside {
            Beside::Rail { rail } => Ok(rail),
            Beside::Plane { normal, .. } => {
                let t = a.jet(true)[1];
                let l = normal.cross(t).normalized().ok_or_else(|| SplineError::Invalid(format!("side {k} stops at its middle")))?;
                // At V_k the hole lies where the left edge leaves, across the bottom one.
                let (t0, into) = (a.jet(false)[1], self.b[(k + N - 1) % N].jet(false)[1]);
                Ok(if t0.cross(into).dot(t.cross(l)) >= 0.0 { l } else { -l })
            }
        }
    }

    fn fill(&self) -> SplineResult<QuadFill> {
        let ms: [Vec3; N] = std::array::from_fn(|k| self.a[k].jet(true)[0]);
        let mut ell = [Vec3::ZERO; N];
        for (k, l) in ell.iter_mut().enumerate() {
            *l = self.leave(k)?;
        }
        // The centre: the middle of each curve between opposite middles (the cubic leaving each
        // along its side's direction, its handles a third of the chord), averaged.
        let mut centre = Vec3::ZERO;
        for k in 0..N / 2 {
            let (p, q) = (ms[k], ms[k + N / 2]);
            centre += (p + q) * 0.5 + (ell[k] + ell[k + N / 2]) * (p.dist(q) / 8.0);
        }
        let centre = centre * (2.0 / N as f64);
        // The spoke tangents at the centre (towards it): the regular star, an affine image of the
        // regular hexagon, closest to the chords from the middles (its first Fourier mode).
        let ang = |k: usize| std::f64::consts::TAU * k as f64 / N as f64;
        let (mut x, mut y) = (Vec3::ZERO, Vec3::ZERO);
        for (k, m) in ms.iter().enumerate() {
            let u = centre - *m;
            x += u * (ang(k).cos() * 2.0 / N as f64);
            y += u * (ang(k).sin() * 2.0 / N as f64);
        }
        let star: [Vec3; N] = std::array::from_fn(|k| x * ang(k).cos() + y * ang(k).sin());
        let normal = star[N - 1].cross(star[0]).normalized().ok_or_else(|| SplineError::Invalid("the spokes' star at the centre is flat".into()))?;
        // The quadratic form at the centre: along each spoke, its height over the tangent plane is
        // the cubic through its ends (h, h' at M; 0, 0 at C), whose h''(1) = 6 h(0) + 2 h'(0) is
        // the form on that spoke's tangent; opposite spokes averaged, three directions solved.
        let e1 = star[0].normalized().ok_or_else(|| SplineError::Invalid("a spoke has no tangent at the centre".into()))?;
        let e2 = normal.cross(e1);
        let mut rows = Vec::with_capacity(9);
        let mut rhs = Vec::with_capacity(3);
        for (k, t) in star.iter().take(N / 2).enumerate() {
            let mut h = 0.0;
            for kk in [k, k + N / 2] {
                let chord = centre.dist(ms[kk]);
                h += 6.0 * normal.dot(ms[kk] - centre) + 2.0 * normal.dot(ell[kk]) * chord;
            }
            let (u, v) = (t.dot(e1), t.dot(e2));
            rows.extend([u * u, 2.0 * u * v, v * v]);
            rhs.push(h / 2.0);
        }
        let form = dense_solve(rows, rhs, 3).ok_or_else(|| SplineError::Invalid("the spokes' directions at the centre do not fix a quadratic form".into()))?;
        let ii = |p: Vec3, q: Vec3| {
            let (pu, pv, qu, qv) = (p.dot(e1), p.dot(e2), q.dot(e1), q.dot(e2));
            form.first().copied().unwrap_or(0.0) * pu * qu
                + form.get(1).copied().unwrap_or(0.0) * (pu * qv + pv * qu)
                + form.get(2).copied().unwrap_or(0.0) * pv * qv
        };
        // Each spoke's second derivative at the centre: the form's normal part and the regular
        // star's tangential half of the tangent (what closes orders 0–2 round the centre).
        let spoke_second_c: [Vec3; N] = std::array::from_fn(|k| star[k] * 0.5 + normal * ii(star[k], star[k]));
        // Each spoke at its middle: its tangent along the side's direction, as long as its chord;
        // its second derivative the quartic's through the rest of its end data, in what the side
        // allows (along the rail; in the plane).
        let mut spoke_m = [[Vec3::ZERO; 3]; N];
        for (k, s) in spoke_m.iter_mut().enumerate() {
            let chord = centre.dist(ms[k]);
            let d1 = ell[k] * chord;
            let raw = (centre - ms[k]) * 12.0 - d1 * 6.0 - star[k] * 6.0 + spoke_second_c[k];
            let d2 = match self.a[k].beside {
                Beside::Rail { rail } => rail * raw.dot(rail),
                Beside::Plane { normal: n, .. } => raw - n * raw.dot(n),
            };
            *s = [ms[k], d1, d2];
        }
        // The rounds' cross derivative speeds a(t) and second derivatives c(t) along each half,
        // and the planar halves' cross derivatives at the middle.
        let mut m_blocks = [ZERO_BLOCK; N];
        let mut round_v: [Option<(SJet, SJet)>; 2 * N] = [None; 2 * N];
        for k in 0..N {
            let a = self.a[k];
            let ja = a.jet(true);
            let [_, d1, d2] = spoke_m[k];
            let mut m = ZERO_BLOCK;
            for (row, jet) in m.iter_mut().zip(ja) {
                row[0] = jet;
            }
            match a.beside {
                Beside::Rail { rail } => {
                    // The face curves' tangent and second derivative along the rail at the two corners.
                    let (fa, fb) = (self.b[(k + N - 1) % N].jet(false), self.a[(k + 1) % N].jet(false));
                    let (av, bv) = (fa[1].dot(rail), fb[1].dot(rail));
                    let (cav, cbv) = (fa[2].dot(rail), fb[2].dot(rail));
                    if !(av > 0.0 && bv > 0.0) {
                        return invalid(format!("side {k}'s face curves do not leave along its rail"));
                    }
                    let am = d1.dot(rail);
                    // The halves' slopes towards M add up to the spoke's speed there (β); the
                    // curvature terms tie as the relation's second order asks.
                    let (sa, sb) = ((am + bv - av) / 2.0, (am + av - bv) / 2.0);
                    let (qa, qb) = ((sa - sb) / 2.0, (sb - sa) / 2.0);
                    let cm = d2.dot(rail);
                    let (ca, cb) = (cm + (cbv - cav) / 2.0, cm + (cav - cbv) / 2.0);
                    let (cqa, cqb) = (ca - cb, cb - ca);
                    let a_m: SJet = [am, sa, qa];
                    let c_m: SJet = [cm, ca, cqa];
                    let start = |v0: f64, at1: SJet| {
                        let [p1, p2] = cubic_start(v0, at1);
                        [v0, p1, p2]
                    };
                    round_v[2 * k] = Some((start(av, a_m), start(cav, c_m)));
                    round_v[2 * k + 1] = Some((start(bv, [am, sb, qb]), start(cbv, [cm, cb, cqb])));
                    for i in 0..3 {
                        m[i][1] = rail * a_m[i];
                        m[i][2] = rail * c_m[i];
                    }
                }
                Beside::Plane { .. } => {
                    m[0][1] = d1;
                    m[1][1] = d1 * (BETA / 2.0);
                    m[0][2] = d2;
                    m[1][2] = d2 * BETA;
                }
            }
            m_blocks[k] = m;
        }
        // The sides' middles' share of the centre's closure (order 4): the sum over the patches of
        // 120 E'(1) + 48 D'(1) + 6 Q'(1) − 24 D(1) − 6 Q(1) must vanish; the planar halves' Q'(1)
        // take what it misses, each in its plane.
        let miss = m_blocks.iter().fold(Vec3::ZERO, |s, m| s + m[1][0] * 120.0 + m[1][1] * 48.0 + m[1][2] * 6.0 - m[0][1] * 24.0 - m[0][2] * 6.0);
        if miss.len() > 0.0 {
            let mut g = [[0.0f64; 3]; 3];
            let mut planes = Vec::new();
            for (k, h) in self.a.iter().enumerate() {
                if let Beside::Plane { normal: n, .. } = h.beside {
                    let nv = [n.x, n.y, n.z];
                    for (r, row) in g.iter_mut().enumerate() {
                        for (c, x) in row.iter_mut().enumerate() {
                            *x += if r == c { 1.0 } else { 0.0 } - nv[r] * nv[c];
                        }
                    }
                    planes.push((k, n));
                }
            }
            let flat: Vec<f64> = g.iter().flatten().copied().collect();
            let target = miss * (-1.0 / 6.0);
            let y = dense_solve(flat, vec![target.x, target.y, target.z], 3)
                .ok_or_else(|| SplineError::Invalid("the hole's planes do not span space".into()))?;
            let yv = Vec3::new(y.first().copied().unwrap_or(0.0), y.get(1).copied().unwrap_or(0.0), y.get(2).copied().unwrap_or(0.0));
            for (k, n) in planes {
                if let Some(m) = m_blocks.get_mut(k) {
                    m[1][2] += yv - n * yv.dot(n);
                }
            }
        }
        // The blocks at the hole's corners, and at the middles on the far side of each spoke.
        let mut blocks = [Blocks { v: ZERO_BLOCK, m: ZERO_BLOCK, p: ZERO_BLOCK, c: ZERO_BLOCK }; N];
        for k in 0..N {
            let prev = (k + N - 1) % N;
            let (bottom, left) = (self.a[k], self.b[prev]);
            let (jb, jl) = (bottom.jet(false), left.jet(false));
            let mut v = ZERO_BLOCK;
            for i in 0..3 {
                v[i][0] = jb[i];
                v[0][i] = jl[i];
            }
            // The round's cross fields' derivatives at V (along its edge), the plane's following.
            if let (Beside::Rail { rail }, Some((a, c))) = (bottom.beside, round_v.get(2 * k).copied().flatten()) {
                v[1][1] = rail * a[1];
                v[2][1] = rail * a[2];
                v[1][2] = rail * c[1];
                v[2][2] = rail * c[2];
            } else if let (Beside::Rail { rail }, Some((a, c))) = (left.beside, round_v.get(2 * prev + 1).copied().flatten()) {
                v[1][1] = rail * a[1];
                v[1][2] = rail * a[2];
                v[2][1] = rail * c[1];
                v[2][2] = rail * c[2];
            } else {
                return invalid(format!("the hole's corner {k} has no round beside it"));
            }
            blocks[k].v = v;
            blocks[k].m = m_blocks[k];
            // Patch k + 1's block at M_k: the spoke's relation from patch k's, its edge row the
            // hole's own.
            let mut p = relation_at_middle(&m_blocks[k]);
            let jn = self.b[k].jet(true);
            for (j, x) in jn.iter().enumerate() {
                p[0][j] = *x;
            }
            blocks[(k + 1) % N].p = p;
        }
        // The centre's blocks: the first patch's orders 0–2 from the centre, the star and the form;
        // its orders 3–4 from the closure and the spokes' quartic rule; the others by the spokes.
        let mut c0 = ZERO_BLOCK;
        let (t5, t0) = (star[N - 1], star[0]);
        c0[0][0] = centre;
        c0[1][0] = t5;
        c0[0][1] = t0;
        c0[2][0] = spoke_second_c[N - 1];
        c0[0][2] = spoke_second_c[0];
        c0[1][1] = (t5 + t0) * 0.5 + normal * ii(t5, t0);
        let around = |c21: Vec3, c12: Vec3, c22: Vec3| -> [Block; N + 1] {
            let mut c = c0;
            c[2][1] = c21;
            c[1][2] = c12;
            c[2][2] = c22;
            let mut out = [ZERO_BLOCK; N + 1];
            out[0] = c;
            for k in 0..N {
                c = relation_to_centre(&m_blocks[k], &c);
                out[k + 1] = c;
            }
            out
        };
        // The nine conditions on (c21, c12, c22), affine: the closure of the order-3 block's
        // c21 row, the alternating sum of the spokes' col_1''(1) misses from their quartics, and
        // the sum of their col_2''(1) misses.
        let conditions = |v: &[f64; 9]| -> [f64; 9] {
            let vec = |i: usize| Vec3::new(v[i], v[i + 1], v[i + 2]);
            let all = around(vec(0), vec(3), vec(6));
            let close = all[N][2][1] - all[0][2][1];
            let (mut alt, mut sum) = (Vec3::ZERO, Vec3::ZERO);
            for (k, (m, c)) in m_blocks.iter().zip(all.iter()).enumerate() {
                let p12 = quartic_second_at_one(m[1], [c[1][0], c[1][1]]);
                let p22 = quartic_second_at_one(m[2], [c[2][0], c[2][1]]);
                let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
                alt += (c[1][2] - p12) * sign;
                sum += c[2][2] - p22;
            }
            [close.x, close.y, close.z, alt.x, alt.y, alt.z, sum.x, sum.y, sum.z]
        };
        let base = conditions(&[0.0; 9]);
        let mut a = vec![0.0; 81];
        for q in 0..9 {
            let mut v = [0.0; 9];
            v[q] = 1.0;
            let r = conditions(&v);
            for (row, (x, y)) in r.iter().zip(base.iter()).enumerate() {
                if let Some(slot) = a.get_mut(row * 9 + q) {
                    *slot = x - y;
                }
            }
        }
        let sol = dense_solve(a, base.iter().map(|x| -x).collect(), 9)
            .ok_or_else(|| SplineError::Invalid("the centre's third derivatives are not fixed by its closure".into()))?;
        let sv = |i: usize| Vec3::new(sol.get(i).copied().unwrap_or(0.0), sol.get(i + 1).copied().unwrap_or(0.0), sol.get(i + 2).copied().unwrap_or(0.0));
        let all = around(sv(0), sv(3), sv(6));
        // The closure, which the rules and the middles' share leave to the input's rounding.
        let scale = all[0].iter().flatten().fold(0.0f64, |m, x| m.max(x.len())).max(1.0);
        let off = (0..3).flat_map(|i| (0..3).map(move |j| (i, j))).map(|(i, j)| all[N][i][j].dist(all[0][i][j])).fold(0.0f64, f64::max);
        if !(off <= 1e-6 * scale) {
            return Err(SplineError::Unsupported(format!("the patches do not close round the centre ({off:e} off; the hole's sides miss each other's conditions)")));
        }
        for (b, c) in blocks.iter_mut().zip(all.iter()) {
            b.c = *c;
        }
        // The nets.
        let knots = vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
        let mut patches = Vec::with_capacity(N);
        for b in &blocks {
            let mut net = [[Vec3::ZERO; 6]; 6];
            for (blk, cu, cv, fu, fv) in [(&b.v, &C0, &C0, false, false), (&b.m, &C1, &C0, true, false), (&b.p, &C0, &C1, false, true), (&b.c, &C1, &C1, true, true)] {
                for ai in 0..3 {
                    for bj in 0..3 {
                        let mut q = Vec3::ZERO;
                        for i in 0..=ai {
                            for j in 0..=bj {
                                q += blk[i][j] * (cu[ai][i] * cv[bj][j]);
                            }
                        }
                        net[if fu { 5 - ai } else { ai }][if fv { 5 - bj } else { bj }] = q;
                    }
                }
            }
            let ctrl: Vec<HPoint> = (0..6).flat_map(|j| (0..6).map(move |i| (i, j))).map(|(i, j)| HPoint::weighted(net[i][j], 1.0)).collect();
            patches.push(NurbsSurface::from_homogeneous([QUAD_DEGREE, QUAD_DEGREE], [knots.clone(), knots.clone()], ctrl)?);
        }
        for (k, patch) in patches.iter().enumerate() {
            // Each patch's bottom row and left column are its hole edges' nets (to rounding).
            let bottom = self.a[k].bez;
            let left = self.b[(k + N - 1) % N].bez;
            let off = (0..6)
                .map(|i| {
                    let pb = patch.point(i, 0).project().unwrap_or(Vec3::new(f64::NAN, 0.0, 0.0));
                    let pl = patch.point(0, i).project().unwrap_or(Vec3::new(f64::NAN, 0.0, 0.0));
                    pb.dist(bottom[i]).max(pl.dist(left[i]))
                })
                .fold(0.0f64, |m, d| if d.is_nan() { f64::INFINITY } else { m.max(d) });
            if !(off <= self.spec.tol) {
                return invalid(format!("patch {k}'s boundary misses its hole edges by {off:e}"));
            }
        }
        let spokes = (0..N)
            .map(|k| {
                let (m, c) = (blocks[k].m, blocks[k].c);
                let (s0, s1) = ([m[0][0], m[0][1], m[0][2]], [c[0][0], c[0][1], c[0][2]]);
                let pts = vec![
                    s0[0],
                    s0[0] + s0[1] * 0.2,
                    s0[0] + s0[1] * 0.4 + s0[2] * 0.05,
                    s1[0] - s1[1] * 0.4 + s1[2] * 0.05,
                    s1[0] - s1[1] * 0.2,
                    s1[0],
                ];
                Spoke { degree: QUAD_DEGREE, knots: knots.clone(), points: pts }
            })
            .collect();
        Ok(QuadFill { patches, spokes, centre })
    }
}

/// Patch `k + 1`'s block at `M_k` from patch `k`'s by the spoke's relation (module docs) at its
/// start: `P[i][0] = M[0][i]`, `P[i][1] = −M[1][i] + β i M[0][i]`,
/// `P[i][2] = M[2][i] − 2β i M[1][i] + (β² i(i − 1) + δ i) M[0][i]`.
fn relation_at_middle(m: &Block) -> Block {
    let mut p = ZERO_BLOCK;
    for i in 0..3 {
        let fi = i as f64;
        p[i][0] = m[0][i];
        p[i][1] = -m[1][i] + m[0][i] * (BETA * fi);
        p[i][2] = m[2][i] - m[1][i] * (2.0 * BETA * fi) + m[0][i] * (BETA * BETA * fi * (fi - 1.0) + DELTA * fi);
    }
    p
}

/// Patch `k + 1`'s block at the centre from patch `k`'s block at `M_k` and at the centre: its
/// spoke-side curves (`∂ᵘᵐS(1, t)`, quintics from those two blocks) carried across the spoke by
/// the relation, read at `t = 1`.
fn relation_to_centre(m: &Block, c: &Block) -> Block {
    // col_m^(d)(1) for d = 0 … 4.
    let d = |mm: usize, k: usize| -> Vec3 {
        let (f0, f1) = (m[mm], c[mm]);
        match k {
            0..=2 => f1[k],
            3 => third_at_one(f0, f1),
            _ => fourth_at_one(f0, f1),
        }
    };
    let (b, dl) = (BETA, DELTA);
    let mut out = ZERO_BLOCK;
    for (i, row) in out.iter_mut().enumerate() {
        let fi = i as f64;
        row[0] = d(0, i);
        row[1] = -d(1, i) + d(0, i + 1) * b + d(0, i) * (b * fi);
        row[2] = d(2, i) - d(1, i + 1) * (2.0 * b) - d(1, i) * (2.0 * b * fi)
            + d(0, i + 2) * (b * b)
            + d(0, i + 1) * (2.0 * b * b * fi)
            + d(0, i) * (b * b * fi * (fi - 1.0))
            + d(0, i + 1) * dl
            + d(0, i) * (dl * fi);
    }
    out
}

#[cfg(test)]
#[path = "fill_quad_tests.rs"]
mod tests;
