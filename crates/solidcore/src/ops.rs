//! Modelling operations, as `solidcraft_smlib`'s: booleans, blends, shell, offset faces, thicken,
//! merge faces and draft. Each is the C API's glue round SMLib's routine (picking edges and faces
//! by points, copying the operands, checking the result, `smlib_ops.cpp`), not ported yet, so
//! each is refused. Only [`Brep::solid_count`] is SMLib's call alone.

use crate::{Brep, Result, error, refused};

/// How [`Brep::boolean`] combines two bodies.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BooleanOp {
    Union,
    /// The first body minus the second.
    Subtract,
    Intersect,
}

/// What [`Brep::blend`] makes of each selected edge.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Blend {
    /// A circular cross-section at a constant radius ([`Brep::fillet`]).
    Round,
    /// A curvature-continuous (G2) cross-section between the contact lines of the round of
    /// that radius.
    Curvature,
    /// A circular cross-section of constant width from contact line to contact line: the size
    /// is that width.
    Chord,
    /// An equal-distance chamfer ([`Brep::chamfer`]): the size is the setback along each face.
    Chamfer,
}

/// The second side of [`Brep::chamfer_sides`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ChamferSecond {
    /// The setback along the other face.
    Distance(f64),
    /// The angle (radians) the chamfer makes with the first face.
    Angle(f64),
}

/// Most bodies [`Brep::union_all`] takes (the C API's limit).
pub const MAX_UNION_BREPS: usize = 100_000;

/// Most points one selection may hold (the C API's limit).
pub const MAX_POINTS: usize = 4096;

/// One edge round an opened corner's hole ([`Brep::blend_setback_open`]): its index in the
/// sheet's edge order, the face beside it in its face order, and whether that face, run
/// counter-clockwise about its outward normal, runs the edge from its start vertex to its end.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct HoleEdge {
    pub edge: usize,
    pub face: usize,
    pub forward: bool,
}

/// The edges round an opened corner's hole, in order: each starts where the one before it ends.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Hole {
    pub edges: Vec<HoleEdge>,
}

impl Brep {
    /// `self` and `other` combined by `op`: not ported (see the module).
    pub fn boolean(&self, other: &Brep, op: BooleanOp) -> Result<Option<Brep>> {
        let _ = (other, op);
        refused()
    }

    /// The union of `breps` in one call: not ported.
    pub fn union_all(breps: &[&Brep]) -> Result<Brep> {
        let _ = breps;
        refused()
    }

    /// Two sheets of planar faces combined as regions of their plane: not ported.
    pub fn boolean2d(&self, other: &Brep, op: BooleanOp) -> Result<Option<Brep>> {
        let _ = (other, op);
        refused()
    }

    /// `self` with each face `tool` meets cut along the intersection: not ported.
    pub fn imprint(&self, tool: &Brep) -> Result<Option<Brep>> {
        let _ = tool;
        refused()
    }

    /// The number of separate solid regions, voids not counted (`SmApiBrepMaterialCensus`).
    pub fn solid_count(&self) -> Result<usize> {
        solidcore::query::api_queries::brep_material_census(&self.0).map(|(solids, _)| solids).map_err(error)
    }

    /// The edges nearest to `edges` rounded at `radius`: not ported.
    pub fn fillet(&self, edges: &[[f64; 3]], tolerance: f64, radius: f64) -> Result<Brep> {
        let _ = (edges, tolerance, radius);
        refused()
    }

    /// The edges nearest to `edges` cut by an equal-distance chamfer: not ported.
    pub fn chamfer(&self, edges: &[[f64; 3]], tolerance: f64, distance: f64) -> Result<Brep> {
        let _ = (edges, tolerance, distance);
        refused()
    }

    /// The edges nearest to `edges` blended as `kind`: not ported.
    pub fn blend(&self, kind: Blend, edges: &[[f64; 3]], sizes: &[f64], tolerance: f64) -> Result<Brep> {
        let _ = (kind, edges, sizes, tolerance);
        refused()
    }

    /// [`Brep::blend`] with corner setbacks: not ported.
    pub fn blend_setback(
        &self,
        kind: Blend,
        edges: &[[f64; 3]],
        sizes: &[f64],
        corners: &[[f64; 3]],
        setbacks: &[f64],
        tolerance: f64,
    ) -> Result<Brep> {
        let _ = (kind, edges, sizes, corners, setbacks, tolerance);
        refused()
    }

    /// [`Brep::blend_setback`] with each set-back corner left open: not ported (the open corner
    /// is the fork's own).
    pub fn blend_setback_open(
        &self,
        kind: Blend,
        edges: &[[f64; 3]],
        sizes: &[f64],
        corners: &[[f64; 3]],
        setbacks: &[f64],
        tolerance: f64,
    ) -> Result<(Brep, Vec<Hole>)> {
        let _ = (kind, edges, sizes, corners, setbacks, tolerance);
        refused()
    }

    /// The edge nearest to `edge` rounded with a varying radius: not ported.
    pub fn fillet_variable(&self, edge: [f64; 3], tolerance: f64, start: [f64; 3], degree: usize, knots: &[f64], radii: &[f64]) -> Result<Brep> {
        let _ = (edge, tolerance, start, degree, knots, radii);
        refused()
    }

    /// The edges nearest to `edges` cut by a two-sided chamfer: not ported.
    pub fn chamfer_sides(&self, edges: &[[f64; 3]], tolerance: f64, distance: f64, second: ChamferSecond, first: [f64; 3]) -> Result<Brep> {
        let _ = (edges, tolerance, distance, second, first);
        refused()
    }

    /// Hollowed to a wall of `thickness`: not ported.
    pub fn shell(&self, open_faces: &[[f64; 3]], tolerance: f64, thickness: f64) -> Result<Brep> {
        let _ = (open_faces, tolerance, thickness);
        refused()
    }

    /// The faces nearest to `faces` moved along their normals: not ported (the fork's face
    /// offset, `SmOffsetFacesExecutive`, is not in solidcore).
    pub fn offset_faces(&self, faces: &[[f64; 3]], tolerance: f64, distance: f64) -> Result<Brep> {
        let _ = (faces, tolerance, distance);
        refused()
    }

    /// A sheet thickened into a solid: not ported.
    pub fn thicken(&self, thickness: f64, symmetric: bool) -> Result<Brep> {
        let _ = (thickness, symmetric);
        refused()
    }

    /// This body with its topological edges and vertices removed: not ported.
    pub fn merge_faces(&self, keep: &[[f64; 3]], tolerance: f64) -> Result<Brep> {
        let _ = (keep, tolerance);
        refused()
    }

    /// The flat faces nearest to `faces` tilted about the neutral plane: not ported.
    pub fn draft(
        &self,
        faces: &[[f64; 3]],
        tolerance: f64,
        neutral_origin: [f64; 3],
        neutral_normal: [f64; 3],
        pull: [f64; 3],
        angle: f64,
    ) -> Result<Brep> {
        let _ = (faces, tolerance, neutral_origin, neutral_normal, pull, angle);
        refused()
    }
}
