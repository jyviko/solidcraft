//! Untrim and extend of a sheet's faces, as `solidcraft_smlib`'s. Both are the fork's C API glue
//! (`smlib_untrim.cpp`, `smlib_extend.cpp`), not ported yet, as is the faces' internal points;
//! all three are refused.

use crate::{Brep, Result, refused};

/// How [`Brep::extend_sheet`] continues a surface past its side.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Continuation {
    /// The surface's own geometry (a B-spline's end spans continued).
    Natural,
    /// A run-off tangent to the surface across the side (G1).
    Tangent,
    /// A run-off that also keeps the surface's curvature across the side (G2).
    Curvature,
}

/// The inner loops [`Brep::untrim_face`] takes away.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Holes<'a> {
    All,
    /// The inner loop nearest each point (within the tolerance, else [`crate::Error::NoMatch`]).
    At(&'a [[f64; 3]]),
}

impl Brep {
    /// The face nearest `point` restored towards its surface: not ported.
    pub fn untrim_face(&self, point: [f64; 3], tolerance: f64, holes: Holes, outer: bool, in_place: bool) -> Result<Brep> {
        let _ = (point, tolerance, holes, outer, in_place);
        refused()
    }

    /// The one-face sheet `self` continued past a side: not ported.
    pub fn extend_sheet(&self, point: [f64; 3], tolerance: f64, distance: f64, how: Continuation) -> Result<Brep> {
        let _ = (point, tolerance, distance, how);
        refused()
    }

    /// A point inside each face: not ported.
    pub fn face_points(&self) -> Result<Vec<[f64; 3]>> {
        refused()
    }
}
