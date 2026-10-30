//! Sheet operations, as `solidcraft_smlib`'s: sewing and a sheet's offset. Both are the C API's
//! glue round SMLib's routines (`smlib_sheets.cpp`), not ported yet, so both are refused.

use crate::{Brep, Result, refused};

impl Brep {
    /// The faces sewn along the vertices and edges that coincide within `tolerance` mm: not
    /// ported.
    pub fn stitch(&self, tolerance: f64) -> Result<(Brep, bool)> {
        let _ = tolerance;
        refused()
    }

    /// This sheet with every face moved onto its surface's offset: not ported.
    pub fn offset_sheet(&self, distance: f64) -> Result<Brep> {
        let _ = distance;
        refused()
    }
}
