//! Split Face and Split Body, as `solidcraft_smlib`'s. Their tools are [`Surface`]s, which this
//! kernel does not make yet, so neither call can be reached; both are refused.

use crate::ssi::Surface;
use crate::{Brep, Result, refused};

impl Brep {
    /// `self` with faces `faces` cut where each of `tools` meets them: not ported.
    pub fn split_faces(&self, faces: &[usize], tools: &[&Surface]) -> Result<Option<(Brep, Vec<usize>)>> {
        let _ = (faces, tools);
        refused()
    }

    /// The manifold solid `self` cut in two by `tool`: not ported.
    pub fn split_body(&self, tool: &Surface) -> Result<Option<(Brep, Brep)>> {
        let _ = tool;
        refused()
    }
}
