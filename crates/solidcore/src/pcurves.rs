//! The kernel's own pcurves for a face's uses that have no closed form, as
//! `solidcraft_smlib::pcurves` makes them: there the face goes to SMLib as a one-face sheet
//! through the flat arrays and comes back with SMLib's UV trims. That bridge is not ported onto
//! solidcore ([`crate::flat`]), so every face is refused.

use solidcraft_geom::brep::Curve2;
use solidcraft_geom::brep::pcurve::{TrimFace, Trims};

/// [`Trims`] by the kernel. The name is `solidcraft_smlib`'s, which the kernel builds against.
#[derive(Clone, Copy, Debug, Default)]
pub struct SmlibTrims;

impl Trims for SmlibTrims {
    fn trims(&self, face: &TrimFace<'_>) -> Result<Vec<Vec<Curve2>>, String> {
        let _ = face;
        Err(format!("solidcore pcurves: {}", crate::Error::NotBuilt))
    }
}
