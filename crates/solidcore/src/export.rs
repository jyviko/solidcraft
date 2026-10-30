//! The kernel's B-reps out to SolidCraft's schema 4 B-rep (`solidcraft_geom::brep`), and the
//! primitive solids beyond the box, as `solidcraft_smlib::export`'s. There a body comes out as
//! SMLib's C API's flat arrays, read and mapped to schema 4; that bridge is not ported onto
//! solidcore ([`crate::flat`]), so no body comes out as a record yet. The primitives are SM_API's
//! own constructors, as the C API calls them (`smlib_primitives.cpp`).

use solidcraft_geom::brep as b;

use crate::history::Order;
use crate::import::Arrays;
use crate::{Brep, Error, Result, error, finite, kernel, point, refused};

use solidcore::primitive::api_primitives as primitives;

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}

impl Brep {
    /// A cylinder about +Z: base centre, radius and height (mm, > 0).
    pub fn make_cylinder(origin: [f64; 3], radius: f64, height: f64) -> Result<Brep> {
        if !finite(&origin) || !positive(radius) || !positive(height) {
            return Err(Error::Argument);
        }
        primitives::create_cylinder(&mut kernel(), point(origin), radius, height).map(Brep).map_err(error)
    }

    /// A sphere: centre and radius (mm, > 0).
    pub fn make_sphere(center: [f64; 3], radius: f64) -> Result<Brep> {
        if !finite(&center) || !positive(radius) {
            return Err(Error::Argument);
        }
        primitives::create_sphere(&mut kernel(), point(center), radius).map(Brep).map_err(error)
    }

    /// A cone or frustum about +Z: base centre, base radius (> 0), top radius (≥ 0, not the base
    /// radius) and height (> 0), in mm.
    pub fn make_cone(origin: [f64; 3], base_radius: f64, top_radius: f64, height: f64) -> Result<Brep> {
        if !finite(&origin)
            || !positive(base_radius)
            || !top_radius.is_finite()
            || !(top_radius >= 0.0)
            || top_radius == base_radius
            || !positive(height)
        {
            return Err(Error::Argument);
        }
        primitives::create_cone(&mut kernel(), point(origin), base_radius, top_radius, height).map(Brep).map_err(error)
    }

    /// A ring torus about +Z: centre, major radius > minor radius > 0 (mm).
    pub fn make_torus(center: [f64; 3], major_radius: f64, minor_radius: f64) -> Result<Brep> {
        if !finite(&center) || !positive(major_radius) || !positive(minor_radius) || !(minor_radius < major_radius) {
            return Err(Error::Argument);
        }
        primitives::create_torus(&mut kernel(), point(center), major_radius, minor_radius).map(Brep).map_err(error)
    }

    /// The flat arrays of this body: not ported (see the module).
    pub fn arrays(&self) -> Result<Arrays> {
        refused()
    }

    /// This body as a schema 4 B-rep, validated: not ported (see the module).
    pub fn to_native(&self) -> std::result::Result<b::Brep, ExportError> {
        self.to_native_report().map(|(brep, _)| brep)
    }

    /// [`Brep::to_native`], with what the mapping kept as the kernel wrote it: not ported.
    pub fn to_native_report(&self) -> std::result::Result<(b::Brep, Report), ExportError> {
        Err(ExportError::Smlib(Error::NotBuilt))
    }

    /// [`Brep::to_native`], with where each face and edge went in the record's canonical order:
    /// not ported.
    pub fn to_native_ordered(&self) -> std::result::Result<(b::Brep, Order), ExportError> {
        Err(ExportError::Smlib(Error::NotBuilt))
    }
}

/// What the mapping could not put in its own kind, or turned.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Report {
    /// NURBS edge curves kept as splines (neither a line nor a circular arc).
    pub nurbs_edges: usize,
    /// Faces whose loops were reversed (the kernel's natural normal against the record's).
    pub reversed_faces: usize,
}

/// Why a body has no schema 4 record. The variant names are `solidcraft_smlib`'s, which the
/// kernel builds against: [`ExportError::Smlib`] is the kernel call's failure.
#[derive(Clone, Debug, PartialEq)]
pub enum ExportError {
    /// The kernel call failed.
    Smlib(Error),
    /// The kernel holds something schema 4 (or this mapping) has no record for, named.
    Unsupported(String),
    /// The arrays are malformed or inconsistent.
    Invalid(String),
}

impl std::fmt::Display for ExportError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ExportError::Smlib(e) => write!(f, "solidcore export: {e}"),
            ExportError::Unsupported(what) => write!(f, "solidcore export: not supported yet: {what}"),
            ExportError::Invalid(what) => write!(f, "solidcore export: {what}"),
        }
    }
}

impl std::error::Error for ExportError {}

/// The schema 4 B-rep of flat arrays: not ported (see the module).
pub fn native(arrays: &Arrays) -> std::result::Result<(b::Brep, Report), ExportError> {
    let _ = arrays;
    Err(ExportError::Smlib(Error::NotBuilt))
}
