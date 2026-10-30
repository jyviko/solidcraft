//! solidcore (tangency-labs/solidcore, Apache-2.0), the pure-Rust port of NVIDIA's SMLib, behind
//! the public API of `solidcraft-smlib`: the same functions, signatures, record types and error
//! kinds, so the kernel builds on either (`solidcraft-kernel`'s `smlib` or `solidcore` feature).
//! Only `crates/kernel` uses this crate.
//!
//! `solidcraft-smlib` is a Rust wrapper over the SolidCraft fork's C API (`capi/` of
//! tangency-labs/smlib), whose C++ glue picks edges and faces by points, lays bodies out as flat
//! arrays and checks every result. A call is implemented here where solidcore has the SMLib
//! routine that glue calls and the glue itself is no more than its argument checks. Every other
//! call is refused with [`Error::NotBuilt`] until that glue is ported: never a panic, never an
//! approximation.
//!
//! Each call runs on its own [`solidcore::Kernel`] ([`kernel`]), so calls share no state and
//! need no lock.
#![deny(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::unimplemented, clippy::todo, clippy::unreachable)]

use std::path::Path;

use solidcore::SolidError;

mod check;
pub mod export;
pub mod flat;
pub mod history;
pub mod import;
mod measure;
mod ops;
pub mod pcurves;
mod sheets;
mod split;
pub mod ssi;
pub mod sweeps;
pub mod trim;

pub use check::{CheckLevel, CheckReport};
pub use measure::{MEASURE_ACCURACY, MassProperties};
pub use ops::{Blend, BooleanOp, ChamferSecond, Hole, HoleEdge, MAX_POINTS, MAX_UNION_BREPS};

/// Why a call failed: the kinds of `solidcraft_smlib::Error`.
#[derive(Clone, Debug, PartialEq)]
pub enum Error {
    /// SMLib's own status (1001 `SM_ERR`, 1007 `SM_ERR_INVALID_INPUT`, …), as solidcore reports
    /// it.
    Status(i32),
    /// An argument rejected before solidcore saw it.
    Argument,
    /// The call succeeded but returned no body (for a boolean: the result is empty).
    Empty,
    /// The call succeeded, but a solid input gave a result that is not a closed solid.
    NotSolid,
    /// Selection point `point` is farther than the tolerance from every edge (face).
    NoMatch { point: usize },
    /// Selection point `point` is equally near to two edges (faces).
    Ambiguous { point: usize },
    /// No room for the operation: an inward offset leaves none inside the body, or a blend is
    /// larger than the faces it rolls on allow.
    NoRoom,
    /// An offset whose walls would meet along a line.
    WallsTouch,
    /// The kernel cannot build this exactly yet.
    Unsupported,
    /// A blend's edges meet at a vertex of a kind no blend corner handles yet.
    CornerUnsupported,
    /// The corner of corner point `corner` is not set back.
    CornerRefused { corner: usize },
    /// The call succeeded but did nothing.
    NoEffect,
    /// The answer is not exact in the input's geometry.
    NotExact,
    /// A drafted extrusion could not close the corner (or wall) at `at` on its profile, for `why`.
    Corner { why: sweeps::CornerRefusal, at: [f64; 3] },
    /// A failure that has no other kind (SMLib's C API: a caught C++ exception).
    Exception,
    /// The kernel is a different version from the one this crate was written against.
    Version(i32),
    /// Not done by this kernel yet: the call's glue (the fork's C API) is not ported onto
    /// solidcore, or solidcore has not ported the SMLib routine it needs.
    NotBuilt,
    /// A B-rep record the kernel cannot take (see [`import::Unsupported`]).
    Record(import::Unsupported),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::Status(1007) => write!(f, "solidcore rejected the input (SM_ERR_INVALID_INPUT)"),
            Error::Status(s) => write!(f, "solidcore failed (status {s})"),
            Error::Argument => write!(f, "invalid argument for solidcore"),
            Error::Empty => write!(f, "solidcore returned an empty body"),
            Error::NotSolid => write!(f, "solidcore's result is not a closed solid"),
            Error::NoMatch { point } => write!(f, "selection point {} is not near any edge or face", point + 1),
            Error::Ambiguous { point } => write!(f, "selection point {} is equally near to two edges or faces", point + 1),
            Error::NoRoom => write!(f, "there is no room for it: an offset leaves none inside the body, or a blend is larger than its faces allow"),
            Error::WallsTouch => write!(f, "the offset's walls would meet along a line, which leaves no solid there"),
            Error::Unsupported => write!(f, "not supported yet: solidcore cannot build this exactly"),
            Error::CornerUnsupported => write!(f, "not supported yet: a corner of this kind is not supported yet by solidcore's blends"),
            Error::CornerRefused { corner } => write!(f, "not supported yet: solidcore does not set back corner {}", corner + 1),
            Error::NoEffect => write!(f, "solidcore made no change: it reported success, but the operation did nothing"),
            Error::NotExact => write!(f, "solidcore cannot build this exactly in the input's geometry"),
            Error::Corner { why, at } => {
                let at = format!("({:.4}, {:.4}, {:.4})", at[0], at[1], at[2]);
                match why {
                    sweeps::CornerRefusal::Apart => write!(f, "the draft's continued curves do not meet at the corner at {at}"),
                    sweeps::CornerRefusal::WallCloses => write!(f, "the draft closes the wall at {at} before the far end"),
                    sweeps::CornerRefusal::NotExact => {
                        write!(f, "the draft cannot continue the walls exactly: the curve from {at} is neither a line nor an arc")
                    }
                }
            }
            Error::Exception => write!(f, "solidcore failed internally"),
            Error::Version(v) => write!(f, "solidcore version {v} is not the one this build expects"),
            Error::NotBuilt => write!(f, "not supported yet: SolidCraft's solidcore kernel does not do this yet"),
            Error::Record(u) => write!(f, "{u}"),
        }
    }
}

impl std::error::Error for Error {}

pub type Result<T> = std::result::Result<T, Error>;

/// The refusal of a call this kernel does not make yet ([`Error::NotBuilt`]).
fn refused<T>() -> Result<T> {
    Err(Error::NotBuilt)
}

/// `SM_ERR`, SMLib's generic failure.
const SM_ERR: i32 = solidcore::error::SM_ERR;

/// solidcore's error as this crate's kind, the one SMLib's C API reports for the same status
/// (`smlib_shim::status`). The full error goes to the log first (target `solidcore`): the kinds
/// keep only the status.
fn error(e: SolidError) -> Error {
    log::warn!(target: "solidcore", "{e}");
    match e {
        SolidError::Unported { .. } => Error::NotBuilt,
        SolidError::Unsupported { .. } => Error::Unsupported,
        SolidError::InvalidInput { .. } => Error::Status(solidcore::error::SM_ERR_INVALID_INPUT),
        e => match e.sm_status() {
            Some(solidcore::error::SM_ERR_OFFSET_WALLS_TOUCH) => Error::WallsTouch,
            Some(s) => Error::Status(s),
            // A degenerate input, a solver that did not converge, a port invariant: SMLib's
            // routines report these as `SM_ERR`.
            None => Error::Status(SM_ERR),
        },
    }
}

/// A new kernel for one call: SMLib 1.0.1's tolerances and defaults (`SmApiGetOrCreateContext`),
/// with the opt-in fixes on whose SolidCraft fork commit is in the SMLib that `solidcraft-smlib`
/// pins (`crates/smlib/SMLIB_REV`), so both kernels run the same SMLib. A fix solidcore names by
/// a fork commit outside that revision, or by no fork commit, stays off.
fn kernel() -> solidcore::Kernel {
    let mut k = solidcore::Kernel::default();
    let f = &mut k.fixes;
    f.fillet_refuse_overlapping_bands = true;
    f.fillet_vertex_error_diagnostics = true;
    f.fillet_corner_ball_center_from_axes = true;
    f.fillet_torus_span_from_start_section = true;
    f.corner_blend_coons_patch_tangent_to_sides = true;
    f.fillet_patch_keeps_fillet_edges = true;
    f.fillet_knife_edge_rounds = true;
    f.offset_extends_cones_as_themselves = true;
    f.edge_check_orders_tangent_faces = true;
    f.torus_check_reads_v_interval = true;
    f.torus_check_spindle_apex = true;
    f.revolve_arc_touching_axis_is_horn_torus = true;
    f.taper_wall_closes_at_apex = true;
    f.revolve_holed_profile_void = true;
    f.local_operation_box_holds_new_surfaces = true;
    f.local_operation_refuses_loose_tangent_edges = true;
    f.revolution_round_offset_collapses = true;
    f.uv_trim_joins_drop_pieces = true;
    f.boolean2d_plane_covers_far_faces = true;
    f.boolean_splits_definition_across_split = true;
    f.seam_point_solutions_are_one = true;
    f.torus_major_radius_from_origin = true;
    f.found_torus_sphere_keep_tested_nurb = true;
    f.precise_properties_body_length_scale = true;
    f.precise_properties_split_at_surface_knots = true;
    f.brep_area_between_solid_and_void = true;
    f.fillet_nx1_opposite_continuation_side_edge = true;
    f.fillet_3x2_mixed_concave_side_edge = true;
    f.offset_analytic_mate_edges = true;
    f.shell_edges_on_analytics_exact = true;
    f.ssi_no_reverse_of_empty_trace = true;
    f.ssi_pole_point_keeps_its_tangent = true;
    f.boolean_skips_faces_apart = true;
    k
}

/// Nothing to set up: every call makes its own [`kernel`]. Kept for the API.
pub fn init() -> Result<()> {
    Ok(())
}

/// Where failing calls leave their inputs for replay: SMLib's C API writes them (`smlib_set_dump_dir`).
/// solidcore keeps no such replay files, so a folder is taken and nothing is written to it: no
/// answer depends on it.
pub fn set_dump_dir(dir: Option<&Path>) -> Result<()> {
    let _ = dir;
    Ok(())
}

/// The call running right now, for crash reports: solidcore keeps no such record (calls are
/// plain Rust, and a panic carries its own location), so none.
pub fn current_call() -> Option<String> {
    None
}

/// A solidcore B-rep. Immutable: operations return new bodies.
pub struct Brep(solidcore::Brep);

impl std::fmt::Debug for Brep {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("solidcore::Brep")
    }
}

fn finite(v: &[f64]) -> bool {
    v.iter().all(|x| x.is_finite())
}

fn point(p: [f64; 3]) -> solidcore::Point3 {
    solidcore::Point3::new(p[0], p[1], p[2])
}

impl Brep {
    /// An axis-aligned box from its lower corner, with positive edge lengths (mm).
    pub fn make_box(origin: [f64; 3], size: [f64; 3]) -> Result<Brep> {
        if !finite(&origin) || !finite(&size) || !size.iter().all(|s| *s > 0.0) {
            return Err(Error::Argument);
        }
        let mut k = kernel();
        solidcore::primitive::api_primitives::create_box(&mut k, point(origin), size[0], size[1], size[2]).map(Brep).map_err(error)
    }

    /// Reads a native SMLib `.smb` file: not bound yet (solidcore reads one through its topology
    /// naming environment, which this crate does not set up yet).
    pub fn read(path: &Path, ascii: bool) -> Result<Brep> {
        let _ = (path, ascii);
        refused()
    }

    /// Writes this body to a native SMLib `.smb` file: not bound yet, as [`Brep::read`].
    pub fn write(&self, path: &Path, ascii: bool) -> Result<()> {
        let _ = (path, ascii);
        refused()
    }

    /// A deep copy: every curve and surface copied with the topology that owns it
    /// (`Brep::clone_standing_for_self`), as `SmApiBrepCopy` gives.
    pub fn try_clone(&self) -> Result<Brep> {
        self.0.clone_standing_for_self().map(Brep).map_err(error)
    }

    /// This body moved by the affine `matrix`: not bound yet (`SmBrep::Transform` edits a copy
    /// that owns its geometry, which needs solidcore's copy with its naming environment).
    pub fn transformed(&self, matrix: &[[f64; 4]; 4]) -> Result<Brep> {
        let _ = matrix;
        refused()
    }

    /// The enclosed volume (mm³) at `relative_accuracy` (`SmApiBrepComputeVolume`, which clamps
    /// it to [1e-8, 1e-1]), integrated on a copy: solidcore's integrator makes the caches it needs
    /// on the body it measures.
    pub fn volume(&self, relative_accuracy: f64) -> Result<f64> {
        if !relative_accuracy.is_finite() {
            return Err(Error::Argument);
        }
        let k = kernel();
        let mut b = self.0.clone_standing_for_self().map_err(error)?;
        let store = solidcore::geom::attribute::AttributeStore::default();
        let attributes =
            solidcore::mass::integrator::properties::PropertyAttributes { store: &store, key_of: &|_| solidcore::geom::attribute::ObjectKey(0) };
        solidcore::mass::integrator::properties::brep_volume(&k, &mut b, relative_accuracy, &attributes).map_err(error)
    }

    /// Face, edge and vertex counts.
    pub fn counts(&self) -> Result<[usize; 3]> {
        Ok([self.0.face_list().len(), self.0.edge_list().len(), self.0.vertex_list().len()])
    }
}

#[cfg(test)]
mod tests;
