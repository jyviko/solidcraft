//! Sweeps, as `solidcraft_smlib::sweeps`'s: extrusion (straight or drafted), revolution, sweep
//! along a path, loft, screw sweep and the sweep through stations. solidcore has SMLib's sweep
//! routines (`sweep::api_brep`, including the fork's drafted corners and sweep placements); the
//! C API's glue (its curve lists made into SMLib curves, the capping, the result checks,
//! `smlib_sweeps.cpp`) is not ported yet, so every sweep is refused.

use crate::{Brep, Result, refused};

/// A 3D curve of a profile, path or section.
#[derive(Clone, Debug, PartialEq)]
pub enum Curve {
    Line {
        start: [f64; 3],
        end: [f64; 3],
    },
    /// Counter-clockwise about the unit `axis` from angle `start` to `end` (radians,
    /// 0 < end - start <= 2 pi), angle 0 along the unit `reference` (square to `axis`).
    Arc {
        center: [f64; 3],
        axis: [f64; 3],
        reference: [f64; 3],
        radius: f64,
        start: f64,
        end: f64,
    },
    /// Control points (x, y, z, weight > 0) and the clamped knot vector (points + degree + 1
    /// knots).
    Nurbs {
        degree: usize,
        points: Vec<[f64; 4]>,
        knots: Vec<f64>,
    },
}

/// How a drafted [`Brep::extrude`] closes the gap its offset opens at a sharp corner of the
/// profile.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum DraftCorner {
    /// Each wall continues as itself until it meets its neighbour.
    #[default]
    Curvature,
    /// Each offset curve continues straight along its end tangent.
    Tangent,
    /// An arc about the corner: a cone about it.
    Round,
}

/// Why [`DraftCorner::Curvature`] could not close a corner ([`crate::Error::Corner`]).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CornerRefusal {
    /// The continued curves do not meet.
    Apart,
    /// A wall closes before the far end.
    WallCloses,
    /// A curve there is neither a line nor an arc: its wall has no exact continuation.
    NotExact,
}

/// How the profile is carried along a [`Brep::sweep`] path.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SweepFrame {
    /// It turns with the path, keeping its attitude to the path's tangent.
    Follow,
    /// It only moves along the path, keeping its orientation.
    Translate,
}

/// How the ends of a [`Brep::sweep_stations`] close.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Ends {
    /// Left open: the swept sheet.
    Open,
    /// Capped by the first and last sections: a closed solid.
    Capped,
    /// Joined: the last section lies on the first, the walls' ends glued to their starts.
    Joined,
}

/// An affine map of a [`Brep::sweep_stations`] station: a point `x` goes to `m · x + t`, `m`
/// given by its columns.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Affine {
    pub m: [[f64; 3]; 3],
    pub t: [f64; 3],
}

/// Most stations of one [`Brep::sweep_stations`] (the C API's limit).
pub const MAX_STATIONS: usize = 4096;
/// Most profile curves times stations of one [`Brep::sweep_stations`] (the C API's limit).
pub const MAX_STATION_CURVES: usize = 400_000;

/// Most curves one list may hold (the C API's limit).
pub const MAX_CURVES: usize = 16_384;

impl Brep {
    /// A sheet of planar faces bounded by `loops`: not ported.
    pub fn planar_faces(loops: &[Curve]) -> Result<Brep> {
        let _ = loops;
        refused()
    }

    /// The planar `profile` extruded `distance` along `direction`: not ported.
    pub fn extrude(profile: &[Curve], direction: [f64; 3], distance: f64, draft: f64, corner: DraftCorner, cap: bool) -> Result<Brep> {
        let _ = (profile, direction, distance, draft, corner, cap);
        refused()
    }

    /// The planar `profile` turned `angle` radians about an axis: not ported.
    pub fn revolve(profile: &[Curve], origin: [f64; 3], direction: [f64; 3], angle: f64, cap: bool) -> Result<Brep> {
        let _ = (profile, origin, direction, angle, cap);
        refused()
    }

    /// `profile` swept along `path`: not ported.
    pub fn sweep(profile: &[Curve], path: &[Curve], frame: SweepFrame, cap: bool) -> Result<Brep> {
        let _ = (profile, path, frame, cap);
        refused()
    }

    /// [`Brep::sweep`], with where the sweep placed the profile: not ported.
    pub fn sweep_placed(profile: &[Curve], path: &[Curve], frame: SweepFrame, cap: bool) -> Result<(Brep, [Affine; 3])> {
        let _ = (profile, path, frame, cap);
        refused()
    }

    /// The `curves` projected along `direction` onto the face nearest `face_point`: not ported.
    pub fn project_curves(
        &self,
        face_point: [f64; 3],
        pick_tolerance: f64,
        curves: &[Curve],
        direction: [f64; 3],
        tolerance: f64,
    ) -> Result<Vec<(Curve, [f64; 2])>> {
        let _ = (face_point, pick_tolerance, curves, direction, tolerance);
        refused()
    }

    /// The section of the face nearest `face_point` by a plane: not ported.
    pub fn face_section(
        &self,
        face_point: [f64; 3],
        pick_tolerance: f64,
        plane_point: [f64; 3],
        plane_normal: [f64; 3],
        tolerance: f64,
    ) -> Result<Vec<(Curve, [f64; 2])>> {
        let _ = (face_point, pick_tolerance, plane_point, plane_normal, tolerance);
        refused()
    }

    /// A loft through `sections`: not ported.
    pub fn loft(sections: &[Vec<Curve>], degree: u32, cap: bool) -> Result<Brep> {
        let _ = (sections, degree, cap);
        refused()
    }

    /// The planar `profile` placed by each of `stations` and skinned: not ported.
    pub fn sweep_stations(profile: &[Curve], stations: &[Affine], ends: Ends) -> Result<Brep> {
        let _ = (profile, stations, ends);
        refused()
    }

    /// A loft from the planar `section` to the point `apex`: not ported.
    pub fn loft_to_point(section: &[Curve], apex: [f64; 3]) -> Result<Brep> {
        let _ = (section, apex);
        refused()
    }

    /// `profile` carried by a screw motion: not ported (the screw sweep is the fork's own).
    pub fn screw_sweep(profile: &[Curve], origin: [f64; 3], direction: [f64; 3], pitch: f64, turns: f64, tolerance: f64) -> Result<Brep> {
        let _ = (profile, origin, direction, pitch, turns, tolerance);
        refused()
    }

    /// [`screw_sweep`](Self::screw_sweep) with the profile moved off the axis as it turns: not
    /// ported.
    pub fn screw_sweep_radial(
        profile: &[Curve],
        origin: [f64; 3],
        direction: [f64; 3],
        pitch: f64,
        turns: f64,
        law: &[(f64, f64)],
        tolerance: f64,
    ) -> Result<Brep> {
        let _ = (profile, origin, direction, pitch, turns, law, tolerance);
        refused()
    }
}
