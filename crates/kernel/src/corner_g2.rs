//! G2 fillet corners (`corner: G2`): the corner where three rounded edges meet, filled with
//! patches that meet every face around it to the second order, instead of the rolling ball's
//! sphere (or SMLib's Coons patch), which is only tangent to the rounds and keeps its own 1 / r
//! across each seam.
//!
//! The corner is set back, and each round eases into its curvature-continuous section before
//! the setback: SMLib runs each round to its cross-section a setback from the vertex (an
//! iso-curve of the eased round, a quintic meeting each face with no curvature across it) and
//! leaves the corner open ([`crate::blend_smlib::fillet_setback_open`]), so the hole's boundary is
//! the rounds' sections and, on each face between two rounds, SMLib's curve from one round's
//! contact line to the next (a quintic meeting both curvature-continuously), each side split at
//! its middle. The hole is filled with six biquintic patches about a centre point
//! ([`crate::fill::quad_patch`], `solidcraft_geom::fill_quad`), built, not fitted: each holds its
//! two half-sides exactly as its boundary rows and is G2 to the faces beside them and to its two
//! neighbouring patches across their spokes along the whole of every seam, its ends included. The
//! faces are sewn into one solid by SMLib ([`crate::stitch`]).

use solidcraft_geom::fill_quad::QuadFill;

use crate::body::Body;
use crate::{KernelError, Result};

/// A G2 corner's default setback is this times the intercept on its edges of the corner's
/// ellipsoid through the points where its rounds' contact lines cross
/// ([`crate::fillet_g2_corner_setback`]; the intercept is `sqrt 2 r` at a box corner, so the
/// setback 1.77 r), less where that and the ease ([`ease`]) would not fit an edge at the corner.
pub const SETBACK_SCALE: f64 = 1.25;

/// The setback of a G2 corner the ellipsoid cannot read (one with a curved face or edge): this
/// times its round's radius along each edge. It must reach past where the rounds meet without
/// one (the radius, at a square corner), and SMLib's own patch bulges past the faces beyond
/// about twice the radius.
pub const DEFAULT_SETBACK_FACTOR: f64 = 1.5;

/// Each round's ease into the curvature-continuous section before a G2 corner's setback: this
/// times its radius along each edge, at every G2 corner (SMLib's `smlib_blend_setback_open_eased`).
/// The ease makes the hole's conditions agree at its corners: the section meets each planar face
/// with no curvature across, and the face's curve leaves the contact line without turning, so the
/// patches can be G2 to the rounds and the faces all the way to the corners. Over the ease the
/// section's middle moves towards the edge by 0.075 r (from the circle's 0.414 r to the quintic's
/// 0.339 r, along the bisector at a square edge), on the ease law's steepest slope 1.875 / length:
/// the surface's normal turns along the edge by `atan(0.141 / EASE_FACTOR)` there, about 10° at
/// 0.8 (30° at a quarter of the radius, a dark streak where each round meets the corner). A round
/// keeps its own section up to where its ease starts, `setback + ease` from the vertex.
pub const EASE_FACTOR: f64 = 0.8;

/// The rounds' ease before a G2 corner's setback, rounds of `radius` ([`EASE_FACTOR`]).
pub fn ease(radius: f64) -> f64 {
    EASE_FACTOR * radius
}

/// How far (this share of the corner's round's curvature 1 / r, the jump a rolling ball makes)
/// the hole's own data may be from the curvatures the patches' construction asks of it: a
/// section's curvature across its plane where it ends, a face curve's where it leaves a contact
/// line, a round's along its rail beside its section (all zero for an eased setback).
pub const CURVATURE_SHARE: f64 = 5e-3;

/// The patches filling one open corner of `sheet` (a sheet from
/// [`crate::blend_smlib::fillet_setback_open`]): `hole` its edges in order, each with whether the
/// face beside it runs it forwards; `radius` the corner's round's radius, which sets the
/// curvature tolerance its data is read to. The patches are a sheet on the sheet's own edges.
pub(crate) fn corner_patch(sheet: &Body, hole: &[(usize, bool)], radius: f64) -> Result<(Body, QuadFill)> {
    if !(radius.is_finite() && radius > 0.0) {
        return Err(KernelError::Invalid("a G2 corner's radius must be finite and positive".into()));
    }
    crate::fill::quad_patch(sheet, hole, CURVATURE_SHARE / radius)
}

#[cfg(test)]
#[path = "corner_g2_tests.rs"]
mod tests;
