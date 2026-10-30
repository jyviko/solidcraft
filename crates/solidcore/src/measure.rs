//! Exact measures, as `solidcraft_smlib`'s: the exact box. Mass properties, area and the
//! edge/face gap are the C API's own sums over SMLib's face integrals (`smlib_measure.cpp`), not
//! ported yet, so they are refused.

use crate::{Brep, Error, Result, error, kernel, refused};

/// The accuracy [`Brep::mass_properties`] and [`Brep::area`] use unless told otherwise: the
/// finest SMLib allows.
pub const MEASURE_ACCURACY: f64 = 1e-8;

/// Mass properties of a body at unit density (mm, mm², mm³, mm⁵). For a density ρ the
/// mass is ρ·`volume` and the inertia scales by ρ.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MassProperties {
    pub area: f64,
    pub volume: f64,
    pub centroid: [f64; 3],
    /// Ixx, Iyy, Izz about axes through the centroid parallel to the world axes:
    /// ∫(y²+z²), ∫(z²+x²), ∫(x²+y²) dV.
    pub moments: [f64; 3],
    /// Iyz, Izx, Ixy about the centroid: ∫yz, ∫zx, ∫xy dV relative to it.
    pub products: [f64; 3],
}

impl MassProperties {
    /// The inertia tensor about the centroid (row-major), whose off-diagonal terms are the
    /// negated products.
    pub fn inertia_tensor(&self) -> [[f64; 3]; 3] {
        let [ixx, iyy, izz] = self.moments;
        let [iyz, izx, ixy] = self.products;
        [[ixx, -ixy, -izx], [-ixy, iyy, -iyz], [-izx, -iyz, izz]]
    }
}

impl Brep {
    /// Mass properties at `relative_accuracy`: not ported (see the module).
    pub fn mass_properties(&self, relative_accuracy: f64) -> Result<MassProperties> {
        let _ = relative_accuracy;
        refused()
    }

    /// The exact axis-aligned box, least corner then greatest: the extrema of every vertex, edge
    /// and face on its exact curve or surface (`SmBrep::CalculateTightBoundingBox`).
    pub fn bounds(&self) -> Result<([f64; 3], [f64; 3])> {
        let b = solidcore::query::api_queries::brep_bounding_box(&kernel(), &self.0, true).map_err(error)?;
        let lo = [b.min_at(0), b.min_at(1), b.min_at(2)];
        let hi = [b.max_at(0), b.max_at(1), b.max_at(2)];
        // As the C API checks it: a box of numbers, least corner first.
        if !crate::finite(&lo) || !crate::finite(&hi) || (0..3).any(|i| lo.get(i) > hi.get(i)) {
            return Err(Error::Status(solidcore::error::SM_ERR_INVALID_INPUT));
        }
        Ok((lo, hi))
    }

    /// How far edge `edge` lies from the surface of face `face`: not ported (see the module).
    pub fn edge_face_distance(&self, edge: usize, face: usize) -> Result<f64> {
        let _ = (edge, face);
        refused()
    }

    /// Total area of every face: not ported (see the module).
    pub fn area(&self, relative_accuracy: f64) -> Result<f64> {
        let _ = relative_accuracy;
        refused()
    }
}
