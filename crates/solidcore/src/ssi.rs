//! Surface/surface intersection, as `solidcraft_smlib::ssi`'s. solidcore has SMLib's
//! intersectors (`Surface::global_surface_intersect`, the advanced and the basic intersector);
//! the C API's surfaces (built from an SSI corpus case's numbers, `smlib_ssi.cpp`) and its result
//! records are not ported yet. So no [`Surface`] is made: every constructor is refused, and a
//! `Surface` value cannot exist.

use solidcraft_geom::brep as b;

use crate::{Result, refused};

/// A surface with its UV domain. This kernel makes none yet (see the module).
pub struct Surface {
    never: std::convert::Infallible,
}

impl std::fmt::Debug for Surface {
    fn fmt(&self, _: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.never {}
    }
}

impl Surface {
    /// `origin + u*x + v*y`: not ported.
    pub fn plane(origin: [f64; 3], xaxis: [f64; 3], yaxis: [f64; 3], u: [f64; 2], v: [f64; 2]) -> Result<Surface> {
        let _ = (origin, xaxis, yaxis, u, v);
        refused()
    }

    /// A cylinder: not ported.
    pub fn cylinder(origin: [f64; 3], xaxis: [f64; 3], yaxis: [f64; 3], radius: f64, angle_deg: [f64; 2], height: [f64; 2]) -> Result<Surface> {
        let _ = (origin, xaxis, yaxis, radius, angle_deg, height);
        refused()
    }

    /// A full cone: not ported.
    pub fn cone(apex: [f64; 3], xaxis: [f64; 3], yaxis: [f64; 3], half_angle_deg: f64, height: [f64; 2]) -> Result<Surface> {
        let _ = (apex, xaxis, yaxis, half_angle_deg, height);
        refused()
    }

    /// A sphere: not ported.
    pub fn sphere(center: [f64; 3], xaxis: [f64; 3], yaxis: [f64; 3], radius: f64, lon_deg: [f64; 2], lat_deg: [f64; 2]) -> Result<Surface> {
        let _ = (center, xaxis, yaxis, radius, lon_deg, lat_deg);
        refused()
    }

    /// A ring torus: not ported.
    pub fn torus(
        center: [f64; 3],
        xaxis: [f64; 3],
        yaxis: [f64; 3],
        major_radius: f64,
        minor_radius: f64,
        major_deg: [f64; 2],
        minor_deg: [f64; 2],
    ) -> Result<Surface> {
        let _ = (center, xaxis, yaxis, major_radius, minor_radius, major_deg, minor_deg);
        refused()
    }

    /// A (rational) B-spline surface: not ported.
    pub fn nurbs(
        degree: [usize; 2],
        nu: usize,
        nv: usize,
        knots_u: &[f64],
        knots_v: &[f64],
        points: &[[f64; 3]],
        weights: Option<&[f64]>,
    ) -> Result<Surface> {
        let _ = (degree, nu, nv, knots_u, knots_v, points, weights);
        refused()
    }

    /// A schema 4 record surface over the parameter box `domain`: not ported.
    pub fn from_record(surface: &b::Surface, domain: [[f64; 2]; 2]) -> Result<Surface> {
        let _ = (surface, domain);
        refused()
    }

    /// The UV domain the intersection runs within, as `[u0, u1, v0, v1]`.
    pub fn domain(&self) -> Result<[f64; 4]> {
        match self.never {}
    }

    /// The surface's point at (u, v).
    pub fn evaluate(&self, u: f64, v: f64) -> Result<[f64; 3]> {
        let _ = (u, v);
        match self.never {}
    }
}

/// Which of SMLib's intersectors runs.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Method {
    /// `SmSurface::GlobalSurfaceIntersect`: analytic special cases, else [`Method::Advanced`].
    Global,
    /// `SmAdvSurfaceIntersector`, the general NURBS path.
    Advanced,
    /// `SmSurfaceIntersector`, the basic intersector.
    Basic,
}

/// A curve's type, SMLib's `SmTsectCurveType`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CurveType {
    /// The method reported no type.
    None,
    Unknown,
    /// Degenerate: one point where the surfaces touch.
    Touching,
    /// Transversal.
    Crossing,
    /// The surfaces are tangent along it.
    Tangent,
    /// Where the surfaces butt up against each other.
    Coincident,
    /// A small angle of intersection.
    NearTangent,
    /// Bounds a region where the surfaces coincide.
    RegionBoundary,
    /// A value this crate does not know.
    Other(i32),
}

/// One intersection curve: its 3D curve and, where the method returned them, its curves in each
/// surface's UV.
#[derive(Clone, Debug)]
pub struct SsiCurve {
    pub curve_type: CurveType,
    /// The kernel's deviation estimate (NaN: none).
    pub deviation: f64,
    /// The 3D curve's natural parameter interval.
    pub interval: [f64; 2],
    pub curve: b::Curve,
    /// On surface a, in a's UV.
    pub uv_a: Option<b::Curve2>,
    /// On surface b, in b's UV.
    pub uv_b: Option<b::Curve2>,
}

/// An isolated point: where a [`CurveType::Touching`] curve starts.
#[derive(Clone, Copy, Debug)]
pub struct SsiPoint {
    pub xyz: [f64; 3],
    pub uv_a: [f64; 2],
    pub uv_b: [f64; 2],
}

/// What an intersection found.
#[derive(Clone, Debug, Default)]
pub struct SsiResult {
    pub curves: Vec<SsiCurve>,
    pub points: Vec<SsiPoint>,
}

/// Intersects `a` and `b` within their domains with `method` at approximation tolerance
/// `approx_tol`.
pub fn intersect(a: &Surface, b: &Surface, method: Method, approx_tol: f64) -> Result<SsiResult> {
    let _ = (b, method, approx_tol);
    match a.never {}
}
