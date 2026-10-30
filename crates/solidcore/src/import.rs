//! Schema 4 records (`solidcraft_geom::brep::Brep`) into the kernel, as `solidcraft_smlib::import`
//! takes them. There they go as SMLib's C API's flat arrays ([`crate::flat`]), which the fork's
//! glue builds into an `SmBrep`; that bridge is not ported onto solidcore, so a record is not
//! taken yet: [`crate::Brep::from_record`] is refused, and so is laying a record out as the
//! arrays, which only that bridge reads.

use std::collections::BTreeMap;

use solidcraft_geom::brep as b;

/// Why a record cannot go to the kernel.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Unsupported {
    /// The record failed `Brep::validate`.
    Invalid(&'static str),
    /// A part of schema 4 this mapping does not carry yet.
    NotYet(&'static str),
}

impl std::fmt::Display for Unsupported {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Unsupported::Invalid(w) => write!(f, "invalid B-rep record: {w}"),
            Unsupported::NotYet(w) => write!(f, "not supported yet: {w} in solidcore"),
        }
    }
}

type R<T> = std::result::Result<T, Unsupported>;

/// The flat arrays of one B-rep, by array id.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Arrays {
    pub ints: BTreeMap<i32, Vec<i64>>,
    pub reals: BTreeMap<i32, Vec<f64>>,
}

/// What this kernel answers for a record: the flat bridge is not ported.
const NOT_PORTED: Unsupported = Unsupported::NotYet("the flat B-rep arrays (SMLib's C API bridge)");

/// The record as flat arrays: not ported (see the module).
pub fn arrays(brep: &b::Brep, trusted: bool) -> R<Arrays> {
    let _ = (brep, trusted);
    Err(NOT_PORTED)
}

/// [`arrays`] with every surface as the record holds it: not ported.
pub fn arrays_keeping_splines(brep: &b::Brep, trusted: bool) -> R<Arrays> {
    let _ = (brep, trusted);
    Err(NOT_PORTED)
}

impl crate::Brep {
    /// The kernel's B-rep of a schema 4 record: not ported (see the module).
    pub fn from_record(brep: &b::Brep, trusted: bool) -> crate::Result<crate::Brep> {
        let _ = (brep, trusted);
        crate::refused()
    }

    /// [`crate::Brep::from_record`] through [`arrays_keeping_splines`]: not ported.
    pub fn from_record_keeping_splines(brep: &b::Brep, trusted: bool) -> crate::Result<crate::Brep> {
        let _ = (brep, trusted);
        crate::refused()
    }
}
