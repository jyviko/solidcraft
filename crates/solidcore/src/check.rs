//! SMLib's validity check (`SmBrep::ValidatePointers`, then `AssertValid`), as
//! `solidcraft_smlib::CheckReport` reports it. solidcore has both routines
//! (`topology::integrity`, `topology::brep::validation`); the C API's report (its text, one line
//! per failed test) is not ported, so [`Brep::check`] is refused.

use crate::{Brep, Result, refused};

/// How many of SMLib's tests to run.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CheckLevel {
    /// The fewest: topology and pointers, no geometry.
    Basic = 0,
    /// Adds the geometry checks (gaps between vertices, edges and faces against tolerances).
    Geometry = 1,
    /// Every test, the slow ones too.
    All = 2,
}

/// What the check found.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CheckReport {
    /// Failed tests; 0 means the body is valid.
    pub failures: usize,
    /// One line per failed test, such as
    /// `SmBrep::Test #82 Tolerance: Vertex tolerance relative to brep tolerance (offset=2, tolerance=0.01) [owner SmVertex]`.
    pub text: String,
}

impl CheckReport {
    pub fn is_valid(&self) -> bool {
        self.failures == 0
    }

    /// The number of each failed test, line by line (`SmBrep::Test #82 …` is 82); a line without
    /// one is left out.
    pub fn failed_tests(&self) -> Vec<u32> {
        self.text
            .lines()
            .filter_map(|l| {
                let rest = l.split_once("Test #")?.1;
                let digits: String = rest.chars().take_while(char::is_ascii_digit).collect();
                digits.parse().ok()
            })
            .collect()
    }
}

impl std::fmt::Display for CheckReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_valid() { f.write_str("valid") } else { write!(f, "{} failed tests:\n{}", self.failures, self.text) }
    }
}

impl Brep {
    /// The validity check at `level`: not ported (see the module).
    pub fn check(&self, level: CheckLevel) -> Result<CheckReport> {
        let _ = level;
        refused()
    }
}
