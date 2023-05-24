//! Compute force constants using [ordinary
//! least-squares](https://en.wikipedia.org/wiki/Ordinary_least_squares)
//! fitting.

use psqs::geom::Geom;
use spectro::{Output, Spectro};
use std::io::Write;
use symm::{Molecule, PointGroup};
use taylor::Taylor;

use super::FreqError;

pub(crate) type AtomicNumbers = Vec<usize>;
type AnpassRes = (Vec<rust_anpass::fc::Fc>, rust_anpass::Bias);
type AnpassError = Box<Result<(Spectro, Output), FreqError>>;

/// Required methods for a least-squares fitted [crate::coord_type::CoordType].
pub trait Fitted {
    /// The error returned by a failing [Self::generate_pts].
    type Error;

    /// Displace `mol` by `step_size` to generate all of the geometries
    /// comprising the QFF. Return the geometries, along with the [Taylor] and
    /// the atomic numbers needed to compute the frequencies.
    fn generate_pts<W: Write>(
        &mut self,
        w: &mut W,
        mol: &Molecule,
        pg: &PointGroup,
        step_size: f64,
    ) -> Result<(Vec<Geom>, Taylor, AtomicNumbers), Self::Error>;

    /// Perform the fitting of the potential energy surface. Return the final
    /// force constants ([rust_anpass::fc::Fc]) and the corresponding
    /// [rust_anpass::Bias].
    fn anpass<W: Write>(
        &self,
        dir: Option<&str>,
        energies: &mut [f64],
        taylor: &Taylor,
        step_size: f64,
        w: &mut W,
    ) -> Result<AnpassRes, AnpassError>;
}
