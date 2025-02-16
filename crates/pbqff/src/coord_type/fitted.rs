//! Compute force constants using [ordinary
//! least-squares](https://en.wikipedia.org/wiki/Ordinary_least_squares)
//! fitting.

use std::{io::Write, path::Path};

use psqs::geom::Geom;
use spectro::{Output, Spectro};
use symm::{Molecule, PointGroup};
use taylor::Taylor;

use super::FreqError;

pub(crate) type AtomicNumbers = Vec<usize>;
type AnpassRes = (Vec<anpass::fc::Fc>, anpass::Bias);
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
        dir: impl AsRef<Path>,
        w: &mut W,
        mol: &Molecule,
        pg: &PointGroup,
        step_size: f64,
    ) -> Result<(Vec<Geom>, Taylor, AtomicNumbers), Self::Error>;

    /// Perform the fitting of the potential energy surface. Return the final
    /// force constants ([anpass::fc::Fc]) and the corresponding
    /// [anpass::Bias].
    fn anpass<W: Write>(
        &self,
        dir: Option<impl AsRef<Path>>,
        energies: &mut [f64],
        taylor: &Taylor,
        step_size: f64,
        w: &mut W,
    ) -> Result<AnpassRes, AnpassError>;
}
