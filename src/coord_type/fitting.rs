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

pub trait Fitted {
    /// The error returned by a failing [Self::generate_points].
    type Error;

    fn generate_pts<W: Write>(
        &mut self,
        w: &mut W,
        mol: &Molecule,
        pg: &PointGroup,
        step_size: f64,
    ) -> Result<(Vec<Geom>, Taylor, AtomicNumbers), Self::Error>;

    fn anpass<W: Write>(
        &self,
        dir: Option<&str>,
        energies: &mut [f64],
        taylor: &Taylor,
        step_size: f64,
        w: &mut W,
    ) -> Result<AnpassRes, AnpassError>;
}
