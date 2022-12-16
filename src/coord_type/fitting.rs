//! least-squares fitted coord types

use psqs::geom::Geom;
use spectro::{Output, Spectro};
use std::io::Write;
use symm::{Molecule, PointGroup};
use taylor::{Disps, Taylor};

use super::FreqError;

pub(crate) type AtomicNumbers = Vec<usize>;

type AnpassRes = (Vec<rust_anpass::fc::Fc>, rust_anpass::Bias);
type AnpassError = Box<Result<(Spectro, Output), FreqError>>;

pub trait Fitted {
    type Prep;
    type Error;

    fn prepare_points<W: Write>(
        &mut self,
        mol: &Molecule,
        step_size: f64,
        pg: &PointGroup,
        w: &mut W,
    ) -> Result<Self::Prep, Self::Error>;

    fn generate_pts<W: Write>(
        &mut self,
        w: &mut W,
        mol: &Molecule,
        pg: &PointGroup,
        step_size: f64,
    ) -> Result<(Vec<Geom>, Taylor, Disps, AtomicNumbers), Self::Error>;

    fn anpass<W: Write>(
        &self,
        dir: &str,
        energies: &mut [f64],
        taylor: &Taylor,
        taylor_disps: &taylor::Disps,
        step_size: f64,
        w: &mut W,
    ) -> Result<AnpassRes, AnpassError>;
}
