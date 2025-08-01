//! Coordinate types for QFF displacements.

use std::{io::Write, path::Path};

use psqs::{program::Program, queue::Queue};
use serde::{Deserialize, Serialize};
use spectro::{Output, Spectro};

use crate::{config::Config, die};

pub use cart::{Cart, CartGeom, Derivative, FirstPart, Nderiv};
pub use sic::Sic;
pub(crate) use sic::{FreqError, make_rel, write_file};
pub mod cart;
pub mod findiff;
pub mod fitted;
pub mod normal;
pub mod sic;

/// default spectro header to write to the input files for use by the fortran
/// version
const SPECTRO_HEADER: [isize; 30] = [
    99, 1, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
];

/// the name for the checkpoint file written here
const CHK_NAME: &str = "res.chk";

const PTS_DIR: &str = "pts";

/// The very high-level description of a coordinate type for running QFFs.
///
/// Given a [Queue], an output destination, and a [Config], return the final
/// [Spectro] and [Output]
pub trait CoordType<
    W: Write,
    Q: Queue<P>,
    P: Program + Clone + Send + Sync + Serialize + for<'a> Deserialize<'a>,
>
{
    /// All of the data necessary to resume from a checkpoint.
    type Resume: Load;

    /// Run a full QFF on `queue`, taking the configuration from `config`, and
    /// return the final [Spectro] and [Output] structs. Log any output to `w`.
    fn run(
        self,
        dir: impl AsRef<Path>,
        w: &mut W,
        queue: &Q,
        config: &Config,
    ) -> (Spectro, Output);

    /// Like [Self::run], but load a checkpoint from `resume`.
    fn resume(
        self,
        dir: impl AsRef<Path>,
        w: &mut W,
        queue: &Q,
        config: &Config,
        resume: Self::Resume,
    ) -> (Spectro, Output);
}

/// Read and write checkpoints for [CoordType::resume].
pub trait Load: Sized + Serialize + for<'a> Deserialize<'a> {
    fn load(p: impl AsRef<Path>) -> Self {
        let Ok(f) = std::fs::File::open(&p) else {
            die!("failed to load checkpoint from `{}`", p.as_ref().display());
        };
        serde_json::from_reader(f).unwrap()
    }

    fn dump(&self, p: impl AsRef<Path>) {
        let f = std::fs::File::create(p).unwrap();
        serde_json::to_writer(f, self).unwrap()
    }
}
