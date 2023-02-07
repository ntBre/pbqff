use std::{io::Write, path::Path};

use psqs::{program::Program, queue::Queue};
use serde::{Deserialize, Serialize};
use spectro::{Output, Spectro};

use crate::config::Config;

#[macro_export]
macro_rules! time {
    ($w:expr, $label:expr, $($s:stmt;)*) => {
        let now = ::std::time::Instant::now();
	$($s)*
        writeln!(
            $w,
            "finished {} after {:.1} sec",
            $label,
            now.elapsed().as_millis() as f64 / 1000.0
        )
        .unwrap();
    };
}

pub use cart::*;
pub use sic::*;
pub mod cart;
pub mod findiff;
pub mod fitting;
pub mod normal;
pub mod sic;

/// default spectro header to write to the input files for use by the fortran
/// version
const SPECTRO_HEADER: [usize; 30] = [
    99, 1, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0,
];

/// the name for the checkpoint file written here
const CHK_NAME: &str = "res.chk";

pub trait CoordType<
    W: Write,
    Q: Queue<P>,
    P: Program + Clone + Send + Sync + Serialize + for<'a> Deserialize<'a>,
>
{
    /// run a full qff, taking the configuration from `config`, the intder
    /// template from `intder`, and the spectro template from `spectro`. Only
    /// the simple internal and symmetry internal coordinates are read from the
    /// intder template. The input options, weights, and curvils are copied from
    /// the spectro template, but the geometry will be updated
    fn run(self, w: &mut W, queue: &Q, config: &Config) -> (Spectro, Output);

    /// contains all of the data necessary to resume from a checkpoint
    type Resume: Load;

    /// resume from a checkpoint and finish the run
    fn resume(
        self,
        w: &mut W,
        queue: &Q,
        config: &Config,
        resume: Self::Resume,
    ) -> (Spectro, Output);
}

pub trait Load: Sized + Serialize + for<'a> Deserialize<'a> {
    fn load(p: impl AsRef<Path>) -> Self {
        let f = std::fs::File::open(p).unwrap();
        serde_json::from_reader(f).unwrap()
    }

    fn dump(&self, p: impl AsRef<Path>) {
        let f = std::fs::File::create(p).unwrap();
        serde_json::to_writer(f, self).unwrap()
    }
}
