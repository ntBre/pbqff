use std::io::Write;

use psqs::{program::mopac::Mopac, queue::Queue};
use spectro::Spectro;
use summarize::Summary;

use crate::config::Config;

pub use cart::*;
pub use sic::*;
pub mod cart;
pub mod sic;

pub trait CoordType<W: Write, Q: Queue<Mopac>> {
    /// run a full qff, taking the configuration from `config`, the intder
    /// template from `intder`, and the spectro template from `spectro`. Only
    /// the simple internal and symmetry internal coordinates are read from the
    /// intder template. The input options, weights, and curvils are copied from
    /// the spectro template, but the geometry will be updated
    fn run(
        &self,
        w: &mut W,
        queue: &Q,
        config: &Config,
        spectro: &Spectro,
    ) -> Summary;
}
