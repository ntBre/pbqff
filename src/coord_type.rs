use std::io::Write;

use psqs::{program::Program, queue::Queue};
use spectro::{Output, Spectro};

use crate::config::Config;

pub use cart::*;
pub use sic::*;
pub mod cart;
pub mod sic;

pub trait CoordType<W: Write, Q: Queue<P>, P: Program + Clone + Send> {
    /// run a full qff, taking the configuration from `config`, the intder
    /// template from `intder`, and the spectro template from `spectro`. Only
    /// the simple internal and symmetry internal coordinates are read from the
    /// intder template. The input options, weights, and curvils are copied from
    /// the spectro template, but the geometry will be updated
    fn run(&self, w: &mut W, queue: &Q, config: &Config) -> (Spectro, Output);
}
