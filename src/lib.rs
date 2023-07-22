pub mod config;
pub mod coord_type;

#[cfg(test)]
mod tests;

use std::path::Path;

pub use intder::Intder;
use psqs::{
    geom::Geom,
    program::{Job, Program, ProgramError, ProgramResult, Template},
    queue::{Check, Queue},
};
use serde::{Deserialize, Serialize};
pub use spectro::{Output, Spectro};

/// clean up from a previous run, emitting warnings on failures
pub fn cleanup(dir: impl AsRef<Path>) {
    for d in ["opt", "pts", "freqs"] {
        let d = dir.as_ref().join(d);
        std::fs::remove_dir_all(&d)
            .unwrap_or_else(|e| eprintln!("failed to remove '{d:?}' with {e}"));
    }
}

pub fn make_check(check_int: usize, dir: impl AsRef<Path>) -> Check {
    if check_int == 0 {
        Check::None
    } else {
        Check::Some {
            check_int,
            check_dir: dir.as_ref().to_str().unwrap().to_owned(),
        }
    }
}

/// Optimize the geometry in `geom` using `template` and the provided `queue`.
/// Creates the `opt` directory if it doesn't already exist.
pub fn optimize<
    Q: Queue<P> + Sync,
    P: Program + Clone + Send + Sync + Serialize + for<'a> Deserialize<'a>,
>(
    dir: impl AsRef<Path>,
    queue: &Q,
    geom: Geom,
    template: Template,
    charge: isize,
) -> Result<ProgramResult, ProgramError> {
    let opt_dir = dir.as_ref().join("opt");
    let _ = std::fs::create_dir(&opt_dir);
    let opt_file = opt_dir.join("opt").to_str().unwrap().to_owned();
    let opt = Job::new(P::new(opt_file, template, charge, geom), 0);
    let mut res = vec![Default::default(); 1];
    let time = queue.energize(opt_dir.to_str().unwrap(), [opt], &mut res)?;
    eprintln!("total optimize time: {time:.1} sec");
    Ok(res.pop().unwrap())
}

/// Like [optimize] but runs a single-point energy instead of an optimization.
/// Also creates the `opt` directory, if it doesn't already exist.
pub fn ref_energy<
    Q: Queue<P> + Sync,
    P: Program + Clone + Send + Sync + Serialize + for<'a> Deserialize<'a>,
>(
    queue: &Q,
    geom: Geom,
    template: Template,
    charge: isize,
) -> f64 {
    let _ = std::fs::create_dir("opt");
    let opt =
        Job::new(P::new("opt/ref".to_string(), template, charge, geom), 0);
    let mut res = vec![0.0; 1];
    let time = queue
        .drain("opt", [opt], &mut res, Check::None)
        .expect("reference energy failed");
    eprintln!("total ref time: {time:.1} sec");
    res.pop().unwrap()
}

// TODO I might have a QFF trait in this package that I can implement on
// psqs::Programs to extend them; pass the common stuff into run and then make
// run a method on this trait
