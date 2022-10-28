#![feature(test)]
pub mod config;
pub mod coord_type;

#[cfg(test)]
mod tests;

pub use intder::Intder;
use psqs::{
    geom::Geom,
    program::{Job, Program, ProgramError, ProgramResult, Template},
    queue::Queue,
};
pub use spectro::Spectro;

/// clean up from a previous run, emitting warnings on failures
pub fn cleanup() {
    std::fs::remove_dir_all("opt")
        .unwrap_or_else(|e| eprintln!("failed to remove 'opt' with {e}"));
    std::fs::remove_dir_all("pts")
        .unwrap_or_else(|e| eprintln!("failed to remove 'pts' with {e}"));
    std::fs::remove_dir_all("freqs")
        .unwrap_or_else(|e| eprintln!("failed to remove 'freqs' with {e}"));
}


/// Some if the optimization succeeds, None otherwise
pub fn optimize<Q: Queue<P>, P: Program + Clone + Send>(
    queue: &Q,
    geom: Geom,
    template: Template,
    charge: isize,
) -> Result<ProgramResult, ProgramError> {
    let _ = std::fs::create_dir("opt");
    let opt =
        Job::new(P::new("opt/opt".to_string(), template, charge, geom), 0);
    let mut res = vec![Default::default(); 1];
    queue.energize("opt", &mut [opt], &mut res)?;
    Ok(res.pop().unwrap())
}

pub fn ref_energy<Q: Queue<P>, P: Program + Clone + Send>(
    queue: &Q,
    geom: Geom,
    template: Template,
    charge: isize,
) -> f64 {
    let _ = std::fs::create_dir("opt");
    let opt =
        Job::new(P::new("opt/ref".to_string(), template, charge, geom), 0);
    let mut res = vec![0.0; 1];
    queue
        .drain("opt", &mut [opt], &mut res)
        .expect("reference energy failed");
    res.pop().unwrap()
}

// TODO I might have a QFF trait in this package that I can implement on
// psqs::Programs to extend them; pass the common stuff into run and then make
// run a method on this trait
