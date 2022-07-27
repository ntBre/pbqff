pub mod config;
pub mod coord_type;

#[cfg(test)]
mod tests;

use std::rc::Rc;

pub use intder::Intder;
use psqs::{
    geom::Geom,
    program::{mopac::Mopac, Job, Template, ProgramResult},
    queue::Queue,
};
pub use spectro::Spectro;

pub fn optimize<Q: Queue<Mopac>>(
    queue: &Q,
    geom: Geom,
    template: Template,
    charge: isize,
) -> ProgramResult {
    let _ = std::fs::create_dir("opt");
    let opt = Job::new(
        Mopac::new(
            "opt/opt".to_string(),
            None,
            Rc::new(geom),
            charge,
            template,
        ),
        0,
    );
    let mut res = vec![Default::default(); 1];
    queue.energize(&mut [opt], &mut res);
    res.pop().unwrap()
}

// TODO this is immensely stupid, just take the energy from the original
// optimization, but that's now how I wrote the queue stuff
pub fn ref_energy<Q: Queue<Mopac>>(
    queue: &Q,
    geom: Geom,
    template: Template,
    charge: isize,
) -> f64 {
    let _ = std::fs::create_dir("opt");
    let opt = Job::new(
        Mopac::new(
            "opt/ref".to_string(),
            None,
            Rc::new(geom),
            charge,
            template,
        ),
        0,
    );
    let mut res = vec![0.0; 1];
    queue.drain(&mut [opt], &mut res);
    res.pop().unwrap()
}

// TODO I might have a QFF trait in this package that I can implement on
// psqs::Programs to extend them; pass the common stuff into run and then make
// run a method on this trait
