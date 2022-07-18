pub mod config;
pub mod coord_type;

#[cfg(test)]
mod tests;

use std::rc::Rc;

pub use intder::Intder;
use psqs::{
    geom::Geom,
    program::{mopac::Mopac, Job, Template},
    queue::Queue,
};
pub use spectro::Spectro;

// TODO take this from a template file
pub(crate) const MOPAC_TMPL: Template = Template::from(
    "scfcrt=1.D-21 aux(precision=14) PM6 external=testfiles/params.dat",
);

pub fn optimize<Q: Queue<Mopac>>(queue: &Q, geom: Geom) -> Geom {
    let _ = std::fs::create_dir("opt");
    let opt = Job::new(
        Mopac::new("opt/opt".to_string(), None, Rc::new(geom), 0, &MOPAC_TMPL),
        0,
    );
    let mut res = vec![Geom::default(); 1];
    queue.optimize(&mut [opt], &mut res);
    res.pop().unwrap()
}

// TODO this is immensely stupid, just take the energy from the original
// optimization, but that's now how I wrote the queue stuff
pub fn ref_energy<Q: Queue<Mopac>>(queue: &Q, geom: Geom) -> f64 {
    // only difference is the addition of `1SCF` to compute a single-point
    // instead of the default optimization
    const MOPAC_TMPL: Template = Template::from(
        "scfcrt=1.D-21 aux(precision=14) PM6 external=testfiles/params.dat",
    );
    let _ = std::fs::create_dir("opt");
    let opt = Job::new(
        Mopac::new("opt/ref".to_string(), None, Rc::new(geom), 0, &MOPAC_TMPL),
        0,
    );
    let mut res = vec![0.0; 1];
    queue.drain(&mut [opt], &mut res);
    res.pop().unwrap()
}

// TODO I might have a QFF trait in this package that I can implement on
// psqs::Programs to extend them; pass the common stuff into run and then make
// run a method on this trait
