#![allow(unused)]

use psqs::{
    program::mopac::Mopac,
    queue::{local::LocalQueue, Queue},
};
use spectro::Spectro;
use summarize::Summary;

use crate::{config::Config, optimize, MOPAC_TMPL};

use super::CoordType;

pub struct Cart;

impl CoordType for Cart {
    fn run(&self, config: &Config, spectro: &Spectro) -> Summary {
        let geom = if config.optimize {
            optimize(config.geometry.clone())
        } else {
            todo!();
        };

        let mut jobs = Mopac::build_jobs(
            todo!(), // &geoms
            None,
            "pts",
            0,
            1.0,
            0,
            config.charge,
            &MOPAC_TMPL,
        );

        // this is not going to be jobs.len(). need a different way to determine
        // the length. since I'm using one vector, the first part needs to be
        // fc2, the second part fc3 and the fourth fc4. calculate these sizes
        // with formula from intder.
        let mut fcs = vec![0.0; jobs.len()];
        LocalQueue {
            dir: "pts".to_string(),
        }
        .drain(&mut jobs, &mut fcs);

	// partition the fcs and call freqs

        todo!()
    }
}
