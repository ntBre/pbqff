pub mod config;

use psqs::program::Template;

// TODO take this from a template file
pub const MOPAC_TMPL: Template = Template::from(
    "scfcrt=1.D-21 aux(precision=14) PM6 external=testfiles/params.dat",
);

#[cfg(test)]
mod tests {
    use std::rc::Rc;

    use intder::Intder;
    use psqs::{
        program::{mopac::Mopac, Job},
        queue::{local::LocalQueue, Queue},
    };

    use crate::{config::Config, MOPAC_TMPL};
    use psqs::atom::Geom;

    #[test]
    fn test_full() {
        let config = Config::load("testfiles/test.toml");
        // optimize the geometry
        let geom = if config.optimize {
            std::fs::create_dir("opt").unwrap();
            let opt = Job::new(
                Mopac::new(
                    "opt/opt".to_string(),
                    None,
                    Rc::new(Geom::Zmat(config.geometry)),
                    0,
                    &MOPAC_TMPL,
                ),
                0,
            );
            let submitter = LocalQueue {
                dir: "opt".to_string(),
            };
            submitter.optimize(opt)
        } else {
            todo!();
        };
        std::fs::remove_dir_all("opt").unwrap();
        println!("{}", geom);

        let intder = Intder::load_file("testfiles/intder.in");
        // DONE convert psqs::Geom -> intder::Geom
        let symm_intder = Intder {
            geom: intder::geom::Geom::from(geom),
            ..intder.clone()
        };
        // TODO do displacements for symmetries
        // TODO convert displacements -> symm::Molecules
	// TODO generate taylor.py input from that
	// TODO run taylor.py to get fcs and disps
        dbg!(intder);
    }
}
