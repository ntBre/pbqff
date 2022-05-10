pub mod config;

use psqs::program::Template;

// TODO take this from a template file
pub const MOPAC_TMPL: Template = Template::from(
    "scfcrt=1.D-21 aux(precision=14) PM6 external=testfiles/params.dat",
);

#[cfg(test)]
mod tests {
    use std::rc::Rc;

    use psqs::{
        program::{mopac::Mopac, Job},
        queue::{local::LocalQueue, Queue},
    };

    use crate::{config::Config, MOPAC_TMPL};
    use psqs::atom::Geom;

    #[test]
    fn test_full() {
        let config = Config::load("testfiles/test.toml");
        if config.optimize {
            match std::fs::create_dir("opt") {
                Ok(_) => (),
                _ => (),
                // TODO use this in the real version
                // Err(e) => {
                //     if e.kind() == ErrorKind::AlreadyExists {
                //         eprintln!("directory `opt` already exists");
                //         std::process::exit(1);
                //     } else {
                //         panic!();
                //     }
                // }
            };
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
            let got = submitter.optimize(opt);
            dbg!(got);
        }
        std::fs::remove_dir_all("opt").unwrap();
    }
}
