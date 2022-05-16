pub mod config;

use psqs::program::Template;

/// step size to take in each symmetry internal coordinate to determine its
/// irrep
pub const TAYLOR_DISP_SIZE: f64 = 0.005;

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
    use symm::Molecule;

    use crate::{config::Config, MOPAC_TMPL, TAYLOR_DISP_SIZE};
    use psqs::geom::Geom;

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

        let mol = Molecule::new(geom.xyz().unwrap().to_vec());
        let atomic_numbers = mol.atomic_numbers();
        let pg = mol.point_group();

        let intder = Intder::load_file("testfiles/intder.in");
        let nsic = intder.symmetry_internals.len();
        // generate a displacement for each SIC
        let mut disps = Vec::new();
        for i in 0..nsic {
            let mut disp = vec![0.0; nsic];
            disp[i] = TAYLOR_DISP_SIZE;
            disps.push(disp);
        }
        let symm_intder = Intder {
            geom: intder::geom::Geom::from(mol),
            disps,
            ..intder.clone()
        };
        // convert them to Cartesian coordinates
        let disps = symm_intder.convert_disps();
        // convert displacements -> symm::Molecules and determine irrep
        let mut irreps = Vec::new();
        for disp in disps {
            let m =
                Molecule::from_slices(atomic_numbers.clone(), disp.as_slice());
            irreps.push(m.irrep(&pg));
        }
        dbg!(irreps);
        // TODO sort by irreps
        // TODO generate taylor.py input from that
        // TODO run taylor.py to get fcs and disps
    }
}
