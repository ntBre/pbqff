use std::rc::Rc;

use intder::Intder;
use psqs::{
    program::{mopac::Mopac, Job},
    queue::{local::LocalQueue, Queue},
};
use symm::Molecule;
use taylor::Taylor;

use crate::{config::Config, *};
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

    // load the initial intder
    let mut intder = Intder::load_file("testfiles/intder.in");
    let nsic = intder.symmetry_internals.len();
    // generate a displacement for each SIC
    let mut disps = Vec::new();
    for i in 0..nsic {
        let mut disp = vec![0.0; nsic];
        disp[i] = TAYLOR_DISP_SIZE;
        disps.push(disp);
    }
    intder.geom = intder::geom::Geom::from(mol);
    intder.disps = disps;
    // convert them to Cartesian coordinates
    let disps = intder.convert_disps();
    // convert displacements -> symm::Molecules and determine irrep
    let mut irreps = Vec::new();
    for (i, disp) in disps.iter().enumerate() {
        let m = Molecule::from_slices(atomic_numbers.clone(), disp.as_slice());
        irreps.push((i, m.irrep(&pg)));
    }
    // sort by irrep symmetry
    irreps.sort_by_key(|k| k.1);
    // generate checks
    let checks = make_taylor_checks(irreps, &pg);
    // run taylor.py to get fcs and disps
    let taylor = Taylor::new(5, nsic, checks.0, checks.1);
    let taylor_disps = taylor.disps();
    intder.disps = disp_to_intder(&taylor_disps);

    // these are the displacements that go in file07, but I'll use them from
    // memory to build the jobs
    let _file07 = intder.convert_disps();

    // TODO build and run the points using psqs

    // TODO use these after running the points
    let energies = vec![];
    let _anpass = taylor_to_anpass(&taylor, &taylor_disps, &energies);

    // TODO intder_geom
    // TODO intder freqs
    // TODO spectro
}
