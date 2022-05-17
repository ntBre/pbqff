use std::rc::Rc;

use intder::Intder;
use psqs::{
    geom::Geom,
    program::{mopac::Mopac, Job},
    queue::{local::LocalQueue, Queue},
};
use spectro::Spectro;
use symm::Molecule;
use taylor::Taylor;

use crate::{config::Config, summary::Summary, *};

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

    let mol = {
        let mut mol = Molecule::new(geom.xyz().unwrap().to_vec());
        mol.normalize();
        mol
    };
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
    intder.geom.to_bohr();
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
    let file07 = intder.convert_disps();

    // TODO build and run the points using psqs
    std::fs::create_dir("pts").unwrap();
    let mut geoms = Vec::with_capacity(file07.len());
    for geom in file07 {
        // this is a bit unsightly, but I also don't want to duplicate the
        // `from_slices` code in psqs
        geoms.push(Rc::new(Geom::from(Molecule::from_slices(
            atomic_numbers.clone(),
            geom.as_slice(),
        ))));
    }
    // TODO switch on Program type eventually

    // TODO take charge from input or template, or don't write charge if
    // self.charge=0 in mopac.write_input
    const CHARGE: isize = 0;
    let mut jobs =
        Mopac::build_jobs(&geoms, None, "pts", 0, 1.0, 0, CHARGE, &MOPAC_TMPL);
    let mut energies = vec![0.0; jobs.len()];
    LocalQueue {
        dir: "pts".to_string(),
    }
    .drain(&mut jobs, &mut energies);
    let min = energies.iter().cloned().reduce(f64::min).unwrap();
    for energy in energies.iter_mut() {
        *energy -= min;
    }
    std::fs::remove_dir_all("pts").unwrap();

    // run anpass
    let anpass = taylor_to_anpass(&taylor, &taylor_disps, &energies);
    let (fcs, long_line) = &anpass.run();

    // intder_geom
    intder.disps = vec![long_line.disp.as_slice().to_vec()];
    let refit_geom = intder.convert_disps();
    let mol =
        Molecule::from_slices(atomic_numbers.clone(), refit_geom[0].as_slice());
    intder.geom = intder::geom::Geom::from(mol.clone());

    // intder freqs
    for fc in fcs {
        // skip zeroth and first derivatives
        if (fc.1, fc.2, fc.3) != (0, 0, 0) {
            intder.add_fc(vec![fc.0, fc.1, fc.2, fc.3], fc.4);
        }
    }

    let (f2, f3, f4) = intder.convert_fcs();
    let _ = std::fs::create_dir("freqs");
    Intder::dump_fcs("freqs", &f2, &f3, &f4);

    // spectro
    let spectro = Spectro {
        geom: mol,
        ..Spectro::load("testfiles/spectro.in")
    };
    spectro.write("freqs/spectro.in").unwrap();

    // run gspectro
    std::process::Command::new(
        "/home/brent/Projects/chemutils/spectro/spectro/spectro",
    )
    .arg("-cmd=/home/brent/Projects/pbqff/bin/spectro")
    .arg("freqs/spectro.in")
    .output()
    .unwrap();

    let summ = Summary::new("freqs/spectro2.out");
    println!("{:?}", summ.fund);
}
