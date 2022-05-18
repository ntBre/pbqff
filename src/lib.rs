pub mod config;
pub mod summary;

#[cfg(test)]
mod tests;

use std::rc::Rc;

use anpass::Anpass;
pub use intder::Intder;
use nalgebra as na;
use psqs::{
    geom::Geom,
    program::{mopac::Mopac, Job, Template},
    queue::{local::LocalQueue, Queue},
};
use rust_anpass as anpass;
pub use spectro::Spectro;
use symm::{Irrep, Molecule, PointGroup};
use taylor::{Checks, Taylor};

use crate::{config::Config, summary::Summary};

/// step size to take in each symmetry internal coordinate to determine its
/// irrep
pub const TAYLOR_DISP_SIZE: f64 = 0.005;

// TODO take this from a template file
pub const MOPAC_TMPL: Template = Template::from(
    "A0 scfcrt=1.D-21 aux(precision=14) PM6 external=testfiles/params.dat",
);

/// generate the Taylor series mod and equivalence checks from `irreps` in `pg`
fn make_taylor_checks(
    irreps: Vec<(usize, Irrep)>,
    pg: &PointGroup,
) -> (Option<Checks>, Option<Checks>) {
    use symm::Irrep::*;
    use symm::PointGroup::*;
    match pg {
        C1 => (None, None),
        C2 { axis: _ } => {
            todo!();
        }
        Cs { plane: _ } => {
            todo!()
        }
        C2v { axis: _, planes: _ } => {
            let mut checks = Checks::new();
            // first one you hit goes in checks.0, second goes in checks.1
            for i in irreps {
                match i.1 {
                    A1 => (),
                    B1 => {
                        if checks[(0, 0)] == 0 {
                            checks[(0, 0)] = i.0 + 1;
                            checks[(0, 1)] = i.0 + 1;
                        } else if i.0 + 1 > checks[(0, 1)] {
                            checks[(0, 1)] = i.0 + 1;
                        }
                    }
                    B2 => {
                        if checks[(1, 0)] == 0 {
                            checks[(1, 0)] = i.0 + 1;
                            checks[(1, 1)] = i.0 + 1;
                        } else if i.0 + 1 > checks[(1, 1)] {
                            checks[(1, 1)] = i.0 + 1;
                        }
                    }
                    A2 => {
                        if checks[(2, 0)] == 0 {
                            checks[(2, 0)] = i.0 + 1;
                            checks[(2, 1)] = i.0 + 1;
                        } else if i.0 + 1 > checks[(2, 1)] {
                            checks[(2, 1)] = i.0 + 1;
                        }
                    }
                    _ => panic!("non-C2v irrep found in C2v point group"),
                }
            }
            (Some(checks.clone()), Some(checks))
        }
    }
}

fn taylor_to_anpass(
    taylor: &Taylor,
    taylor_disps: &Vec<Vec<isize>>,
    energies: &[f64],
) -> Anpass {
    let mut disps = Vec::new();
    for disp in taylor_disps {
        for coord in disp {
            disps.push(*coord as f64 * TAYLOR_DISP_SIZE);
        }
    }
    let tdl = taylor_disps.len();
    let fl = taylor.forces.len();
    let mut fs = Vec::new();
    for row in &taylor.forces {
        for c in row {
            fs.push(*c as i32);
        }
    }
    Anpass {
        disps: anpass::Dmat::from_row_slice(tdl, disps.len() / tdl, &disps),
        energies: anpass::Dvec::from_row_slice(energies),
        exponents: na::DMatrix::from_column_slice(
            taylor.forces[0].len(),
            fl,
            &fs,
        ),
        bias: None,
    }
}

fn disp_to_intder(disps: &Vec<Vec<isize>>) -> Vec<Vec<f64>> {
    let mut ret = Vec::new();
    for disp in disps {
        let disp: Vec<_> =
            disp.iter().map(|i| *i as f64 * TAYLOR_DISP_SIZE).collect();
        ret.push(disp);
    }
    ret
}

pub fn optimize(geom: Geom) -> Geom {
    // TODO handle error
    let _ = std::fs::create_dir("opt");
    let opt = Job::new(
        Mopac::new("opt/opt".to_string(), None, Rc::new(geom), 0, &MOPAC_TMPL),
        0,
    );
    // TODO pass in submitter, via `run`, actually have to pass in the Program
    // too I think
    let submitter = LocalQueue {
        dir: "opt".to_string(),
    };
    submitter.optimize(opt)
}

type TaylorDisps = Vec<Vec<isize>>;
type AtomicNumbers = Vec<usize>;

pub fn generate_pts(
    geom: Geom,
    intder: &mut Intder,
) -> (Vec<Rc<Geom>>, Taylor, TaylorDisps, AtomicNumbers) {
    let mol = {
        let mut mol = Molecule::new(geom.xyz().unwrap().to_vec());
        mol.normalize();
        mol
    };
    let atomic_numbers = mol.atomic_numbers();
    let pg = mol.point_group();

    // load the initial intder
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

    // build and run the points using psqs
    // TODO handle error
    let _ = std::fs::create_dir("pts");
    let mut geoms = Vec::with_capacity(file07.len());
    for geom in file07 {
        // this is a bit unsightly, but I also don't want to duplicate the
        // `from_slices` code in psqs
        geoms.push(Rc::new(Geom::from(Molecule::from_slices(
            atomic_numbers.clone(),
            geom.as_slice(),
        ))));
    }
    (geoms, taylor, taylor_disps, atomic_numbers)
}

pub fn freqs(
    mut energies: Vec<f64>,
    intder: &mut Intder,
    taylor: &Taylor,
    taylor_disps: &TaylorDisps,
    atomic_numbers: &AtomicNumbers,
    spectro: &Spectro,
    gspectro_cmd: &String,
    spectro_cmd: &String,
    summary_cmd: &String,
) -> Summary {
    let min = energies.iter().cloned().reduce(f64::min).unwrap();
    for energy in energies.iter_mut() {
        *energy -= min;
    }

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
    // TODO handle error
    let _ = std::fs::create_dir("freqs");
    Intder::dump_fcs("freqs", &f2, &f3, &f4);

    // spectro
    let mut spectro = spectro.clone();
    spectro.geom = mol;
    spectro.write("freqs/spectro.in").unwrap();

    // run gspectro
    let spectro_arg = String::from("-cmd=") + spectro_cmd;
    std::process::Command::new(gspectro_cmd.clone())
        .arg(spectro_arg)
        .arg("freqs/spectro.in")
        .output()
        .unwrap();

    Summary::new(summary_cmd, "freqs/spectro2.out")
}

// TODO I might have a QFF trait in this package that I can implement on
// psqs::Programs to extend them; pass the common stuff into run and then make
// run a method on this trait

/// run a full qff, taking the configuration from `config`, the intder template
/// from `intder`, and the spectro template from `spectro`. Only the simple
/// internal and symmetry internal coordinates are read from the intder
/// template. The input options, weights, and curvils are copied from the
/// spectro template, but the geometry will be updated
pub fn run(config: &Config, intder: &Intder, spectro: &Spectro) -> Summary {
    // optimize the geometry
    let geom = if config.optimize {
        optimize(config.geometry.clone())
    } else {
        todo!();
    };

    let mut intder = intder.clone();
    let (geoms, taylor, taylor_disps, atomic_numbers) =
        generate_pts(geom, &mut intder);

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

    freqs(
        energies,
        &mut intder,
        &taylor,
        &taylor_disps,
        &atomic_numbers,
        spectro,
        &config.gspectro_cmd,
	&config.spectro_cmd,
	&config.summary_cmd,
    )
}
