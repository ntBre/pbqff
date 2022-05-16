pub mod config;

use psqs::program::Template;
use symm::{Irrep, PointGroup};
use taylor::Checks;

/// step size to take in each symmetry internal coordinate to determine its
/// irrep
pub const TAYLOR_DISP_SIZE: f64 = 0.005;

// TODO take this from a template file
pub const MOPAC_TMPL: Template = Template::from(
    "scfcrt=1.D-21 aux(precision=14) PM6 external=testfiles/params.dat",
);

/// generate the Taylor series mod and equivalence checks from `irreps` in `pg`
pub fn make_taylor_checks(
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

#[cfg(test)]
mod tests {
    use std::rc::Rc;

    use intder::Intder;
    use psqs::{
        program::{mopac::Mopac, Job},
        queue::{local::LocalQueue, Queue},
    };
    use symm::Molecule;
    use taylor::Taylor;

    use crate::{
        config::Config, make_taylor_checks, MOPAC_TMPL, TAYLOR_DISP_SIZE,
    };
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
        for (i, disp) in disps.iter().enumerate() {
            let m =
                Molecule::from_slices(atomic_numbers.clone(), disp.as_slice());
            irreps.push((i, m.irrep(&pg)));
        }
        // sort by irrep symmetry
        irreps.sort_by_key(|k| k.1);
        // generate checks
        let checks = make_taylor_checks(irreps, &pg);
        // run taylor.py to get fcs and disps
        let taylor = Taylor::new(5, nsic, checks.0, checks.1);
        let _disps = taylor.disps();
        let _anpass_bot = taylor.forces;
        // put these into intder and anpass
    }
}
