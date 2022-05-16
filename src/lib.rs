#![allow(dead_code)]
pub mod config;

#[cfg(test)]
mod tests;

use anpass::Anpass;
use nalgebra as na;
use psqs::program::Template;
use rust_anpass as anpass;
use symm::{Irrep, PointGroup};
use taylor::{Checks, Taylor};

/// step size to take in each symmetry internal coordinate to determine its
/// irrep
pub const TAYLOR_DISP_SIZE: f64 = 0.005;

// TODO take this from a template file
pub const MOPAC_TMPL: Template = Template::from(
    "scfcrt=1.D-21 aux(precision=14) PM6 external=testfiles/params.dat",
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
        exponents: na::DMatrix::from_row_slice(fl, taylor.forces[0].len(), &fs),
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
