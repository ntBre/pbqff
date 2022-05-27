#![allow(unused)]
use std::rc::Rc;

use intder::ANGBOHR;
use psqs::{
    geom::Geom,
    program::{mopac::{Mopac, KCALHT}, Job, Template},
    queue::{local::LocalQueue, Queue},
};
use spectro::Spectro;
use summarize::Summary;
use symm::Atom;

use nalgebra as na;

use crate::{config::Config, optimize};

use super::CoordType;

const MOPAC_TMPL: Template = Template::from(
    "scfcrt=1.D-21 aux(precision=14) PM6 external=testfiles/params.dat",
);

pub struct Cart;

fn atom_parts(atoms: &Vec<Atom>) -> (Vec<&str>, Vec<f64>) {
    let mut names = Vec::new();
    let mut coords = Vec::new();
    for atom in atoms {
        names.push(atom.label());
        coords.extend(atom.coord());
    }
    (names, coords)
}

struct CartGeom {
    geom: Geom,
    coeff: f64,
    index: (usize, usize, usize, usize),
}

/// geom is None if no displacement is required, i.e. this is the reference
/// geometry
struct Proto {
    geom: Option<Geom>,
    coeff: f64,
}

fn zip_atoms(names: Vec<&str>, coords: na::DVector<f64>) -> Vec<Atom> {
    // this makes sure they match and that coords is divisible by 3
    assert!(3 * names.len() == coords.len());
    names
        .iter()
        .zip(coords.as_slice().chunks(3))
        .map(|(name, coord)| {
            Atom::new_from_label(name, coord[0], coord[1], coord[2])
        })
        .collect()
}

fn new_geom(
    names: Vec<&str>,
    coords: na::DVector<f64>,
    step_size: f64,
    steps: Vec<isize>,
) -> Geom {
    let mut v = vec![0.0; coords.len()];
    for step in steps {
        if step < 1 {
            v[(-step - 1) as usize] -= step_size;
        } else {
            v[(step - 1) as usize] += step_size;
        }
    }
    let coords = coords + na::DVector::from(v);
    Geom::Xyz(zip_atoms(names, coords))
}

macro_rules! geom {
    ( $names: ident, $coords: ident, $step_size: ident, $( $steps: expr ),* ) => {
        Some(new_geom($names.clone(), $coords.clone(), $step_size, vec![$($steps),*]))
    };
}

macro_rules! proto {
    ( $names: ident, $coords: ident, $step_size: ident, $scale: expr, $( $steps: expr ),* ) => {
	Proto {
	    geom: geom!($names, $coords, $step_size, $($steps),*),
	    coeff: $scale,
	}
    };
    ( $coeff: expr ) => {
	Proto {
	    geom: None,
	    coeff: $coeff,
	}
    }
}

fn make2d(
    names: Vec<&str>,
    coords: na::DVector<f64>,
    step_size: f64,
    i: usize,
    j: usize,
) -> Vec<Proto> {
    let scale = ANGBOHR * ANGBOHR / (4.0 * step_size * step_size);
    let i = i as isize;
    let j = j as isize;
    if i == j {
        vec![
            proto!(names, coords, step_size, scale, i, i),
            proto!(-2. * scale),
            proto!(names, coords, step_size, scale, -i, -i),
        ]
    } else {
        vec![
            proto!(names, coords, step_size, scale, i, j),
            proto!(names, coords, step_size, -scale, i, -j),
            proto!(names, coords, step_size, -scale, -i, j),
            proto!(names, coords, step_size, scale, -i, -j),
        ]
    }
}

fn make3d(
    names: Vec<&str>,
    coords: na::DVector<f64>,
    step_size: f64,
    i: usize,
    j: usize,
    k: usize,
) -> Vec<Proto> {
    let scale =
        ANGBOHR * ANGBOHR * ANGBOHR / (8.0 * step_size * step_size * step_size);
    let i = i as isize;
    let j = j as isize;
    let k = k as isize;

    let make3d_2_1 = |i, _, k| {
        vec![
            proto!(names, coords, step_size, scale, i, i, k),
            proto!(names, coords, step_size, -2. * scale, k),
            proto!(names, coords, step_size, scale, -i, -i, k),
            proto!(names, coords, step_size, -scale, i, i, -k),
            proto!(names, coords, step_size, 2. * scale, -k),
            proto!(names, coords, step_size, -scale, -i, -i, -k),
        ]
    };
    if i == j && i == k {
        vec![
            proto!(names, coords, step_size, scale, i, i, i),
            proto!(names, coords, step_size, -3. * scale, i),
            proto!(names, coords, step_size, 3. * scale, -i),
            proto!(names, coords, step_size, -scale, -i, -i, -i),
        ]
    } else if i == j {
        make3d_2_1(i, j, k)
    } else if i == k {
        make3d_2_1(i, k, j)
    } else if j == k {
        make3d_2_1(j, k, i)
    } else {
        vec![
            proto!(names, coords, step_size, 1. * scale, i, j, k),
            proto!(names, coords, step_size, -1. * scale, i, -j, k),
            proto!(names, coords, step_size, -1. * scale, -i, j, k),
            proto!(names, coords, step_size, 1. * scale, -i, -j, k),
            proto!(names, coords, step_size, -1. * scale, i, j, -k),
            proto!(names, coords, step_size, 1. * scale, i, -j, -k),
            proto!(names, coords, step_size, 1. * scale, -i, j, -k),
            proto!(names, coords, step_size, -1. * scale, -i, -j, -k),
        ]
    }
}

fn make4d(
    names: Vec<&str>,
    coords: na::DVector<f64>,
    step_size: f64,
    i: usize,
    j: usize,
    k: usize,
    l: usize,
) -> Vec<Proto> {
    let scale = ANGBOHR * ANGBOHR * ANGBOHR * ANGBOHR
        / (16.0 * step_size * step_size * step_size * step_size);
    let i = i as isize;
    let j = j as isize;
    let k = k as isize;
    let l = l as isize;

    let make4d_3_1 = |i, _, _, l| {
        vec![
            proto!(names, coords, step_size, 1. * scale, i, i, i, l),
            proto!(names, coords, step_size, -3. * scale, i, l),
            proto!(names, coords, step_size, 3. * scale, -i, l),
            proto!(names, coords, step_size, -1. * scale, -i, -i, -i, l),
            proto!(names, coords, step_size, -1. * scale, i, i, i, -l),
            proto!(names, coords, step_size, 3. * scale, i, -l),
            proto!(names, coords, step_size, -3. * scale, -i, -l),
            proto!(names, coords, step_size, 1. * scale, -i, -i, -i, -l),
        ]
    };

    let make4d_2_2 = |i, _, k, _| {
        vec![
            proto!(names, coords, step_size, 1. * scale, i, i, k, k),
            proto!(names, coords, step_size, 1. * scale, -i, -i, -k, -k),
            proto!(names, coords, step_size, 1. * scale, -i, -i, k, k),
            proto!(names, coords, step_size, 1. * scale, i, i, -k, -k),
            proto!(names, coords, step_size, -2. * scale, i, i),
            proto!(names, coords, step_size, -2. * scale, k, k),
            proto!(names, coords, step_size, -2. * scale, -i, -i),
            proto!(names, coords, step_size, -2. * scale, -k, -k),
            proto!(4. * scale),
        ]
    };

    let make4d_2_1_1 = |i, _, k, l| {
        vec![
            proto!(names, coords, step_size, 1. * scale, i, i, k, l),
            proto!(names, coords, step_size, -2. * scale, k, l),
            proto!(names, coords, step_size, 1. * scale, -i, -i, k, l),
            proto!(names, coords, step_size, -1. * scale, i, i, -k, l),
            proto!(names, coords, step_size, 2. * scale, -k, l),
            proto!(names, coords, step_size, -1. * scale, -i, -i, -k, l),
            proto!(names, coords, step_size, -1. * scale, i, i, k, -l),
            proto!(names, coords, step_size, 2. * scale, k, -l),
            proto!(names, coords, step_size, -1. * scale, -i, -i, k, -l),
            proto!(names, coords, step_size, 1. * scale, i, i, -k, -l),
            proto!(names, coords, step_size, -2. * scale, -k, -l),
            proto!(names, coords, step_size, 1. * scale, -i, -i, -k, -l),
        ]
    };

    if i == j && i == k && i == l {
        vec![
            proto!(names, coords, step_size, 1. * scale, i, i, i, i),
            proto!(names, coords, step_size, -4. * scale, i, i),
            proto!(6. * scale),
            proto!(names, coords, step_size, -4. * scale, -i, -i),
            proto!(names, coords, step_size, 1. * scale, -i, -i, -i, -i),
        ]
    // 3 and 1
    } else if i == j && i == k {
        make4d_3_1(i, j, k, l)
    } else if i == j && i == l {
        make4d_3_1(i, j, l, k)
    } else if i == k && i == l {
        make4d_3_1(i, k, l, j)
    } else if j == k && j == l {
        make4d_3_1(j, k, l, i)
    // 2 and 2
    } else if i == j && k == l {
        make4d_2_2(i, j, k, l)
    } else if i == k && j == l {
        make4d_2_2(i, k, j, l)
    } else if i == l && j == k {
        make4d_2_2(i, l, j, k)
    // 2 and 1 and 1
    } else if i == j {
        make4d_2_1_1(i, j, k, l)
    } else if i == k {
        make4d_2_1_1(i, j, k, l)
    } else if i == l {
        make4d_2_1_1(i, j, k, l)
    } else if j == k {
        make4d_2_1_1(i, j, k, l)
    } else if j == l {
        make4d_2_1_1(i, j, k, l)
    } else if k == l {
        make4d_2_1_1(i, j, k, l)
    } else {
        vec![
            proto!(names, coords, step_size, 1. * scale, i, j, k, l),
            proto!(names, coords, step_size, -1. * scale, i, -j, k, l),
            proto!(names, coords, step_size, -1. * scale, -i, j, k, l),
            proto!(names, coords, step_size, 1. * scale, -i, -j, k, l),
            proto!(names, coords, step_size, -1. * scale, i, j, -k, l),
            proto!(names, coords, step_size, 1. * scale, i, -j, -k, l),
            proto!(names, coords, step_size, 1. * scale, -i, j, -k, l),
            proto!(names, coords, step_size, -1. * scale, -i, -j, -k, l),
            proto!(names, coords, step_size, -1. * scale, i, j, k, -l),
            proto!(names, coords, step_size, 1. * scale, i, -j, k, -l),
            proto!(names, coords, step_size, 1. * scale, -i, j, k, -l),
            proto!(names, coords, step_size, -1. * scale, -i, -j, k, -l),
            proto!(names, coords, step_size, 1. * scale, i, j, -k, -l),
            proto!(names, coords, step_size, -1. * scale, i, -j, -k, -l),
            proto!(names, coords, step_size, -1. * scale, -i, j, -k, -l),
            proto!(names, coords, step_size, 1. * scale, -i, -j, -k, -l),
        ]
    }
}

impl Cart {
    fn build_points(
        &self,
        geom: Geom,
        step_size: f64,
        ref_energy: f64,
        dest: &mut [f64],
    ) -> Vec<CartGeom> {
        let atoms = geom.xyz().unwrap();
        let (names, coords) = atom_parts(atoms);
        let ncoords = coords.len();
        let coords = na::DVector::from(coords);

        let mut ret = Vec::new();
        // start at 1 so that k = l = 0 indicates second derivative
        for i in 1..=ncoords {
            for j in 1..=i {
                for k in 0..=j {
                    for l in 0..=k {
                        let protos = match (k, l) {
                            (0, 0) => make2d(
                                names.clone(),
                                coords.clone(),
                                step_size,
                                i,
                                j,
                            ),
                            // (_, 0) => make3d(
                            //     names.clone(),
                            //     coords.clone(),
                            //     step_size,
                            //     i,
                            //     j,
                            //     k,
                            // ),
                            // (_, _) => make4d(
                            //     names.clone(),
                            //     coords.clone(),
                            //     step_size,
                            //     i,
                            //     j,
                            //     k,
                            //     l,
                            // ),
                            (_, _) => vec![],
                        };

                        // TODO handle symmetry
                        let idx = (i, j, k, l);
                        for p in protos {
                            if p.geom.is_none() {
                                let (i, j, k, l) = idx;
                                dest[index(ncoords, i, j, k, l)] += p.coeff * ref_energy;
                            } else {
                                ret.push(CartGeom {
                                    geom: p.geom.unwrap(),
                                    coeff: p.coeff,
                                    index: idx,
                                });
                            }
                        }
                    }
                }
            }
        }
        ret
    }

    fn freqs(&self, _spectro: &Spectro) -> Summary {
        Summary::new("spectro2.out")
    }
}

/// compute the index in the force constant array
fn index(n: usize, a: usize, b: usize, c: usize, d: usize) -> usize {
    match (c, d) {
        (0, 0) => intder::fc2_index(n, a, b),
        (_, 0) => intder::fc3_index(a, b, c),
        (_, _) => intder::fc4_index(a, b, c, d),
    }
}

impl CoordType for Cart {
    fn run(&self, config: &Config, spectro: &Spectro) -> Summary {
        let geom = if config.optimize {
            // TODO take the reference energy from the same calculation if
            // optimizing anyway
            optimize(config.geometry.clone())
        } else {
            todo!();
        };

	println!("{}", &geom);

        // 3 * #atoms
        let n = 3 * geom.xyz().unwrap().len();
        let nfc2 = n * n;
        let nfc3 = n * (n + 1) * (n + 2) / 6;
        let nfc4 = n * (n + 1) * (n + 2) * (n + 3) / 24;
        let mut fcs = vec![0.0; nfc2 + nfc3 + nfc4];

        let ref_energy = 0.12660293116764660226E+03 / KCALHT; // TODO
        let geoms =
            self.build_points(geom, config.step_size, ref_energy, &mut fcs);

        let mut jobs = {
            let dir = "pts";
            let mut job_num = 0;
            let mut jobs = Vec::new();
            for mol in geoms {
                let filename = format!("{dir}/job.{:08}", job_num);
                job_num += 1;
                let (i, j, k, l) = mol.index;
                let index = index(n, i, j, k, l);
                let mut job = Job::new(
                    Mopac::new(
                        filename,
                        None,
                        Rc::new(mol.geom),
                        config.charge,
                        &MOPAC_TMPL,
                    ),
                    index,
                );
                job.coeff = mol.coeff;
                jobs.push(job);
            }
            jobs
        };

        LocalQueue {
            dir: "pts".to_string(),
        }
        .drain(&mut jobs, &mut fcs);

        // mirror symmetric quadratic fcs
        let mut fc2 = intder::DMat::zeros(n, n);
        for i in 0..n {
            fc2[(i, i)] = fcs[intder::fc2_index(n, i + 1, i + 1)];
            for j in 0..i {
                let f = fcs[intder::fc2_index(n, i + 1, j + 1)];
                fc2[(i, j)] = f;
                fc2[(j, i)] = f;
            }
        }
        println!("{:.8}", &fc2);

        intder::Intder::dump_fcs(".", &fc2, &fcs[nfc2..nfc3], &fcs[nfc3..]);

        todo!();

        self.freqs(spectro)
    }
}
