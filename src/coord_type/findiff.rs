use super::{CartGeom, CoordType, Derivative, DEBUG};
use bighash::Index;
use intder::ANGBOHR;
use nalgebra as na;
use psqs::{geom::Geom, program::Program, queue::Queue};
use std::{cmp::min, io::Write, marker::Sync};
use symm::{Atom, Molecule};

pub mod bighash;

/// geom is None if no displacement is required, i.e. this is the reference
/// geometry
#[derive(Debug, PartialEq)]
struct Proto {
    geom: Option<Geom>,
    coeff: f64,
}

fn zip_atoms(names: &[&str], coords: na::DVector<f64>) -> Vec<Atom> {
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
    names: &[&str],
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
    names: &[&str],
    coords: &na::DVector<f64>,
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
    names: &[&str],
    coords: &na::DVector<f64>,
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
        make3d_2_1(i, k, j) // unreachable
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

// there is an issue with the all same case => index 0, also an issue with the
// 2-2 case => index 2
fn make4d(
    names: &[&str],
    coords: &na::DVector<f64>,
    step_size: f64,
    i: usize,
    j: usize,
    k: usize,
    l: usize,
) -> Vec<Proto> {
    let scale = ANGBOHR.powi(4) / (16.0 * step_size.powi(4));
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
        make4d_3_1(i, j, l, k) // unreachable
    } else if i == k && i == l {
        make4d_3_1(i, k, l, j) // unreachable
    } else if j == k && j == l {
        make4d_3_1(j, k, l, i)
    // 2 and 2
    } else if i == j && k == l {
        make4d_2_2(i, j, k, l)
    } else if i == k && j == l {
        make4d_2_2(i, k, j, l) // unreachable
    } else if i == l && j == k {
        make4d_2_2(i, l, j, k) // unreachable
                               // 2 and 1 and 1, first two are the equal ones
    } else if i == j {
        make4d_2_1_1(i, j, k, l)
    } else if i == k {
        make4d_2_1_1(i, k, j, l) // unreachable
    } else if i == l {
        make4d_2_1_1(i, l, j, k) // unreachable
    } else if j == k {
        make4d_2_1_1(j, k, i, l)
    } else if j == l {
        make4d_2_1_1(j, l, i, k) // unreachable
    } else if k == l {
        make4d_2_1_1(k, l, i, j)
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

fn atom_parts(atoms: &Vec<Atom>) -> (Vec<&str>, Vec<f64>) {
    let mut names = Vec::new();
    let mut coords = Vec::new();
    for atom in atoms {
        names.push(atom.label());
        coords.extend(atom.coord());
    }
    (names, coords)
}

/// compute the index in the force constant array
fn index(n: usize, a: usize, b: usize, c: usize, d: usize) -> usize {
    match (c, d) {
        (0, 0) => intder::fc2_index(n, a, b),
        (_, 0) => intder::fc3_index(a, b, c),
        (_, _) => intder::fc4_index(a, b, c, d),
    }
}

/// a trait for CoordTypes that rely on finite differences instead of fittings
pub trait FiniteDifference<
    W: Write,
    Q: Queue<P>,
    P: Program + Clone + Send + Sync,
>: CoordType<W, Q, P>
{
    fn build_points(
        &self,
        geom: Geom,
        step_size: f64,
        ref_energy: f64,
        deriv: Derivative,
        fcs: &mut [f64],
        map: &mut bighash::BigHash,
    ) -> Vec<CartGeom> {
        let atoms = geom.xyz().unwrap();
        let (names, coords) = atom_parts(atoms);
        let ncoords = coords.len();
        let coords = nalgebra::DVector::from(coords);

        let (nfc2, nfc3, k_max, l_max) = match deriv {
            Derivative::Harmonic(x) => (x, 0, 0, 0),
            Derivative::Cubic(x, y) => (x, y, ncoords, 0),
            Derivative::Quartic(x, y, _) => (x, y, ncoords, ncoords),
        };

        let mut geoms = Vec::new();
        // counter is the index into the the energies array that the jobs will
        // be run into
        let mut counter = 0;
        // start at 1 so that k = l = 0 indicates second derivative
        for i in 1..=ncoords {
            for j in 1..=i {
                for k in 0..=min(j, k_max) {
                    for l in 0..=min(k, l_max) {
                        let protos = match (k, l) {
                            (0, 0) => make2d(&names, &coords, step_size, i, j),
                            (_, 0) => {
                                make3d(&names, &coords, step_size, i, j, k)
                            }
                            (_, _) => {
                                make4d(&names, &coords, step_size, i, j, k, l)
                            }
                        };
                        let idx = (i, j, k, l);
                        for p in protos {
                            let (i, j, k, l) = idx;
                            // index is the index into the force constant array
                            let mut index = index(ncoords, i, j, k, l);
                            match (k, l) {
                                (0, 0) => (),
                                (_, 0) => index += nfc2,
                                (_, _) => index += nfc2 + nfc3,
                            };
                            if let Some(geom) = p.geom.clone() {
                                let mol =
                                    Molecule::new(geom.xyz().unwrap().to_vec());
                                // if the structure has already been seen, push
                                // this target index to its list of Targets
                                if let Some(result) = map.get_mut(&mol) {
                                    result.indices.push(bighash::Index {
                                        index,
                                        coeff: p.coeff,
                                    });
                                } else {
                                    // otherwise, mark the structure as seen by
                                    // inserting it into the map and return a
                                    // job to be run
                                    let i = map.insert(
                                        mol,
                                        bighash::Target {
                                            source_index: counter,
                                            indices: vec![Index {
                                                index,
                                                coeff: p.coeff,
                                            }],
                                        },
                                    );
                                    assert_eq!(i, None);
                                    geoms.push(CartGeom {
                                        geom: p.geom.unwrap(),
                                        coeff: p.coeff,
                                        index: counter,
                                    });
                                    counter += 1;
                                }
                            } else {
                                // reference energy, handle it directly in fcs
                                fcs[index] += p.coeff * ref_energy;
                            }
                        }
                    }
                }
            }
        }
        geoms
    }

    /// use `target_map` to symmetrize the force constants and return them in
    /// the form wanted by spectro
    fn make_fcs<'a>(
        &self,
        target_map: &mut bighash::BigHash,
        energies: &[f64],
        fcs: &'a mut [f64],
        n: usize,
        deriv: Derivative,
        dir: &str,
    ) -> (nalgebra::DMatrix<f64>, &'a [f64], &'a [f64]) {
        // copy energies into all of the places they're needed
        for target in target_map.values() {
            if DEBUG == "fcs" {
                eprintln!("source index: {}", target.source_index);
            }
            let energy = energies[target.source_index];
            for idx in &target.indices {
                if DEBUG == "fcs" {
                    eprintln!(
                        "\tfcs[{}] += {:12.8} * {:12.8}",
                        idx.index, idx.coeff, energy
                    );
                }
                fcs[idx.index] += idx.coeff * energy;
            }
        }
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

        let _ = std::fs::create_dir(dir);
        let Derivative::Quartic(nfc2, nfc3, _) = deriv else { todo!() };
        intder::Intder::dump_fcs(
            dir,
            &fc2,
            &fcs[nfc2..nfc2 + nfc3],
            &fcs[nfc2 + nfc3..],
        );

        (fc2, &fcs[nfc2..nfc2 + nfc3], &fcs[nfc2 + nfc3..])
    }
}
