use self::bighash::BigHash;

use super::{cart::DEBUG, CartGeom, Derivative};
use bighash::Index;
use intder::ANGBOHR;
use nalgebra as na;
use psqs::geom::Geom;
use std::cmp::min;
use symm::{Atom, Molecule};

pub mod bighash;

/// geom is None if no displacement is required, i.e. this is the reference
/// geometry
#[derive(Debug, PartialEq)]
pub struct Proto {
    pub geom: Option<Geom>,
    pub coeff: f64,
}

pub(crate) fn zip_atoms(names: &[&str], coords: na::DVector<f64>) -> Vec<Atom> {
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

macro_rules! geom {
    ($s:ident, $names:ident, $coords:ident, $step_size:ident, $($steps: expr),* ) => {
	Some($s.new_geom($names.clone(), $coords.clone(), $step_size, vec![$($steps),*]))
    };
}

pub(crate) use geom;

macro_rules! proto {
    ($s:ident, $names: ident, $coords: ident, $step_size: ident, $scale: expr, $( $steps: expr ),* ) => {
	Proto {
	    geom: $crate::coord_type::findiff::geom!($s, $names, $coords, $step_size, $($steps),*),
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

pub(crate) use proto;

pub(crate) fn atom_parts(atoms: &Vec<Atom>) -> (Vec<&str>, Vec<f64>) {
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

pub type Idx = (usize, usize, usize, usize);

/// a trait for CoordTypes that rely on finite differences instead of fittings
pub trait FiniteDifference {
    fn new_geom(
        &self,
        names: &[&str],
        coords: na::DVector<f64>,
        step_size: f64,
        steps: Vec<isize>,
    ) -> Geom;

    /// generate the points for `self` using the reference `geom` and reference
    /// energy. `fcs` and `map` are out parameters. `fcs` will be filled
    /// directly with derivatives that only rely on the reference energy, and
    /// `map` will describe how where the rest of the energies in the return
    /// value should end up
    #[allow(clippy::too_many_arguments)]
    fn build_points(
        &self,
        geom: Geom,
        step_size: f64,
        ref_energy: f64,
        deriv: Derivative,
        fcs: &mut [f64],
        map: &mut BigHash,
        ncoords: usize,
    ) -> Vec<CartGeom> {
        let atoms = geom.xyz().unwrap();
        let (names, coords) = atom_parts(atoms);
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
                            (0, 0) => {
                                self.make2d(&names, &coords, step_size, i, j)
                            }
                            (_, 0) => {
                                self.make3d(&names, &coords, step_size, i, j, k)
                            }
                            (_, _) => self.make4d(
                                &names,
                                &coords,
                                step_size,
                                (i, j, k, l),
                            ),
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

    /// use `target_map` to map `energies` to their proper locations in `fcs` to
    /// yield the force constants as one vector
    fn map_energies(
        &self,
        target_map: &BigHash,
        energies: &[f64],
        fcs: &mut [f64],
    ) {
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
    }

    /// use `target_map` to symmetrize the force constants and return them in
    /// the form wanted by spectro
    fn make_fcs<'a>(
        &self,
        target_map: &mut BigHash,
        energies: &[f64],
        fcs: &'a mut [f64],
        n: usize,
        deriv: Derivative,
        dir: &str,
    ) -> (nalgebra::DMatrix<f64>, &'a [f64], &'a [f64]) {
        self.map_energies(target_map, energies, fcs);
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
        let (nfc2, nfc3) = match deriv {
            Derivative::Harmonic(n) => (n, 0),
            Derivative::Cubic(n, m) => (n, m),
            Derivative::Quartic(n, m, _) => (n, m),
        };
        intder::Intder::dump_fcs(
            dir,
            &fc2,
            &fcs[nfc2..nfc2 + nfc3],
            &fcs[nfc2 + nfc3..],
        );

        (fc2, &fcs[nfc2..nfc2 + nfc3], &fcs[nfc2 + nfc3..])
    }

    fn scale(&self, nderiv: usize, step_size: f64) -> f64 {
        match nderiv {
            2 => ANGBOHR * ANGBOHR / (4.0 * step_size * step_size),
            3 => {
                ANGBOHR * ANGBOHR * ANGBOHR
                    / (8.0 * step_size * step_size * step_size)
            }
            4 => {
                ANGBOHR * ANGBOHR * ANGBOHR * ANGBOHR
                    / (16.0 * step_size * step_size * step_size * step_size)
            }
            _ => panic!("unrecognized derivative level"),
        }
    }

    fn make2d(
        &self,
        names: &[&str],
        coords: &na::DVector<f64>,
        step_size: f64,
        i: usize,
        j: usize,
    ) -> Vec<Proto> {
        let scale = self.scale(2, step_size);
        let i = i as isize;
        let j = j as isize;
        if i == j {
            vec![
                proto!(self, names, coords, step_size, scale, i, i),
                proto!(-2. * scale),
                proto!(self, names, coords, step_size, scale, -i, -i),
            ]
        } else {
            vec![
                proto!(self, names, coords, step_size, scale, i, j),
                proto!(self, names, coords, step_size, -scale, i, -j),
                proto!(self, names, coords, step_size, -scale, -i, j),
                proto!(self, names, coords, step_size, scale, -i, -j),
            ]
        }
    }

    fn make3d(
        &self,
        names: &[&str],
        coords: &na::DVector<f64>,
        step_size: f64,
        i: usize,
        j: usize,
        k: usize,
    ) -> Vec<Proto> {
        let scale = self.scale(3, step_size);
        let i = i as isize;
        let j = j as isize;
        let k = k as isize;

        let make3d_2_1 = |i, _, k| {
            vec![
                proto!(self, names, coords, step_size, scale, i, i, k),
                proto!(self, names, coords, step_size, -2. * scale, k),
                proto!(self, names, coords, step_size, scale, -i, -i, k),
                proto!(self, names, coords, step_size, -scale, i, i, -k),
                proto!(self, names, coords, step_size, 2. * scale, -k),
                proto!(self, names, coords, step_size, -scale, -i, -i, -k),
            ]
        };
        if i == j && i == k {
            vec![
                proto!(self, names, coords, step_size, scale, i, i, i),
                proto!(self, names, coords, step_size, -3. * scale, i),
                proto!(self, names, coords, step_size, 3. * scale, -i),
                proto!(self, names, coords, step_size, -scale, -i, -i, -i),
            ]
        } else if i == j {
            make3d_2_1(i, j, k)
        } else if i == k {
            make3d_2_1(i, k, j) // unreachable
        } else if j == k {
            make3d_2_1(j, k, i)
        } else {
            vec![
                proto!(self, names, coords, step_size, 1. * scale, i, j, k),
                proto!(self, names, coords, step_size, -1. * scale, i, -j, k),
                proto!(self, names, coords, step_size, -1. * scale, -i, j, k),
                proto!(self, names, coords, step_size, 1. * scale, -i, -j, k),
                proto!(self, names, coords, step_size, -1. * scale, i, j, -k),
                proto!(self, names, coords, step_size, 1. * scale, i, -j, -k),
                proto!(self, names, coords, step_size, 1. * scale, -i, j, -k),
                proto!(self, names, coords, step_size, -1. * scale, -i, -j, -k),
            ]
        }
    }

    fn make4d(
        &self,
        names: &[&str],
        coords: &na::DVector<f64>,
        step: f64,
        idx: Idx,
    ) -> Vec<Proto> {
        let scale = self.scale(4, step);
        let (i, j, k, l) = idx;
        let i = i as isize;
        let j = j as isize;
        let k = k as isize;
        let l = l as isize;

        let make4d_3_1 = |i, _, _, l| {
            vec![
                proto!(self, names, coords, step, 1. * scale, i, i, i, l),
                proto!(self, names, coords, step, -3. * scale, i, l),
                proto!(self, names, coords, step, 3. * scale, -i, l),
                proto!(self, names, coords, step, -1. * scale, -i, -i, -i, l),
                proto!(self, names, coords, step, -1. * scale, i, i, i, -l),
                proto!(self, names, coords, step, 3. * scale, i, -l),
                proto!(self, names, coords, step, -3. * scale, -i, -l),
                proto!(self, names, coords, step, 1. * scale, -i, -i, -i, -l),
            ]
        };

        let make4d_2_2 = |i, _, k, _| {
            vec![
                proto!(self, names, coords, step, 1. * scale, i, i, k, k),
                proto!(self, names, coords, step, 1. * scale, -i, -i, -k, -k),
                proto!(self, names, coords, step, 1. * scale, -i, -i, k, k),
                proto!(self, names, coords, step, 1. * scale, i, i, -k, -k),
                proto!(self, names, coords, step, -2. * scale, i, i),
                proto!(self, names, coords, step, -2. * scale, k, k),
                proto!(self, names, coords, step, -2. * scale, -i, -i),
                proto!(self, names, coords, step, -2. * scale, -k, -k),
                proto!(4. * scale),
            ]
        };

        let make4d_2_1_1 = |i, _, k, l| {
            vec![
                proto!(self, names, coords, step, 1. * scale, i, i, k, l),
                proto!(self, names, coords, step, -2. * scale, k, l),
                proto!(self, names, coords, step, 1. * scale, -i, -i, k, l),
                proto!(self, names, coords, step, -1. * scale, i, i, -k, l),
                proto!(self, names, coords, step, 2. * scale, -k, l),
                proto!(self, names, coords, step, -1. * scale, -i, -i, -k, l),
                proto!(self, names, coords, step, -1. * scale, i, i, k, -l),
                proto!(self, names, coords, step, 2. * scale, k, -l),
                proto!(self, names, coords, step, -1. * scale, -i, -i, k, -l),
                proto!(self, names, coords, step, 1. * scale, i, i, -k, -l),
                proto!(self, names, coords, step, -2. * scale, -k, -l),
                proto!(self, names, coords, step, 1. * scale, -i, -i, -k, -l),
            ]
        };

        if i == j && i == k && i == l {
            vec![
                proto!(self, names, coords, step, 1. * scale, i, i, i, i),
                proto!(self, names, coords, step, -4. * scale, i, i),
                proto!(6. * scale),
                proto!(self, names, coords, step, -4. * scale, -i, -i),
                proto!(self, names, coords, step, 1. * scale, -i, -i, -i, -i),
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
                proto!(self, names, coords, step, 1. * scale, i, j, k, l),
                proto!(self, names, coords, step, -1. * scale, i, -j, k, l),
                proto!(self, names, coords, step, -1. * scale, -i, j, k, l),
                proto!(self, names, coords, step, 1. * scale, -i, -j, k, l),
                proto!(self, names, coords, step, -1. * scale, i, j, -k, l),
                proto!(self, names, coords, step, 1. * scale, i, -j, -k, l),
                proto!(self, names, coords, step, 1. * scale, -i, j, -k, l),
                proto!(self, names, coords, step, -1. * scale, -i, -j, -k, l),
                proto!(self, names, coords, step, -1. * scale, i, j, k, -l),
                proto!(self, names, coords, step, 1. * scale, i, -j, k, -l),
                proto!(self, names, coords, step, 1. * scale, -i, j, k, -l),
                proto!(self, names, coords, step, -1. * scale, -i, -j, k, -l),
                proto!(self, names, coords, step, 1. * scale, i, j, -k, -l),
                proto!(self, names, coords, step, -1. * scale, i, -j, -k, -l),
                proto!(self, names, coords, step, -1. * scale, -i, j, -k, -l),
                proto!(self, names, coords, step, 1. * scale, -i, -j, -k, -l),
            ]
        }
    }
}
