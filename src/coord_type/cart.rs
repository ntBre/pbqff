use std::{
    collections::{hash_map::Values, HashMap},
    hash::Hash,
    io,
    path::Path,
    rc::Rc,
};

use intder::ANGBOHR;
use psqs::{
    geom::Geom,
    program::{mopac::Mopac, Job, Template},
    queue::Queue,
};
use spectro::{Output, Spectro};
use symm::{Atom, Molecule, PointGroup};

use nalgebra as na;

use crate::{config::Config, optimize, ref_energy};

use super::CoordType;

/// debugging options. currently supported options: disp, fcs, none
pub(crate) static DEBUG: &str = "none";

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

#[derive(Debug)]
pub struct CartGeom {
    pub geom: Geom,
    pub coeff: f64,
    pub index: usize,
}

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

/// the fields are options because it's possible for `detect_buddies` to fail
/// and find that not all atoms are matched. if this happens, we want to skip
/// over those entirely instead of filling the Hash with junk geometries
#[derive(Clone, Debug)]
enum Buddy {
    C1,
    C2 {
        axis: Option<Vec<usize>>,
    },
    Cs {
        plane: Option<Vec<usize>>,
    },
    C2v {
        axis: Option<Vec<usize>>,
        plane0: Option<Vec<usize>>,
        plane1: Option<Vec<usize>>,
    },
    D2h {
        axes: Vec<Option<Vec<usize>>>,
        planes: Vec<Option<Vec<usize>>>,
    },
}

impl Buddy {
    /// apply a Buddy to a Molecule, generating all of the symmetry-equivalent
    /// molecules
    fn apply(&self, mol: &Molecule) -> Vec<Molecule> {
        let mut ret = Vec::new();
        match self {
            Buddy::C1 => (),
            Buddy::C2 { axis } => {
                if let Some(axis) = axis {
                    ret.push(Molecule::new(
                        axis.iter().map(|i| mol.atoms[*i]).collect(),
                    ));
                }
            }
            Buddy::Cs { plane } => {
                if let Some(plane) = plane {
                    ret.push(Molecule::new(
                        plane.iter().map(|i| mol.atoms[*i]).collect(),
                    ));
                }
            }
            Buddy::C2v {
                axis,
                plane0,
                plane1,
            } => {
                if let Some(axis) = axis {
                    ret.push(Molecule::new(
                        axis.iter().map(|i| mol.atoms[*i]).collect(),
                    ));
                }
                if let Some(plane0) = plane0 {
                    ret.push(Molecule::new(
                        plane0.iter().map(|i| mol.atoms[*i]).collect(),
                    ));
                }
                if let Some(plane1) = plane1 {
                    ret.push(Molecule::new(
                        plane1.iter().map(|i| mol.atoms[*i]).collect(),
                    ));
                }
            }
            Buddy::D2h { axes, planes } => {
                for axis in axes.iter().flatten() {
                    ret.push(Molecule::new(
                        axis.iter().map(|i| mol.atoms[*i]).collect(),
                    ));
                }
                for plane in planes.iter().flatten() {
                    ret.push(Molecule::new(
                        plane.iter().map(|i| mol.atoms[*i]).collect(),
                    ));
                }
            }
        }
        ret
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Key {
    atom: usize,
    x: isize,
    y: isize,
    z: isize,
}

#[derive(Clone, Debug)]
pub struct BigHash {
    map: HashMap<Vec<Key>, Target>,
    pg: PointGroup,
    // buddies are pairs of atoms that are interchanged across symmetry
    // operations
    buddy: Buddy,
}

impl BigHash {
    /// NOTE: assumes mol is already normalized
    pub fn new(mol: Molecule, pg: PointGroup) -> Self {
        let buddy = match &pg {
            PointGroup::C1 => Buddy::C1,
            PointGroup::C2 { axis } => Buddy::C2 {
                axis: mol.try_detect_buddies(&mol.rotate(180.0, axis), 1e-8),
            },
            PointGroup::Cs { plane } => Buddy::Cs {
                plane: mol.try_detect_buddies(&mol.reflect(plane), 1e-8),
            },
            PointGroup::C2v { axis, planes } => {
                let axis =
                    mol.try_detect_buddies(&mol.rotate(180.0, axis), 1e-8);
                let plane0 =
                    mol.try_detect_buddies(&mol.reflect(&planes[0]), 1e-8);
                let plane1 =
                    mol.try_detect_buddies(&mol.reflect(&planes[1]), 1e-8);
                Buddy::C2v {
                    axis,
                    plane0,
                    plane1,
                }
            }
            PointGroup::D2h { axes, planes } => {
                let mut new_axes = Vec::new();
                for axis in axes {
                    new_axes.push(
                        mol.try_detect_buddies(&mol.rotate(180.0, axis), 1e-8),
                    );
                }
                let mut new_planes = Vec::new();
                for plane in planes {
                    new_planes
                        .push(mol.try_detect_buddies(&mol.reflect(plane), 1e-8))
                }
                Buddy::D2h {
                    axes: new_axes,
                    planes: new_planes,
                }
            }
            PointGroup::C3v { axis: _, plane: _ } => todo!(),
        };
        Self {
            map: HashMap::<Vec<Key>, Target>::new(),
            pg,
            buddy,
        }
    }

    fn to_keys(mol: &Molecule) -> Vec<Key> {
        let mut ret = Vec::with_capacity(mol.atoms.len());
        for atom in &mol.atoms {
            ret.push(Key {
                atom: atom.atomic_number,
                x: (atom.x * 1e8).round() as isize,
                y: (atom.y * 1e8).round() as isize,
                z: (atom.z * 1e8).round() as isize,
            })
        }
        ret
    }

    fn get_mut(&mut self, orig: &Molecule) -> Option<&mut Target> {
        use symm::PointGroup::{C2v, Cs, C1, C2};
        // first check the original structure
        let mol = &Self::to_keys(orig);
        if self.map.contains_key(mol) {
            return Some(self.map.get_mut(mol).unwrap());
        }

        // TODO DRY this out
        match &self.pg {
            C1 => (),
            C2 { axis } => {
                // check C2 axis
                let mol = Self::to_keys(&orig.rotate(180.0, axis));
                if self.map.contains_key(&mol) {
                    return Some(self.map.get_mut(&mol).unwrap());
                }
                for buddy in self.buddy.apply(orig) {
                    // check C2 axis
                    let mol = Self::to_keys(&buddy.rotate(180.0, axis));
                    if self.map.contains_key(&mol) {
                        return Some(self.map.get_mut(&mol).unwrap());
                    }
                }
            }
            Cs { plane } => {
                // check first mirror plane
                let mol = Self::to_keys(&orig.reflect(plane));
                if self.map.contains_key(&mol) {
                    return Some(self.map.get_mut(&mol).unwrap());
                }
                for buddy in self.buddy.apply(orig) {
                    // check first mirror plane
                    let mol = Self::to_keys(&buddy.reflect(plane));
                    if self.map.contains_key(&mol) {
                        return Some(self.map.get_mut(&mol).unwrap());
                    }
                }
            }
            C2v { axis, planes } => {
                // check C2 axis
                let mol = Self::to_keys(&orig.rotate(180.0, axis));
                if self.map.contains_key(&mol) {
                    return Some(self.map.get_mut(&mol).unwrap());
                }
                // check first mirror plane
                let mol = Self::to_keys(&orig.reflect(&planes[0]));
                if self.map.contains_key(&mol) {
                    return Some(self.map.get_mut(&mol).unwrap());
                }
                // check second mirror plane
                let mol = Self::to_keys(&orig.reflect(&planes[1]));
                if self.map.contains_key(&mol) {
                    return Some(self.map.get_mut(&mol).unwrap());
                }
                for buddy in self.buddy.apply(orig) {
                    // check C2 axis
                    let mol = &buddy.rotate(180.0, axis);
                    let key = Self::to_keys(mol);
                    if self.map.contains_key(&key) {
                        return Some(self.map.get_mut(&key).unwrap());
                    }
                    // check first mirror plane
                    let mol = Self::to_keys(&buddy.reflect(&planes[0]));
                    if self.map.contains_key(&mol) {
                        return Some(self.map.get_mut(&mol).unwrap());
                    }
                    // check second mirror plane
                    let mol = Self::to_keys(&buddy.reflect(&planes[1]));
                    if self.map.contains_key(&mol) {
                        return Some(self.map.get_mut(&mol).unwrap());
                    }
                }
            }
            PointGroup::D2h { axes, planes } => {
                // check the C2 axes
                for axis in axes {
                    let mol = &orig.rotate(180.0, axis);
                    let key = Self::to_keys(mol);
                    if self.map.contains_key(&key) {
                        return Some(self.map.get_mut(&key).unwrap());
                    }
                }
                // check the mirror planes
                for plane in planes {
                    let mol = Self::to_keys(&orig.reflect(plane));
                    if self.map.contains_key(&mol) {
                        return Some(self.map.get_mut(&mol).unwrap());
                    }
                }
                for buddy in self.buddy.apply(orig) {
                    for axis in axes {
                        let mol = &buddy.rotate(180.0, axis);
                        let key = Self::to_keys(mol);
                        if self.map.contains_key(&key) {
                            return Some(self.map.get_mut(&key).unwrap());
                        }
                    }
                    // check the mirror planes
                    for plane in planes {
                        let mol = Self::to_keys(&buddy.reflect(plane));
                        if self.map.contains_key(&mol) {
                            return Some(self.map.get_mut(&mol).unwrap());
                        }
                    }
                }
            }
            PointGroup::C3v { axis: _, plane: _ } => todo!(),
        }
        None
    }

    fn insert(&mut self, key: Molecule, value: Target) -> Option<Target> {
        self.map.insert(Self::to_keys(&key), value)
    }

    pub fn values(&self) -> Values<Vec<Key>, Target> {
        self.map.values()
    }

    pub fn len(&self) -> usize {
        self.map.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

impl Cart {
    #[allow(clippy::too_many_arguments)]
    pub fn build_points(
        &self,
        geom: Geom,
        step_size: f64,
        ref_energy: f64,
        nfc2: usize,
        nfc3: usize,
        fcs: &mut [f64],
        map: &mut BigHash,
    ) -> Vec<CartGeom> {
        let atoms = geom.xyz().unwrap();
        let (names, coords) = atom_parts(atoms);
        let ncoords = coords.len();
        let coords = na::DVector::from(coords);

        let mut geoms = Vec::new();
        // counter is the index into the the energies array that the jobs will
        // be run into
        let mut counter = 0;
        // start at 1 so that k = l = 0 indicates second derivative
        for i in 1..=ncoords {
            for j in 1..=i {
                for k in 0..=j {
                    for l in 0..=k {
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
                                    result.indices.push(Index {
                                        index,
                                        coeff: p.coeff,
                                    });
                                } else {
                                    // otherwise, mark the structure as seen by
                                    // inserting it into the map and return a
                                    // job to be run
                                    let i = map.insert(
                                        mol,
                                        Target {
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
}

pub fn freqs(dir: &str, spectro: &Spectro, mol: &Molecule) -> Output {
    let mut spectro = spectro.clone();
    let mut mol = mol.clone();
    spectro.geom = {
        mol.to_bohr();
        mol
    };

    // write input
    let input = format!("{}/spectro.in", dir);
    spectro.write(&input).unwrap();

    let spectro = Spectro::load(&input);

    let dir = Path::new(dir);
    let (output, _) = spectro.run(
        dir.join("fort.15"),
        dir.join("fort.30"),
        dir.join("fort.40"),
    );
    output
}

/// compute the index in the force constant array
fn index(n: usize, a: usize, b: usize, c: usize, d: usize) -> usize {
    match (c, d) {
        (0, 0) => intder::fc2_index(n, a, b),
        (_, 0) => intder::fc3_index(a, b, c),
        (_, _) => intder::fc4_index(a, b, c, d),
    }
}

#[derive(Clone, Debug, PartialEq)]
struct Index {
    index: usize,
    coeff: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Target {
    /// into the energy array drain is called on
    pub source_index: usize,

    /// index into the fc array with a coefficient
    indices: Vec<Index>,
}

impl<W: io::Write, Q: Queue<Mopac>> CoordType<W, Q> for Cart {
    fn run(
        &self,
        w: &mut W,
        queue: &Q,
        config: &Config,
        spectro: &Spectro,
    ) -> Output {
        let template = Template::from(&config.template);
        let (geom, ref_energy) = if config.optimize {
            let res = optimize(
                queue,
                config.geometry.clone(),
                template.clone(),
                config.charge,
            )
            .expect("optimization failed");
            (Geom::Xyz(res.cart_geom), res.energy)
        } else {
            let ref_energy = ref_energy(
                queue,
                config.geometry.clone(),
                template.clone(),
                config.charge,
            );
            (config.geometry.clone(), ref_energy)
        };

        let geom = geom.xyz().expect("expected an XYZ geometry, not Zmat");
        // 3 * #atoms
        let n = 3 * geom.len();
        let nfc2 = n * n;
        let nfc3 = n * (n + 1) * (n + 2) / 6;
        let nfc4 = n * (n + 1) * (n + 2) * (n + 3) / 24;
        let mut fcs = vec![0.0; nfc2 + nfc3 + nfc4];

        let mut mol = Molecule::new(geom.to_vec());
        mol.normalize();
        let pg = mol.point_group();
        writeln!(w, "normalized geometry:\n{}", mol).unwrap();
        let mut target_map = BigHash::new(mol.clone(), pg);

        let geoms = self.build_points(
            Geom::Xyz(mol.atoms.clone()),
            config.step_size,
            ref_energy,
            nfc2,
            nfc3,
            &mut fcs,
            &mut target_map,
        );

        let mut jobs = {
            let dir = "pts";
            let mut jobs = Vec::new();
            for (job_num, mol) in geoms.into_iter().enumerate() {
                let filename = format!("{dir}/job.{:08}", job_num);
                let job = Job::new(
                    Mopac::new(
                        filename,
                        None,
                        Rc::new(mol.geom),
                        config.charge,
                        template.clone(),
                    ),
                    mol.index,
                );
                jobs.push(job);
            }
            jobs
        };
        writeln!(
            w,
            "{n} Cartesian coordinates requires {} points",
            jobs.len()
        )
        .unwrap();

        // drain into energies
        let mut energies = vec![0.0; jobs.len()];
        queue
            .drain(&mut jobs, &mut energies)
            .expect("single-point calculations failed");

        make_fcs(&mut target_map, &energies, &mut fcs, n, nfc2, nfc3, "freqs");

        freqs("freqs", spectro, &mol)
    }
}

/// use `target_map` to symmetrize the force constants and dump them to the
/// files expected by spectro in `dir`
pub fn make_fcs(
    target_map: &mut BigHash,
    energies: &[f64],
    fcs: &mut [f64],
    n: usize,
    nfc2: usize,
    nfc3: usize,
    dir: &str,
) {
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
    intder::Intder::dump_fcs(
        dir,
        &fc2,
        &fcs[nfc2..nfc2 + nfc3],
        &fcs[nfc2 + nfc3..],
    );
}

extern crate test;
#[bench]
fn bench_to_keys(b: &mut test::Bencher) {
    use std::str::FromStr;
    let mol = Molecule::from_str(
        "
    C        0.000000   -0.888844    0.000000
    C       -0.662697    0.368254    0.000000
    C        0.662697    0.368254    0.000000
    H       -1.595193    0.906925    0.000000
    H        1.595193    0.906925    0.000000
",
    )
    .unwrap();
    b.iter(|| BigHash::to_keys(&mol))
}
