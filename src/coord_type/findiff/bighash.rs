use std::{
    collections::hash_map::Values, fmt::Display, num::ParseIntError,
    str::FromStr,
};

use serde::{Deserialize, Serialize};
use symm::PointGroup;

use rustc_hash::FxHashMap;

use symm::Molecule;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub(crate) struct Index {
    pub(crate) index: usize,
    pub(crate) coeff: f64,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Target {
    /// into the energy array drain is called on
    pub source_index: usize,

    /// index into the fc array with a coefficient
    pub(crate) indices: Vec<Index>,
}

#[cfg(test)]
impl approx::AbsDiffEq for Target {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.source_index == other.source_index
            && self.indices.len() == other.indices.len()
            && self.indices.iter().zip(other.indices.iter()).all(|(a, b)| {
                a.index == b.index && a.coeff.abs_diff_eq(&b.coeff, epsilon)
            })
    }
}

/// the fields are options because it's possible for `detect_buddies` to fail
/// and find that not all atoms are matched. if this happens, we want to skip
/// over those entirely instead of filling the Hash with junk geometries
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub(crate) enum Buddy {
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
    pub(crate) fn apply(&self, mol: &Molecule) -> Vec<Molecule> {
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
    pub(crate) atom: usize,
    pub(crate) x: isize,
    pub(crate) y: isize,
    pub(crate) z: isize,
}

impl Display for Key {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let Self { atom, x, y, z } = self;
        write!(f, "{atom} {x} {y} {z}")
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Hash, Deserialize)]
#[serde(from = "String")]
pub struct KeyChain(Vec<Key>);

impl From<String> for KeyChain {
    fn from(value: String) -> Self {
        Self::from_str(&value).unwrap()
    }
}

impl FromStr for KeyChain {
    type Err = ParseIntError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut ret = Vec::new();
        for line in s.lines() {
            let sp: Vec<_> = line.split_ascii_whitespace().collect();
            ret.push(Key {
                atom: sp[0].parse()?,
                x: sp[1].parse()?,
                y: sp[2].parse()?,
                z: sp[3].parse()?,
            });
        }
        Ok(Self(ret))
    }
}

impl Serialize for KeyChain {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.serialize_str(&self.to_string())
    }
}

impl Display for KeyChain {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for key in &self.0 {
            writeln!(f, "{key}")?
        }
        Ok(())
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct BigHash {
    pub(crate) map: FxHashMap<KeyChain, Target>,
    pub(crate) pg: PointGroup,
    // buddies are pairs of atoms that are interchanged across symmetry
    // operations
    pub(crate) buddy: Buddy,
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
            PointGroup::C3v { .. } => todo!(),
            PointGroup::C5v { .. } => todo!(),
            // could use c2v subgroup here
            PointGroup::D3h { .. } => todo!(),
            PointGroup::D5h { .. } => todo!(),
        };
        Self {
            map: FxHashMap::<KeyChain, Target>::default(),
            pg,
            buddy,
        }
    }

    pub(crate) fn to_keys(mol: &Molecule) -> KeyChain {
        let mut ret = Vec::with_capacity(mol.atoms.len());
        for atom in &mol.atoms {
            ret.push(Key {
                atom: atom.atomic_number,
                x: (atom.x * 1e8).round() as isize,
                y: (atom.y * 1e8).round() as isize,
                z: (atom.z * 1e8).round() as isize,
            })
        }
        KeyChain(ret)
    }

    pub(crate) fn get_mut(&mut self, orig: &Molecule) -> Option<&mut Target> {
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
            PointGroup::C3v { .. } => todo!(),
            PointGroup::C5v { .. } => todo!(),
            PointGroup::D3h { .. } => todo!(),
            PointGroup::D5h { .. } => todo!(),
        }
        None
    }

    pub(crate) fn insert(
        &mut self,
        key: Molecule,
        value: Target,
    ) -> Option<Target> {
        self.map.insert(Self::to_keys(&key), value)
    }

    pub fn values(&self) -> Values<KeyChain, Target> {
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
