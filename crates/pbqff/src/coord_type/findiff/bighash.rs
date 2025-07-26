use std::{fmt::Display, num::ParseIntError, str::FromStr};

use serde::{Deserialize, Serialize};
use symm::{Axis, Plane, PointGroup};

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
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
pub(crate) struct Buddy {
    axes: Vec<Vec<usize>>,
    planes: Vec<Vec<usize>>,
}

impl Buddy {
    /// apply a Buddy to a Molecule, generating all of the symmetry-equivalent
    /// molecules
    pub(crate) fn apply(&self, mol: &Molecule) -> Vec<Molecule> {
        let Buddy { axes, planes } = self;
        let mut ret = Vec::new();
        for axis in axes {
            ret.push(Molecule::new(
                axis.iter().map(|i| mol.atoms[*i]).collect(),
            ));
        }
        for plane in planes {
            ret.push(Molecule::new(
                plane.iter().map(|i| mol.atoms[*i]).collect(),
            ));
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

// TODO combine this with check_axes?

/// like [check_axes], but instead of checking multiple axes for a single
/// rotation, check a single axis for multiple rotations by stepping over the
/// range 1..order. this is most useful for order > 2
fn check_degs(
    mol: &Molecule,
    axis: &Axis,
    order: usize,
    eps: f64,
) -> Vec<Vec<usize>> {
    let t = 360.0 / order as f64;
    (1..order)
        .flat_map(|d| {
            mol.try_detect_buddies(&mol.rotate(t * d as f64, axis), eps)
        })
        .collect()
}

// using impl IntoIterator so you can pass a single axis/plane as Some(...)
// instead of having to allocate a vector

/// try to generate [Buddy]s for rotations of `deg` around each [Axis] in `axes`
fn check_axes<'a>(
    mol: &Molecule,
    axes: impl IntoIterator<Item = &'a Axis>,
    deg: f64,
    eps: f64,
) -> Vec<Vec<usize>> {
    axes.into_iter()
        .flat_map(|axis| mol.try_detect_buddies(&mol.rotate(deg, axis), eps))
        .collect()
}

fn check_planes<'a>(
    mol: &Molecule,
    planes: impl IntoIterator<Item = &'a Plane>,
    eps: f64,
) -> Vec<Vec<usize>> {
    planes
        .into_iter()
        .flat_map(|plane| mol.try_detect_buddies(&mol.reflect(plane), eps))
        .collect()
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
        const EPS: f64 = 1e-8;
        let buddy = match &pg {
            PointGroup::C1 => Buddy::default(),
            PointGroup::C2 { axis } => Buddy {
                axes: check_axes(&mol, Some(axis), 180.0, EPS),
                planes: vec![],
            },
            PointGroup::C3 { axis } => Buddy {
                axes: check_axes(&mol, Some(axis), 120.0, EPS),
                planes: vec![],
            },
            PointGroup::Cs { plane } => Buddy {
                axes: vec![],
                planes: check_planes(&mol, Some(plane), EPS),
            },
            PointGroup::C2v { axis, planes } => Buddy {
                axes: check_axes(&mol, Some(axis), 180.0, EPS),
                planes: check_planes(&mol, planes, EPS),
            },
            PointGroup::D2h { axes, planes } => Buddy {
                axes: check_axes(&mol, axes, 180.0, EPS),
                planes: check_planes(&mol, planes, EPS),
            },
            PointGroup::C3v { axis, plane } => Buddy {
                axes: check_degs(&mol, axis, 3, EPS),
                planes: check_planes(&mol, Some(plane), EPS),
            },

            PointGroup::C5v { axis, plane } => Buddy {
                axes: check_degs(&mol, axis, 5, EPS),
                planes: check_planes(&mol, Some(plane), EPS),
            },
            PointGroup::D3h { c3, c2, sh, sv } => {
                let mut axes = check_degs(&mol, c3, 3, EPS);
                axes.extend(check_degs(&mol, c2, 2, EPS));
                let planes = check_planes(&mol, [sh, sv], EPS);
                Buddy { axes, planes }
            }
            PointGroup::D5h { c5, c2, sh, sv } => {
                let mut axes = check_degs(&mol, c5, 5, EPS);
                axes.extend(check_degs(&mol, c2, 2, EPS));
                let planes = check_planes(&mol, [sh, sv], EPS);
                Buddy { axes, planes }
            }
            PointGroup::C2h { axis, plane } => Buddy {
                axes: check_axes(&mol, [axis], 180.0, EPS),
                planes: check_planes(&mol, [plane], EPS),
            },
            PointGroup::C6h { c6, sh } => {
                let mut axes = check_degs(&mol, c6, 6, EPS);
                axes.extend(check_degs(&mol, c6, 2, EPS));
                axes.extend(check_degs(&mol, c6, 3, EPS));
                let planes = check_planes(&mol, [sh], EPS);
                Buddy { axes, planes }
            }
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
        use symm::PointGroup::{C1, C2, C2v, Cs};
        // first check the original structure
        let mol = &Self::to_keys(orig);
        if self.map.contains_key(mol) {
            return Some(self.map.get_mut(mol).unwrap());
        }

        macro_rules! dnh {
            ($self:expr, $orig:expr, $pairs:expr, $planes:expr) => {
                let buddies = $self.buddy.apply($orig);
                for buddy in Some($orig).into_iter().chain(buddies.iter()) {
                    for (axis, order) in $pairs {
                        let deg = 360.0 / order as f64;
                        for d in 1..order {
                            let mol = &buddy.rotate(deg * d as f64, axis);
                            let key = Self::to_keys(mol);
                            if self.map.contains_key(&key) {
                                return Some($self.map.get_mut(&key).unwrap());
                            }
                        }
                    }
                    // check the mirror planes
                    for plane in $planes {
                        let mol = Self::to_keys(&buddy.reflect(plane));
                        if self.map.contains_key(&mol) {
                            return Some($self.map.get_mut(&mol).unwrap());
                        }
                    }
                }
            };
        }

        /// helper macro for Cₙᵥ point groups
        macro_rules! cnv {
            ($self:expr, $orig:expr, $axis:expr, $order:expr, $planes:expr) => {
                dnh!($self, $orig, [($axis, $order)], $planes);
            };
        }

        match &self.pg {
            C1 => (),
            C2 { axis } => {
                cnv!(self, orig, axis, 2, []);
            }
            PointGroup::C3 { axis } => {
                cnv!(self, orig, axis, 3, []);
            }
            Cs { plane } => {
                let buddies = self.buddy.apply(orig);
                for buddy in Some(orig).into_iter().chain(buddies.iter()) {
                    // check first mirror plane
                    let mol = Self::to_keys(&buddy.reflect(plane));
                    if self.map.contains_key(&mol) {
                        return Some(self.map.get_mut(&mol).unwrap());
                    }
                }
            }
            C2v { axis, planes } => {
                cnv!(self, orig, axis, 2, planes);
            }
            PointGroup::C3v { axis, plane } => {
                cnv!(self, orig, axis, 3, [plane]);
            }
            PointGroup::C5v { axis, plane } => {
                cnv!(self, orig, axis, 5, [plane]);
            }
            PointGroup::D2h { axes, planes } => {
                dnh!(self, orig, axes.iter().zip([2, 2, 2]), planes);
            }
            PointGroup::D3h { c3, c2, sh, sv } => {
                dnh!(self, orig, [c3, c2].iter().zip([3, 2]), [sh, sv]);
            }
            PointGroup::D5h { c5, c2, sh, sv } => {
                dnh!(self, orig, [c5, c2].iter().zip([5, 2]), [sh, sv]);
            }
            PointGroup::C2h { axis, plane } => {
                cnv!(self, orig, axis, 2, [plane]);
            }
            PointGroup::C6h { c6, sh } => {
                dnh!(self, orig, [c6, c6, c6].iter().zip([6, 3, 2]), [sh]);
            }
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

    pub fn values(self) -> Vec<Target> {
        self.map.into_values().collect()
    }

    pub fn len(&self) -> usize {
        self.map.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}
