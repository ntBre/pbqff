use serde::{Deserialize, Serialize};
use symm::Molecule;

/// A `DummyVal` is either a literal `Value` or the index of the real `Atom` to
/// take the value from
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub(crate) enum DummyVal {
    Value(f64),
    Atom(usize),
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Dummy {
    pub(crate) x: DummyVal,
    pub(crate) y: DummyVal,
    pub(crate) z: DummyVal,
}

impl Dummy {
    pub(crate) fn get_vals(&self, mol: &Molecule) -> [f64; 3] {
        [
            match self.x {
                DummyVal::Value(v) => v,
                DummyVal::Atom(i) => mol.atoms[i].x,
            },
            match self.y {
                DummyVal::Value(v) => v,
                DummyVal::Atom(i) => mol.atoms[i].y,
            },
            match self.z {
                DummyVal::Value(v) => v,
                DummyVal::Atom(i) => mol.atoms[i].z,
            },
        ]
    }
}

impl From<Vec<DummyVal>> for Dummy {
    fn from(vals: Vec<DummyVal>) -> Self {
        assert_eq!(vals.len(), 3);
        Self {
            x: vals[0],
            y: vals[1],
            z: vals[2],
        }
    }
}
