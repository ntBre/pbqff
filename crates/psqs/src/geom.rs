use serde::{Deserialize, Serialize};
use std::{fmt::Display, str::FromStr};
use symm::atom::Atom;

#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub enum Geom {
    Xyz(Vec<Atom>),
    Zmat(String),
}

impl Default for Geom {
    fn default() -> Self {
        Self::Xyz(Default::default())
    }
}

impl Display for Geom {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Geom::Xyz(atoms) => {
                for atom in atoms {
                    writeln!(
                        f,
                        "{:5}{:15.10}{:15.10}{:15.10}",
                        atom.label(),
                        atom.x,
                        atom.y,
                        atom.z,
                    )?
                }
            }
            Geom::Zmat(g) => write!(f, "{g}")?,
        }
        Ok(())
    }
}

impl From<symm::Molecule> for Geom {
    fn from(mol: symm::Molecule) -> Self {
        Geom::Xyz(mol.atoms)
    }
}

impl FromStr for Geom {
    type Err = std::string::ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut atoms = Vec::new();
        let mut skip = 0;
        for line in s.lines() {
            let fields = line.split_whitespace().collect::<Vec<_>>();
            if skip > 0 {
                skip -= 1;
                continue;
            } else if fields.is_empty() {
                continue;
            } else if fields.len() == 1 {
                // one field, all letters => zmat
                if fields[0].chars().all(char::is_alphabetic) {
                    return Ok(Geom::Zmat(String::from(s)));
                } else {
                    // else, start of XYZ with comment line
                    skip = 1;
                    continue;
                }
            } else {
                atoms.push(line.parse().unwrap());
            }
        }
        Ok(Geom::Xyz(atoms))
    }
}

impl Geom {
    pub fn xyz(&self) -> Option<&Vec<Atom>> {
        match &self {
            Geom::Xyz(x) => Some(x),
            Geom::Zmat(_) => None,
        }
    }
    pub fn zmat(&self) -> Option<&String> {
        match &self {
            Geom::Zmat(x) => Some(x),
            Geom::Xyz(_) => None,
        }
    }

    pub fn is_xyz(&self) -> bool {
        matches!(self, Geom::Xyz(_))
    }
    pub fn is_zmat(&self) -> bool {
        matches!(self, Geom::Zmat(_))
    }
}

pub fn geom_string(geom: &Geom) -> String {
    use std::fmt::Write;
    match geom {
        Geom::Xyz(geom) => {
            let mut ret = String::with_capacity(50 * geom.len());
            for g in geom {
                writeln!(
                    ret,
                    "{} {:.12} {:.12} {:.12}",
                    g.label(),
                    g.x,
                    g.y,
                    g.z
                )
                .unwrap();
            }
            ret
        }
        Geom::Zmat(geom) => geom.to_string(),
    }
}
