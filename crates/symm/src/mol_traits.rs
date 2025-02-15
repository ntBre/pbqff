use crate::{Atom, Molecule, NUMBER_TO_SYMBOL};
use approx::AbsDiffEq;
use std::{
    collections::HashMap,
    fmt::Display,
    ops::{Add, Sub},
    str::FromStr,
    string::ParseError,
};

impl std::fmt::Debug for Molecule {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{self}")
    }
}

/// A Molecule is AbsDiffEq if each of its Atoms is
impl AbsDiffEq for Molecule {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-8
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        assert!(self.atoms.len() == other.atoms.len());
        let mut theirs = other.atoms.clone();
        if self.atoms.len() != theirs.len() {
            return false;
        }
        for atom in &self.atoms {
            let mut pops = Vec::new();
            let mut found = false;
            for (i, btom) in theirs.iter().enumerate() {
                if atom.abs_diff_eq(btom, epsilon) {
                    pops.push(i);
                    found = true;
                    break;
                }
            }
            if !found {
                return false;
            }
            // remove high indices first
            pops.sort();
            pops.reverse();
            for p in pops {
                theirs.remove(p);
            }
        }
        true
    }
}

impl PartialEq for Molecule {
    /// compare molecules for equality, irrespective of order. try to find an
    /// atom in other that equals the current atom in self. If found, remove it,
    /// so it can't be double-counted.
    fn eq(&self, other: &Self) -> bool {
        let mut theirs = other.atoms.clone();
        if self.atoms.len() != theirs.len() {
            return false;
        }
        for atom in &self.atoms {
            let mut pops = Vec::new();
            let mut found = false;
            for (i, btom) in theirs.iter().enumerate() {
                if *atom == *btom {
                    pops.push(i);
                    found = true;
                    break;
                }
            }
            if !found {
                return false;
            }
            // remove high indices first
            pops.sort();
            pops.reverse();
            for p in pops {
                theirs.remove(p);
            }
        }
        true
    }
}

impl FromStr for Molecule {
    type Err = ParseError;

    /// parse lines like
    ///      O           0.000000000    0.000000000   -0.124238453
    ///      H           0.000000000    1.431390207    0.986041184
    ///      H           0.000000000   -1.431390207    0.986041184
    /// into a molecule
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut ret = Self::default();
        let atomic_symbols: HashMap<_, _> = NUMBER_TO_SYMBOL
            .iter()
            .enumerate()
            .map(|(i, s)| (s.to_string(), i))
            .collect();
        for line in s.lines() {
            let fields = line.split_whitespace().collect::<Vec<_>>();
            if fields.len() == 4 {
                let sym = if let Some(&s) = atomic_symbols.get(fields[0]) {
                    s
                } else {
                    panic!(
                        "atomic symbol '{}' not found, tell Brent!",
                        fields[0]
                    );
                };
                ret.atoms.push(Atom::new(
                    sym,
                    fields[1].parse().unwrap(),
                    fields[2].parse().unwrap(),
                    fields[3].parse().unwrap(),
                ));
            }
        }
        Ok(ret)
    }
}

impl Add<Vec<f64>> for Molecule {
    type Output = Self;

    /// panics if the size of `rhs` doesn't align with the size of `self.atoms`
    fn add(mut self, rhs: Vec<f64>) -> Self::Output {
        if 3 * self.atoms.len() != rhs.len() {
            panic!(
                "{} atoms but {} displacements",
                self.atoms.len(),
                rhs.len()
            );
        }
        // panic above ensures rhs is exactly divisble by 3
        for (i, chunk) in rhs.chunks_exact(3).enumerate() {
            self.atoms[i].x += chunk[0];
            self.atoms[i].y += chunk[1];
            self.atoms[i].z += chunk[2];
        }
        self
    }
}

impl Sub<Molecule> for Molecule {
    type Output = Self;

    fn sub(self, rhs: Molecule) -> Self::Output {
        assert_eq!(self.atoms.len(), rhs.atoms.len());
        let mut ret = Molecule::default();
        for i in 0..self.atoms.len() {
            assert_eq!(self.atoms[i].atomic_number, rhs.atoms[i].atomic_number);
            ret.atoms.push(Atom::new(
                self.atoms[i].atomic_number,
                self.atoms[i].x - rhs.atoms[i].x,
                self.atoms[i].y - rhs.atoms[i].y,
                self.atoms[i].z - rhs.atoms[i].z,
            ));
        }
        ret
    }
}

impl Display for Molecule {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let precision = f.precision().unwrap_or(8);
        let width = f.width().unwrap_or(precision + 4);
        writeln!(f)?;
        for atom in &self.atoms {
            writeln!(
                f,
                "{:5}{:w$.p$}{:w$.p$}{:w$.p$}",
                NUMBER_TO_SYMBOL[atom.atomic_number],
                atom.x,
                atom.y,
                atom.z,
                w = width,
                p = precision,
            )?;
        }
        Ok(())
    }
}
