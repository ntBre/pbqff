use std::{
    fmt::Display,
    io,
    ops::{Add, AddAssign, Neg},
    str::FromStr,
};

use approx::AbsDiffEq;
use serde::{Deserialize, Serialize};

use crate::{weights::WEIGHTS, Vec3};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Atom {
    pub atomic_number: usize,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    #[serde(default)]
    pub weight: Option<f64>,
}

impl PartialEq for Atom {
    fn eq(&self, other: &Self) -> bool {
        let eps = 1e-8;
        let close = |a: f64, b: f64| (a - b).abs() < eps;
        self.atomic_number == other.atomic_number
            && close(self.x, other.x)
            && close(self.y, other.y)
            && close(self.z, other.z)
    }
}

impl AbsDiffEq for Atom {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        1e-8
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        let close = |a: f64, b: f64| (a - b).abs() < epsilon;
        self.atomic_number == other.atomic_number
            && close(self.x, other.x)
            && close(self.y, other.y)
            && close(self.z, other.z)
    }
}

impl Neg for Atom {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            atomic_number: self.atomic_number,
            x: -self.x,
            y: -self.y,
            z: -self.z,
            weight: self.weight,
        }
    }
}

impl Add<Vec3> for Atom {
    type Output = Atom;

    fn add(self, rhs: Vec3) -> Self::Output {
        Atom {
            x: self.x + rhs[0],
            y: self.y + rhs[1],
            z: self.z + rhs[2],
            ..self
        }
    }
}

impl AddAssign<Vec3> for Atom {
    fn add_assign(&mut self, rhs: Vec3) {
        *self = *self + rhs
    }
}

impl Display for Atom {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{:2} {:15.10} {:15.10} {:15.10}",
            self.label(),
            self.x,
            self.y,
            self.z
        )
    }
}

impl FromStr for Atom {
    type Err = io::Error;

    /// parse an Atom from a line like
    ///  C 1.0 1.0 1.0
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let fields: Vec<_> = s.split_whitespace().collect();
        if fields.len() != 4 {
            return Err(io::Error::other("wrong number of fields in Atom"));
        }
        let coord = fields[1..].iter().map(|s| s.parse());
        if coord.clone().any(|s| s.is_err()) {
            return Err(io::Error::other(
                "failed to parse coordinate field as f64",
            ));
        }
        let coord: Vec<_> = coord.flatten().collect();
        Ok(Self::new_from_label(
            fields[0], coord[0], coord[1], coord[2],
        ))
    }
}

pub const NUMBER_TO_SYMBOL: [&str; 55] = [
    "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
    "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
    "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe",
];

fn symbol_to_number(s: &str) -> Option<usize> {
    NUMBER_TO_SYMBOL.iter().position(|&x| x == s)
}

fn titlecase(s: &str) -> String {
    let cs: Vec<_> = s.chars().collect();
    let mut ret = String::from(cs[0]).to_uppercase();
    for c in cs.iter().skip(1) {
        ret.push_str(&c.to_lowercase().to_string());
    }
    ret
}

impl Atom {
    pub fn new(atomic_number: usize, x: f64, y: f64, z: f64) -> Self {
        Self {
            atomic_number,
            x,
            y,
            z,
            weight: None,
        }
    }

    pub fn new_from_label(atomic_symbol: &str, x: f64, y: f64, z: f64) -> Self {
        let sym = match symbol_to_number(atomic_symbol) {
            Some(s) => s,
            None => symbol_to_number(&titlecase(atomic_symbol)).unwrap_or_else(
                || panic!("failed to locate atomic symbol {atomic_symbol}"),
            ),
        };
        Self::new(sym, x, y, z)
    }

    #[inline]
    pub const fn label(&self) -> &str {
        debug_assert!(self.atomic_number != 0 && self.atomic_number < 55);
        NUMBER_TO_SYMBOL[self.atomic_number]
    }

    pub fn coord(&self) -> Vec<f64> {
        vec![self.x, self.y, self.z]
    }

    pub fn weight(&self) -> f64 {
        self.weight.unwrap_or(WEIGHTS[self.atomic_number])
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn titlecase() {
        assert_eq!(super::titlecase("AL"), "Al");
        assert_eq!(super::titlecase("Al"), "Al");
        assert_eq!(super::titlecase("al"), "Al");
        assert_eq!(super::titlecase("H"), "H");
        assert_eq!(super::titlecase("h"), "H");
    }
}
