use std::{fmt::Display, str::FromStr};

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Params {
    pub names: Vec<String>,
    pub atoms: Vec<String>,
    pub values: Vec<f64>,
}

impl Default for Params {
    fn default() -> Self {
        Self {
            names: Default::default(),
            atoms: Default::default(),
            values: vec![0.; 0],
        }
    }
}

impl FromStr for Params {
    type Err = std::string::ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut names = Vec::new();
        let mut atoms = Vec::new();
        let mut values = Vec::new();
        for line in s.lines() {
            if !line.is_empty() {
                let fields: Vec<_> = line.split_whitespace().collect();
                names.push(fields[0].to_string());
                atoms.push(fields[1].to_string());
                values.push(fields[2].parse().unwrap());
            }
        }
        Ok(Self {
            names,
            atoms,
            values,
        })
    }
}

#[test]
fn test_from_str() {
    let got: Params = "USS            H    -11.246958000000
ZS             H      1.268641000000
BETAS          H     -8.352984000000
GSS            H     14.448686000000
USS            C    -51.089653000000
UPP            C    -39.937920000000
ZS             C      2.047558000000
ZP             C      1.702841000000
BETAS          C    -15.385236000000
BETAP          C     -7.471929000000
GSS            C     13.335519000000
GPP            C     10.778326000000
GSP            C     11.528134000000
GP2            C      9.486212000000
HSP            C      0.717322000000
FN11           C      0.046302000000"
        .parse()
        .unwrap();
    let want = Params::from_literal(
        vec![
            "USS", "ZS", "BETAS", "GSS", "USS", "UPP", "ZS", "ZP", "BETAS",
            "BETAP", "GSS", "GPP", "GSP", "GP2", "HSP", "FN11",
        ],
        vec![
            "H", "H", "H", "H", "C", "C", "C", "C", "C", "C", "C", "C", "C",
            "C", "C", "C",
        ],
        vec![
            -11.246958000000,
            1.268641000000,
            -8.352984000000,
            14.448686000000,
            -51.089653000000,
            -39.937920000000,
            2.047558000000,
            1.702841000000,
            -15.385236000000,
            -7.471929000000,
            13.335519000000,
            10.778326000000,
            11.528134000000,
            9.486212000000,
            0.717322000000,
            0.046302000000,
        ],
    );
    assert_eq!(got, want);
}

impl Display for Params {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (i, n) in self.names.iter().enumerate() {
            writeln!(f, "{} {} {:.12}", n, self.atoms[i], self.values[i])?;
        }
        Ok(())
    }
}

impl PartialEq for Params {
    fn eq(&self, other: &Self) -> bool {
        if self.names.len() != other.names.len() {
            return false;
        }
        if self.atoms.len() != other.atoms.len() {
            return false;
        }
        if self.values.len() != other.values.len() {
            return false;
        }
        for (i, n) in self.names.iter().enumerate() {
            if *n != other.names[i] {
                #[cfg(test)]
                eprintln!("{}: {} != {}", i, *n, other.names[i]);
                return false;
            }
            if self.atoms[i] != other.atoms[i] {
                #[cfg(test)]
                eprintln!("{}: {} != {}", i, self.atoms[i], other.atoms[i]);
                return false;
            }
            let diff = (self.values[i] - other.values[i]).abs();
            if diff >= 1e-12 {
                #[cfg(test)]
                eprintln!(
                    "{}: {} != {}, diff = {}",
                    i, self.values[i], other.values[i], diff
                );
                return false;
            }
        }
        true
    }
}

impl Params {
    pub fn new(
        names: Vec<String>,
        atoms: Vec<String>,
        values: Vec<f64>,
    ) -> Self {
        Self {
            names,
            atoms,
            values,
        }
    }
    pub fn from(
        names: Vec<String>,
        atoms: Vec<String>,
        values: Vec<f64>,
    ) -> Self {
        Self {
            names,
            atoms,
            values,
        }
    }

    pub fn from_literal(
        names: Vec<&str>,
        atoms: Vec<&str>,
        values: Vec<f64>,
    ) -> Self {
        Self {
            names: names.iter().map(|s| s.to_string()).collect(),
            atoms: atoms.iter().map(|s| s.to_string()).collect(),
            values,
        }
    }

    pub fn len(&self) -> usize {
        assert_eq!(self.names.len(), self.atoms.len());
        assert_eq!(self.names.len(), self.values.len());
        self.names.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}
