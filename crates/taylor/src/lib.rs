// could I use const generics for these for m and n?
use anpass::Anpass;
use nalgebra as na;
use serde::{Deserialize, Serialize};
use symm::{Irrep, PointGroup};

pub use checks::*;
pub mod checks;

#[cfg(test)]
mod tests;

/// a taylor series expansion of f(x1, x2, ... n) of order m-1
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct Taylor {
    pub forces: Vec<Vec<u8>>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Disps(Vec<Vec<i8>>);

impl Disps {
    pub fn to_intder(&self, step_size: f64) -> Vec<Vec<f64>> {
        let mut ret = Vec::new();
        for disp in &self.0 {
            let disp: Vec<_> =
                disp.iter().map(|i| *i as f64 * step_size).collect();
            ret.push(disp);
        }
        ret
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    #[cfg(test)]
    pub(crate) fn sort(&mut self) {
        self.0.sort()
    }
}

impl IntoIterator for Disps {
    type Item = Vec<i8>;

    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl Taylor {
    /// generate the Taylor series mod and equivalence checks from `irreps` in
    /// `pg`
    pub fn make_checks(
        irreps: Vec<(usize, Irrep)>,
        pg: &PointGroup,
    ) -> (Option<Checks>, Option<Checks>) {
        use symm::Irrep::*;
        use symm::PointGroup::*;
        match pg {
            C1 => (None, None),
            Cs { .. } => c2_cs_checks(&irreps, App),
            C2 { .. } => c2_cs_checks(&irreps, B),
            C2v { axis: _, planes: _ } => {
                let mut checks = Checks::default();
                // first one you hit goes in checks.0, second goes in checks.1
                for i in irreps {
                    match i.1 {
                        A1 => (),
                        B1 => {
                            checks.0[1].push(i.0 + 1);
                        }
                        B2 => {
                            checks.0[0].push(i.0 + 1);
                        }
                        A2 => {
                            checks.0[2].push(i.0 + 1);
                        }
                        _ => eprintln!(
                            "warning: non-C2v irrep found in C2v point group"
                        ),
                    }
                }
                (Some(checks.clone()), Some(checks))
            }
            // if we don't know how to handle the symmetry, just print a warning
            // and don't
            _ => {
                eprintln!("don't know how to handle taylor for {pg}");
                (None, None)
            }
        }
    }

    /// returns the directly-derived Cartesian product row, where index is the
    /// desired row index, n is the truncation order and m is the number of
    /// variables in the Taylor series expansion. This corresponds to Algorithm
    /// 3 from Thackston18 with the meanings of n and m reversed to actually
    /// work
    fn row(mut index: usize, n: usize, m: usize) -> Vec<usize> {
        let mut ret = Vec::new();
        for i in (0..n).rev() {
            let ni = m.pow(i as u32);
            let di = index / ni;
            ret.push(di);
            index -= di * ni;
        }
        ret
    }

    /// takes an invalid row of the Cartesian product, the number of variables
    /// n, and the truncation order m and returns the index of the next valid
    /// row. This corresponds to Algorithm 4 in Thackston18
    fn next_row(mut row: Vec<usize>, n: usize, m: usize) -> usize {
        for i in (0..n).rev() {
            if row[i] > 0 {
                row[i] = 0;
                if i > 0 {
                    row[i - 1] += 1;
                }
                break;
            }
        }
        let mut index = 0;
        let lr = row.len();
        for i in (0..n).rev() {
            index += row[lr - i - 1] * m.pow(i as u32);
        }
        index
    }

    pub fn new(
        m: usize,
        n: usize,
        modchecks: Option<Checks>,
        eqchecks: Option<Checks>,
    ) -> Self {
        let last_index = m.pow(n as u32);
        let mut forces = Vec::new();
        let mut i = 0;
        while i < last_index {
            let row = Self::row(i, n, m);
            let s: usize = row.iter().sum();
            if s < m {
                let mc = if let Some(checks) = &modchecks {
                    checks.mod_check(&row)
                } else {
                    true
                };
                let ec = if let Some(checks) = &eqchecks {
                    checks.eq_check(&row)
                } else {
                    true
                };
                if (modchecks.is_none() && !ec)
                    || (eqchecks.is_none() && !mc)
                    || (!ec && !mc)
                {
                    i += 1;
                    continue;
                }
                forces.push(row.iter().map(|&r| r as u8).collect());
                i += 1;
            } else {
                i = Self::next_row(row, n, m);
            }
        }
        Self { forces }
    }

    /// CartProd returns the Cartesian product of the elements in prods.
    /// Implementation adapted from
    /// https://docs.python.org/3/library/itertools.html#itertools.product
    fn cart_prod(pools: Vec<Vec<i8>>) -> Vec<Vec<i8>> {
        let mut result = vec![vec![]];
        for pool in pools {
            let mut tmp = Vec::new();
            for x in result {
                for y in &pool {
                    let l = x.len() + 1;
                    let mut a = x.clone();
                    a.resize(l, 0);
                    a[l - 1] = *y;
                    tmp.push(a);
                }
            }
            result = tmp;
        }
        result
    }

    /// return the displacements associated with the expansion described by
    /// `self`
    pub fn disps(&self) -> Disps {
        let mut disps = Vec::new();
        for row in &self.forces {
            let mut indices = Vec::new();
            let mut values = Vec::new();
            for (i, digit) in row.iter().enumerate() {
                if *digit != 0 {
                    indices.push(i);
                    values.push(digit);
                }
            }
            if values.is_empty() {
                disps.push(row.iter().map(|u| *u as i8).collect());
                continue;
            }
            let mut prods = Vec::new();
            for digit in values {
                let digit = *digit as i8;
                let mut tmp = Vec::new();
                for j in (-digit..=digit).step_by(2) {
                    tmp.push(j);
                }
                prods.push(tmp);
            }
            let new_rows = Self::cart_prod(prods);
            for nrow in new_rows {
                let mut r: Vec<_> = row.iter().map(|u| *u as i8).collect();
                for (i, index) in indices.iter().enumerate() {
                    r[*index] = nrow[i];
                }
                disps.push(r);
            }
        }
        // sort -u on disps
        disps.sort();
        disps.dedup();
        Disps(disps)
    }

    pub fn to_anpass(
        &self,
        taylor_disps: &Disps,
        energies: &[f64],
        step_size: f64,
    ) -> Anpass {
        let mut disps = Vec::new();
        for disp in &taylor_disps.0 {
            for coord in disp {
                disps.push(*coord as f64 * step_size);
            }
        }
        let tdl = taylor_disps.len();
        let fl = self.forces.len();
        let mut fs = Vec::new();
        for row in &self.forces {
            for c in row {
                fs.push(*c as i32);
            }
        }
        Anpass {
            disps: na::DMatrix::from_row_slice(tdl, disps.len() / tdl, &disps),
            energies: na::DVector::from_row_slice(energies),
            exponents: na::DMatrix::from_column_slice(
                self.forces[0].len(),
                fl,
                &fs,
            ),
            bias: None,
        }
    }
}

/// helper function for generating the checks for C2 and Cs point groups.
/// `target` is the irrep to treat specially. it should be `B` for C2 and `App`
/// for Cs
fn c2_cs_checks(
    irreps: &Vec<(usize, Irrep)>,
    target: Irrep,
) -> (Option<Checks>, Option<Checks>) {
    let mut checks = Checks::default();
    for i in irreps {
        if i.1 == target {
            checks.0[0].push(i.0 + 1);
        }
    }
    (Some(checks.clone()), Some(checks))
}
