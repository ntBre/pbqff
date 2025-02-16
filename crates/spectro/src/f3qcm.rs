use std::{fmt::Display, ops::Index};

use serde::{Deserialize, Serialize};

use crate::utils::find3;

#[derive(Debug, Serialize, Deserialize)]
pub struct F3qcm(Vec<f64>);

impl F3qcm {
    pub(crate) fn with_capacity(cap: usize) -> Self {
        Self(Vec::with_capacity(cap))
    }

    pub(crate) fn push(&mut self, val: f64) {
        self.0.push(val)
    }

    pub fn new(v: Vec<f64>) -> Self {
        Self(v)
    }
}

impl Index<(usize, usize, usize)> for F3qcm {
    type Output = f64;

    fn index(&self, index: (usize, usize, usize)) -> &Self::Output {
        let (a, b, c) = index;
        let index = find3(a, b, c);
        &self.0[index]
    }
}

impl From<F3qcm> for nalgebra::DVector<f64> {
    fn from(v: F3qcm) -> Self {
        Self::from(v.0)
    }
}

impl Display for F3qcm {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let width = f.width().unwrap_or(16);
        let prec = f.precision().unwrap_or(8);
        for v in &self.0 {
            writeln!(f, "{v:width$.prec$}")?;
        }
        Ok(())
    }
}
