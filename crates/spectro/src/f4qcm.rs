use std::{
    fmt::Display,
    ops::{Index, IndexMut},
};

use serde::{Deserialize, Serialize};

use crate::utils::find4;

#[derive(Debug, Deserialize, Serialize)]
pub struct F4qcm(Vec<f64>);

impl F4qcm {
    pub fn new(v: Vec<f64>) -> Self {
        Self(v)
    }
}

#[macro_export]
macro_rules! f4qcm {
    ($elem:expr; $n:expr) => {
        F4qcm::new(vec![$elem; $n])
    };
}

impl Index<(usize, usize, usize, usize)> for F4qcm {
    type Output = f64;

    fn index(&self, index: (usize, usize, usize, usize)) -> &Self::Output {
        let (a, b, c, d) = index;
        let index = find4(a, b, c, d);
        &self.0[index]
    }
}

impl IndexMut<(usize, usize, usize, usize)> for F4qcm {
    fn index_mut(
        &mut self,
        index: (usize, usize, usize, usize),
    ) -> &mut Self::Output {
        let (a, b, c, d) = index;
        let index = find4(a, b, c, d);
        &mut self.0[index]
    }
}

impl From<F4qcm> for nalgebra::DVector<f64> {
    fn from(v: F4qcm) -> Self {
        Self::from(v.0)
    }
}

impl Display for F4qcm {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let width = f.width().unwrap_or(16);
        let prec = f.precision().unwrap_or(8);
        for v in &self.0 {
            writeln!(f, "{v:width$.prec$}")?;
        }
        Ok(())
    }
}
