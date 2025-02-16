use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum Mode {
    I1(usize),
    I2(usize, usize),
    I3(usize, usize, usize),
}

impl Mode {
    /// return the count of each type of mode in `modes`. these are referred to
    /// in the Fortran code as `n1dm`, `n2dm`, and `n3dm`
    pub fn count(modes: &[Self]) -> (usize, usize, usize) {
        let mut ret = (0, 0, 0);
        for m in modes {
            match m {
                Mode::I1(_) => ret.0 += 1,
                Mode::I2(_, _) => ret.1 += 1,
                Mode::I3(_, _, _) => ret.2 += 1,
            }
        }
        ret
    }

    /// return vectors of the separated singly-degenerate, doubly-degenerate,
    /// and triply-degenerate modes. these are referrred to in the Fortran code
    /// as `i1mode`, `i2mode`, and `i3mode`.
    #[allow(clippy::type_complexity)]
    pub fn partition(
        modes: &[Self],
    ) -> (Vec<usize>, Vec<(usize, usize)>, Vec<(usize, usize, usize)>) {
        let mut ret = (vec![], vec![], vec![]);
        for m in modes {
            match m {
                Mode::I1(i) => ret.0.push(*i),
                Mode::I2(i, j) => {
                    ret.1.push((*i, *j));
                }
                Mode::I3(i, j, k) => {
                    ret.2.push((*i, *j, *k));
                }
            }
        }
        ret
    }
}
