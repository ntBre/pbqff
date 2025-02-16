use std::ops::{Index, IndexMut};

#[derive(Clone, Debug, Default, PartialEq)]
pub struct Checks(pub [Vec<usize>; 3]);

impl Index<(usize, usize)> for Checks {
    type Output = usize;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        &self.0[index.0][index.1]
    }
}

impl IndexMut<(usize, usize)> for Checks {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        &mut self.0[index.0][index.1]
    }
}

impl Checks {
    /// ModCheck computes a mod check of one or more subsets of digits. I'm
    /// honestly not too sure what it means, but it does something in taylor.py.
    /// Also, taylor.py takes modchecks as a dict of {2: [][]int}, so I've
    /// omitted the variable k=2 and hard-coded it since that's all we usually
    /// use.
    pub(crate) fn mod_check(&self, row: &[usize]) -> bool {
        for check in &self.0 {
            if check.iter().map(|&i| row[i - 1]).sum::<usize>() % 2 != 0 {
                return false;
            }
        }
        true
    }

    /// EqCheck computes an equivalence check of one or more subsets of digits.
    /// Not sure what this means either, but it does something in taylor.py.
    /// Like ModCheck, this takes a dict of {1: eqchecks} in the Python version,
    /// so I've ommitted the variable for the 1 since that's all we use.
    pub(crate) fn eq_check(&self, row: &[usize]) -> bool {
        for check in &self.0 {
            if check.iter().map(|&i| row[i - 1]).sum::<usize>() != 1 {
                return false;
            }
        }
        true
    }
}
