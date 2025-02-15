use std::collections::HashMap;

pub(crate) struct Ifrm2(HashMap<(usize, usize), usize>);

impl Ifrm2 {
    pub(crate) fn new() -> Self {
        Self(HashMap::new())
    }

    pub(crate) fn insert(
        &mut self,
        k: (usize, usize),
        v: usize,
    ) -> Option<usize> {
        self.0.insert(k, v)
    }

    pub(crate) fn get(&self, k: &(usize, usize)) -> Option<&usize> {
        self.0.get(k)
    }

    /// like [get] but also check for the key with the elements in the other
    /// order
    pub(crate) fn get_either(&self, k: &(usize, usize)) -> Option<&usize> {
        if let Some(e) = self.0.get(k) {
            Some(e)
        } else {
            self.0.get(&(k.1, k.0))
        }
    }

    /// check if `k` is contained in `self` and its value is `v`
    pub(crate) fn check(&self, k: (usize, usize), v: usize) -> bool {
        let tmp = self.get(&k);
        tmp.is_some() && *tmp.unwrap() == v
    }

    /// like [check] but also check for the key with the elements in the other
    /// order
    pub(crate) fn check_either(&self, k: (usize, usize), v: usize) -> bool {
        let tmp = self.get(&k);
        if tmp.is_some() && *tmp.unwrap() == v {
            return true;
        }
        let tmp = self.get(&(k.1, k.0));
        tmp.is_some() && *tmp.unwrap() == v
    }
}

impl std::fmt::Debug for Ifrm2 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f)?;
        for ((a, b), v) in &self.0 {
            writeln!(f, "({a:5}, {b:5}) => {v:5}")?;
        }
        Ok(())
    }
}
