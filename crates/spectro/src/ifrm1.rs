use std::collections::HashMap;

#[derive(Debug)]
pub(crate) struct Ifrm1(HashMap<usize, usize>);

impl Ifrm1 {
    pub(crate) fn new() -> Self {
        Self(HashMap::new())
    }

    pub(crate) fn insert(&mut self, k: usize, v: usize) -> Option<usize> {
        self.0.insert(k, v)
    }

    pub(crate) fn get(&self, k: &usize) -> Option<&usize> {
        self.0.get(k)
    }

    /// check if `k` is contained in `self` and its value is `v`
    pub(crate) fn check(&self, k: usize, v: usize) -> bool {
        let tmp = self.get(&k);
        tmp.is_some() && *tmp.unwrap() == v
    }
}
