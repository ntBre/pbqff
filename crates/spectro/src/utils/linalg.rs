use std::cmp::Ordering;

use crate::Dmat;
use crate::Dvec;
use nalgebra::SymmetricEigen;

/// compute the eigen decomposition of the symmetric matrix `mat` and return
/// both the sorted eigenvalues and the corresponding eigenvectors in descending
/// order. the implementation is taken from RSP and the subroutines called
/// therein
pub fn symm_eigen_decomp(mat: Dmat, reverse: bool) -> (Dvec, Dmat) {
    let SymmetricEigen {
        eigenvectors: vecs,
        eigenvalues: vals,
    } = SymmetricEigen::new(mat);
    let mut pairs: Vec<_> = vals.iter().enumerate().collect();
    if reverse {
        pairs.sort_by(|(_, a), (_, b)| {
            b.partial_cmp(a).unwrap_or(Ordering::Equal)
        });
    } else {
        pairs.sort_by(|(_, a), (_, b)| {
            a.partial_cmp(b).unwrap_or(Ordering::Equal)
        });
    }
    let (rows, cols) = vecs.shape();
    let mut ret = Dmat::zeros(rows, cols);
    for i in 0..cols {
        ret.set_column(i, &vecs.column(pairs[i].0));
    }
    (
        Dvec::from_iterator(vals.len(), pairs.iter().map(|a| *a.1)),
        ret,
    )
}
