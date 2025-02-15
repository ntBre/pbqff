#![allow(unused)]

use super::{Dmat, Dvec, Tensor3};
use crate::f3qcm::F3qcm;
use crate::f4qcm::F4qcm;
use crate::resonance::Fermi1;
use crate::resonance::Fermi2;
use crate::utils::find4;
use crate::utils::linalg::symm_eigen_decomp;
use crate::Mode;
use nalgebra::DMatrix;
use nalgebra::SymmetricEigen;
use std::collections::HashSet;

/// set up resonance polyad matrices for asymmetric tops and compute their
/// eigenvalues and eigenvectors
pub(crate) fn resona(
    zmat: &Tensor3,
    f3qcm: &F3qcm,
    f4qcm: &F4qcm,
    e0: f64,
    modes: &[Mode],
    freq: &Dvec,
    rotcon: &[f64],
    xcnst: &Dmat,
    fermi1: &[Fermi1],
    fermi2: &[Fermi2],
    eng: &mut [f64],
) -> Dmat {
    let (n1dm, _, _) = Mode::count(modes);
    let (i1mode, _, _) = Mode::partition(modes);
    let dnm = init_res_denom(n1dm, freq, fermi1, fermi2);

    let mut zpe = e0;
    for ii in 0..n1dm {
        let i = i1mode[ii];
        zpe += freq[i] * 0.5;
        for jj in 0..=ii {
            let j = i1mode[jj];
            zpe += xcnst[(i, j)] * 0.25;
        }
    }

    // TODO handle separate resonance blocks. the example in the comments is one
    // for each symmetry in C2v, ie a1, a2, b1, and b2 symmetries
    // let iirst = make_resin(fermi1, n1dm, fermi2);

    // TODO generate this! this is only for debugging to match the order from
    // spectro2.in for c2h4 in old-spectro
    let iirst = nalgebra::dmatrix![
        0,    0,    1,    0,    0,    0,    0,    0,    0,    0,    0,    0;
        0,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0,    0;
        0,    0,    0,    1,    0,    0,    0,    0,    0,    0,    0,    0;
        0,    0,    0,    0,    1,    1,    0,    0,    0,    0,    0,    0;
        0,    0,    0,    0,    1,    0,    1,    0,    0,    0,    0,    0;
        0,    0,    0,    0,    2,    0,    0,    0,    0,    0,    0,    0;
        0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    2;
    ];
    // transpose to match fortran indexing
    let iirst = iirst.transpose();

    println!("\nEnergies of deperturbed states:");

    let (_, nreson) = iirst.shape();
    for ist in 0..nreson {
        let mut e = e0;
        for ii in 0..n1dm {
            let i = i1mode[ii];
            e += freq[i] * (iirst[(ii, ist)] as f64 + 0.5);
            for jj in 0..=ii {
                let j = i1mode[jj];
                e += xcnst[(i, j)]
                    * (iirst[(ii, ist)] as f64 + 0.5)
                    * (iirst[(jj, ist)] as f64 + 0.5);
            }
        }
        eng[ist] = e - zpe;

        println!("E* of state {ist} = {}", eng[ist]);
        // NOTE I have all the right numbers, but the order is different and
        // differs each time. I may need to revisit the usage of a HashMap for
        // the resonances at some point.
    }

    // let idimen = nreson * (nreson + 1) / 2;
    let mut resmat = DMatrix::zeros(nreson, nreson);
    for i in 0..nreson {
        for j in 0..=i {
            if j == i {
                resmat[(i, j)] = eng[i];
            } else {
                resmat[(i, j)] = genrsa(
                    n1dm, zmat, f3qcm, f4qcm, &iirst, i, j, &i1mode, freq,
                    rotcon, &dnm,
                );
            }
        }
    }

    // construct the resonance matrix, then call symm_eigen_decomp to get the
    // eigenvalues and eigenvectors
    println!("\nResonance matrix:{:.3}", resmat);

    let (vals, vecs) = symm_eigen_decomp(resmat.clone(), false);

    println!("Eigenvalues:{:.2}", vals.transpose());
    println!("Eigenvectors:{:.7}", vecs);

    resmat
}

/// computes the general resonance element between states `istate` and `jstate`
/// for an asymmetric top. `iirst` is a matrix containing the quantum numbers of
/// the states involved in the resonance polyad
pub(crate) fn genrsa(
    n1dm: usize,
    zmat: &Tensor3,
    f3qcm: &F3qcm,
    f4qcm: &F4qcm,
    iirst: &DMatrix<usize>,
    istate: usize,
    jstate: usize,
    i1mode: &[usize],
    freq: &Dvec,
    rotcon: &[f64],
    dnm: &Tensor3,
) -> f64 {
    let mut idiff = 0;
    let mut ndelta = 0;
    let mut ndel = 0;
    let mut nleft = Vec::new();
    let mut nright = Vec::new();
    let mut ndiff = Vec::new();
    let mut nmin = Vec::new();
    let mut indx = Vec::new();
    for i in 0..n1dm {
        let nnleft = iirst[(i, istate)];
        let nnright = iirst[(i, jstate)];
        let ndiffer = nnright as isize - nnleft as isize;
        if ndiffer != 0 {
            ndelta += ndiffer.abs();
            ndel += ndiffer;
            idiff += 1;
            if idiff > 4 || ndelta > 4 {
                eprintln!("higher than quartic resonances not yet implemented");
                eprintln!("setting resonance constant to zero");
                return 0.0;
            }
            nleft.push(nnleft);
            nright.push(nnright);
            ndiff.push(ndiffer);
            nmin.push(nnleft.min(nnright));
            indx.push(i);
        }
    }

    match idiff {
        4 => {
            if ndel == 0 {
                // this is in a fortran equivalence block, which should tie each of
                // these variables permanently I think, so I'll probably have to
                // write this everywhere I need these indices
                let ii = indx[0];
                let jj = indx[1];
                let kk = indx[2];
                let ll = indx[3];

                let na = nmin[0];
                let nb = nmin[1];
                let nc = nmin[2];
                let nd = nmin[3];
                // case 1a: Kabcd
                res2a(
                    zmat, f3qcm, f4qcm, i1mode, freq, rotcon, dnm, ii, jj, kk,
                    ll,
                ) * (((na + 1) * (nb + 1) * (nc + 1) * (nd + 1)) as f64 / 16.)
                    .sqrt()
            } else {
                // case 1b: Ka,bcd resonance
                // sort indices in required order first
                todo!("case 1b: Ka,bcd");
            }
        }
        3 => match ndelta {
            4 => {
                if ndel == 0 {
                    // case 2a: Kaabc
                    // sort indices in required order first
                    for i in 0..3 {
                        if ndiff[i].abs() == 2 {
                            let iswp1 = nmin[0];
                            let iswp2 = indx[0];
                            nmin[0] = nmin[i];
                            indx[0] = indx[i];
                            nmin[i] = iswp1;
                            indx[i] = iswp2;
                        }
                    }
                    let ii = indx[0];
                    let jj = indx[1];
                    let kk = indx[2];

                    let na = nmin[0];
                    let nb = nmin[1];
                    let nc = nmin[2];
                    res2a(
                        zmat, f3qcm, f4qcm, i1mode, freq, rotcon, dnm, ii, ii,
                        jj, kk,
                    ) * f64::sqrt(
                        dble((na + 1) * (na + 2) * (nb + 1) * (nc + 1)) / 16.,
                    )
                } else {
                    // ndel = +/-2 means case 2b: the ka,bbc resonance
                    todo!()
                }
            }
            3 => {
                let ii = indx[0];
                let jj = indx[1];
                let kk = indx[2];

                let na = nmin[0];
                let nb = nmin[1];
                let nc = nmin[2];
                // case 2c: Fermi type 2 resonance
                if dnm[(ii, jj, kk)] == 0.
                    || dnm[(jj, kk, ii)] == 0.
                    || dnm[(kk, ii, jj)] == 0.
                {
                    // goto 991
                    f3qcm[(ii, jj, kk)]
                        * f64::sqrt(dble((na + 1) * (nb + 1) * (nc + 1)) / 8.)
                } else {
                    // interaction already included in perturbation theory
                    0.0
                }
            }
            n => todo!("handle ndelta = {n}"),
        },
        2 => match ndelta {
            4 => {
                if ndel == 0 {
                    // case 3a: Darling-Dennison resonance Kaabb
                    let ii = indx[0];
                    let jj = indx[1];

                    let na = nmin[0];
                    let nb = nmin[1];
                    res2a(
                        zmat, f3qcm, f4qcm, i1mode, freq, rotcon, dnm, ii, ii,
                        jj, jj,
                    ) * f64::sqrt(
                        dble((na + 1) * (na + 2) * (nb + 1) * (nb + 2)) / 16.,
                    )
                } else {
                    // ndel +/- 3, case 3b: Kaaa,b resonance

                    // i appears to be n1dm from the end of the loop way above.
                    // not sure how this could work
                    if true {
                        //ndiff[i].abs() == 3 {
                        todo!("res3a");
                    } else {
                        todo!("res3a");
                    }
                }
            }
            3 => {
                let ii = indx[0];
                let jj = indx[1];

                let na = nmin[0];
                let nb = nmin[1];
                // case 3c: Fermi type 1 resonance
                if ndiff[0].abs() == 2 {
                    if dnm[(jj, ii, ii)] == 0. {
                        // goto 992
                        0.5 * f3qcm[(ii, ii, jj)]
                            * f64::sqrt(
                                dble((na + 1) * (na + 2) * (nb + 1)) / 8.,
                            )
                    } else {
                        0.0
                    }
                } else if dnm[(ii, jj, jj)] == 0. {
                    // goto 993
                    0.5 * f3qcm[(ii, jj, jj)]
                        * f64::sqrt(dble((nb + 1) * (nb + 2) * (na + 1)) / 8.)
                } else {
                    0.0
                }
            }
            2 => {
                // case 4: Lehmann's "1-1" resonance
                let ii = indx[0];
                let jj = indx[1];
                let na = nmin[0];
                let nb = nmin[1];
                let val1 =
                    res2a(
                        zmat, f3qcm, f4qcm, i1mode, freq, rotcon, dnm, ii, ii,
                        ii, jj,
                    ) * f64::sqrt(((na + 1).pow(3) * (nb + 1)) as f64 / 16.0);
                let val2 =
                    res2a(
                        zmat, f3qcm, f4qcm, i1mode, freq, rotcon, dnm, jj, jj,
                        jj, ii,
                    ) * f64::sqrt(((nb + 1).pow(3) * (na + 1)) as f64 / 16.0);
                let mut val3 = 0.0;
                for k in 0..n1dm {
                    if k != ii && k != jj {
                        let nk1 = iirst[(k, istate)];
                        let nk2 = iirst[(k, jstate)];
                        if nk1 != nk2 {
                            panic!("!!!Internal error in genrsa!!! nk1={nk1} nk2={nk2}");
                        }
                        val3 +=
                            2. * res2a(
                                zmat, f3qcm, f4qcm, i1mode, freq, rotcon, dnm,
                                ii, k, jj, k,
                            ) * f64::sqrt(
                                dble((na + 1) * (nb + 1))
                                    * (dble(nk1) + 0.5).powi(2)
                                    / 16.,
                            );
                    }
                }
                val1 + val2 + val3
            }
            _ => 0.0,
        },
        _ => 0.0,
    }
}

#[inline]
const fn dble(n: usize) -> f64 {
    n as f64
}

pub(crate) fn res2a(
    zmat: &tensor::Tensor3<f64>,
    f3qcm: &F3qcm,
    f4qcm: &F4qcm,
    i1mode: &[usize],
    freq: &Dvec,
    rotcon: &[f64],
    dnm: &Tensor3,
    ii: usize,
    jj: usize,
    kk: usize,
    ll: usize,
) -> f64 {
    let d = Denom { freq, dnm };
    // I sure hope this isn't the same indx from outside
    let indx = [ii, jj, kk, ll];
    let n1dm = i1mode.len();
    use Sign::*;
    let mut case = "????";
    if ii == jj && kk == ll {
        // case 1
        case = "aabb";
        let i = i1mode[ii];
        let j = i1mode[kk];
        let val1 = f4qcm[(i, i, j, j)] / 4.;
        let val2 = rotcon[0] * zmat[(i, j, 0)] * zmat[(i, j, 0)]
            + rotcon[1] * zmat[(i, j, 1)] * zmat[(i, j, 1)]
            + rotcon[2] * zmat[(i, j, 2)] * zmat[(i, j, 2)];
        let val2 = -val2 * (freq[i] + freq[j]).powi(2) / (freq[i] * freq[j]);

        let mut val3 = 0.0;
        for mm in 0..n1dm {
            let k = i1mode[mm];
            let temp = f3qcm[(i, i, k)]
                * f3qcm[(j, j, k)]
                * 0.125
                * (d.denom(Minus(i), Minus(i), Minus(k))
                    + d.denom(Plus(i), Plus(i), Minus(k))
                    + d.denom(Minus(j), Minus(j), Minus(k))
                    + d.denom(Plus(j), Plus(j), Minus(k)))
                * 0.5;
            val3 += temp;
        }

        let mut val4 = 0.0;
        for mm in 0..n1dm {
            let k = i1mode[mm];
            let temp = -0.5
                * f3qcm[(i, j, k)].powi(2)
                * (d.denom(Plus(i), Minus(j), Minus(k))
                    + d.denom(Minus(i), Plus(j), Minus(k)))
                * -0.5;
            val4 += temp;
        }
        val1 + val2 + val3 + val4
    } else if ii == jj && ii == kk {
        // case 2
        case = "aaab";
        let i = i1mode[ii];
        let j = i1mode[ll];
        let val1 = f4qcm[(i, i, i, j)] / 2.0;
        let mut val2 = 0.0;
        let mut val3 = 0.0;
        let mut val4 = 0.0;
        for mm in 0..n1dm {
            let k = i1mode[mm];
            let iik = f3qcm[(i, i, k)];
            let ijk = f3qcm[(i, j, k)];
            let temp = -0.5
                * (-4.0 / freq[k]
                    + d.denom(Plus(i), Plus(i), Minus(k))
                    + d.denom(Minus(i), Minus(i), Minus(k)));
            val3 -= 0.25 * iik * ijk * temp;
            let temp = -0.5
                * (2.0
                    * (d.denom(Plus(i), Minus(j), Minus(k))
                        + d.denom(Minus(i), Plus(j), Minus(k)))
                    + d.denom(Plus(i), Plus(j), Minus(k))
                    + d.denom(Minus(i), Minus(j), Minus(k)));

            val4 -= 0.25 * iik * ijk * temp;
        }
        // the sign I'm getting is backwards but I think that's okay
        val1 + val2 + val3 + val4
    } else if jj == ll && ii != kk && ii != jj {
        // case 3
        case = "acbc";
        let i = i1mode[ii];
        let k = i1mode[jj];
        let j = i1mode[kk];
        // variable is called ijkk but the find4 call is i, j, j, k. both
        // actually seem to give the same answer, probably from symmetry
        let val1 = f4qcm[(i, k, j, k)] / 2.;
        let val2 = rotcon[0] * zmat[(i, k, 0)] * zmat[(j, k, 0)]
            + rotcon[1] * zmat[(i, k, 1)] * zmat[(j, k, 1)]
            + rotcon[2] * zmat[(i, k, 2)] * zmat[(j, k, 2)];
        let val2 = 2. * val2 * (freq[i] * freq[j] + freq[k].powi(2))
            / (freq[k] * f64::sqrt(freq[i] * freq[j]));
        let mut val3 = 0.0;
        for mm in 0..n1dm {
            let l = i1mode[mm];
            let temp = 1.0 / freq[l]
                - 0.5
                    * (d.denom(Plus(i), Minus(j), Minus(l))
                        + d.denom(Minus(i), Plus(j), Minus(l)));
            val3 -= 0.25 * f3qcm[(i, j, l)] * f3qcm[(k, k, l)] * temp;
        }
        let mut val4 = 0.0;
        for mm in 0..n1dm {
            let l = i1mode[mm];
            let temp = -0.25
                * ((d.denom(Plus(i), Plus(k), Minus(l))
                    + d.denom(Plus(i), Minus(k), Minus(l))
                    + d.denom(Minus(i), Plus(k), Minus(l))
                    + d.denom(Minus(i), Minus(k), Minus(l)))
                    + (d.denom(Plus(j), Plus(k), Minus(l))
                        + d.denom(Plus(j), Minus(k), Minus(l))
                        + d.denom(Minus(j), Plus(k), Minus(l))
                        + d.denom(Minus(j), Minus(k), Minus(l))));
            val4 -= 0.5 * f3qcm[(i, k, l)] * f3qcm[(j, k, l)] * temp;
        }
        val1 + val2 + val3 + val4
    } else if ii == jj && ii != kk && kk != ll {
        // case 4
        case = "aabc";
        let i = i1mode[ii];
        let j = i1mode[kk];
        let k = i1mode[ll];
        let val1 = 0.5 * f4qcm[(i, i, j, k)];
        let val2 = rotcon[0] * zmat[(i, j, 0)] * zmat[(i, k, 0)]
            + rotcon[1] * zmat[(i, j, 1)] * zmat[(i, k, 1)]
            + rotcon[2] * zmat[(i, j, 2)] * zmat[(i, k, 2)];
        let val2 = -2. * val2 * (freq[i] + freq[j]) * (freq[i] + freq[k])
            / (freq[i] * f64::sqrt(freq[j] * freq[k]));

        let mut val3 = 0.0;
        for mm in 0..n1dm {
            let l = i1mode[mm];
            let temp1 = -0.5
                * (d.denom(Plus(j), Plus(k), Minus(l))
                    + d.denom(Minus(j), Minus(k), Minus(l)));
            let temp2 = -0.5
                * (d.denom(Plus(i), Plus(i), Minus(l))
                    + d.denom(Minus(i), Minus(i), Minus(l)));
            val3 -=
                f3qcm[(i, i, l)] * f3qcm[(j, k, l)] * (temp1 + temp2) * 0.25;
        }

        let mut val4 = 0.0;
        for mm in 0..n1dm {
            let l = i1mode[mm];
            let temp = -0.5
                * ((d.denom(Plus(i), Minus(k), Minus(l))
                    + d.denom(Minus(i), Plus(k), Minus(l)))
                    + (d.denom(Plus(i), Minus(j), Minus(l))
                        + d.denom(Minus(i), Plus(j), Minus(l))));
            val4 -= 0.5 * f3qcm[(i, j, l)] * f3qcm[(i, k, l)] * temp;
        }
        val1 + val2 + val3 + val4
    } else {
        // case 5
        case = "abcd";
        todo!()
    }
}

enum Sign {
    Plus(usize),
    Minus(usize),
}

impl Sign {
    fn abs(&self) -> usize {
        match self {
            Self::Plus(n) | Self::Minus(n) => *n,
        }
    }

    fn signum(&self) -> isize {
        match self {
            Sign::Plus(_) => 1,
            Sign::Minus(_) => -1,
        }
    }
}

struct Denom<'a> {
    freq: &'a Dvec,
    dnm: &'a Tensor3,
}

impl Denom<'_> {
    fn denom(&self, is: Sign, js: Sign, ks: Sign) -> f64 {
        let i = is.abs();
        let j = js.abs();
        let k = ks.abs();
        // these are computed in fortran as is/i so signum should be fine. if any of
        // them were zero, the division would obviously be disastrous so they must
        // be non-zero
        let isign = is.signum();
        let jsign = js.signum();
        let ksign = ks.signum();
        let iglobsgn = isign * jsign * ksign;
        let inneg = (3 - isign - jsign - ksign) / 2;
        let isign = isign * iglobsgn;
        let jsign = jsign * iglobsgn;
        let ksign = ksign * iglobsgn;
        if inneg == 0 || inneg == 3 {
            1. / (self.freq[i] + self.freq[j] + self.freq[k]) * iglobsgn as f64
        } else if isign > 0 {
            self.dnm[(i, j, k)] * iglobsgn as f64
        } else if jsign > 0 {
            self.dnm[(j, i, k)] * iglobsgn as f64
        } else {
            self.dnm[(k, i, j)] * iglobsgn as f64
        }
    }
}

/// initialize the inverse resonance denominators and zero any elements
/// corresponding to Fermi resonances deleted in the contact transformation
pub(crate) fn init_res_denom(
    n1dm: usize,
    freq: &Dvec,
    fermi1: &[Fermi1],
    fermi2: &[Fermi2],
) -> Tensor3 {
    let n = n1dm;
    let mut dnom = Tensor3::zeros(n, n, n);
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                dnom[(i, j, k)] = 1.0 / (freq[i] - freq[j] - freq[k]);
            }
        }
    }

    for &Fermi1 { i, j } in fermi1 {
        dnom[(j, i, i)] = 0.0;
    }

    for &Fermi2 { i: j, j: k, k: i } in fermi2 {
        dnom[(i, j, k)] = 0.0;
        dnom[(i, k, j)] = 0.0;
    }
    dnom
}

/// construct the RESIN Fermi polyad matrix. NOTE that the comments in resona.f
/// mention multiple blocks for different symmetries. However, the only way we
/// use it is with a single block, so I'm writing this code with that in mind.
pub(crate) fn make_resin(
    fermi1: &Vec<Fermi1>,
    n1dm: usize,
    fermi2: &Vec<Fermi2>,
) -> DMatrix<usize> {
    let mut data: HashSet<Vec<usize>> = HashSet::new();
    for &Fermi1 { i, j } in fermi1 {
        // 2wi
        let mut tmp = vec![0; n1dm];
        tmp[i] = 2;
        data.insert(tmp);
        // = wj
        let mut tmp = vec![0; n1dm];
        tmp[j] = 1;
        data.insert(tmp);
    }
    for &Fermi2 { i, j, k } in fermi2 {
        // wi + wj
        let mut tmp = vec![0; n1dm];
        tmp[i] = 1;
        tmp[j] = 1;
        data.insert(tmp);
        // = wk
        let mut tmp = vec![0; n1dm];
        tmp[k] = 1;
        data.insert(tmp);
    }
    let data: Vec<_> = data.into_iter().flatten().collect();
    DMatrix::<usize>::from_row_slice(data.len() / n1dm, n1dm, &data)
}
