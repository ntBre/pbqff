//! compute the vibrationally-averaged rotational constants for symmetric tops

use std::collections::HashMap;

use crate::consts::ALPHA_CONST;
use crate::f3qcm::F3qcm;
use crate::resonance::Coriolis;
use crate::rotor::Rotor;
use crate::state::{State, StatePartition};
use crate::{Dmat, Dvec, Spectro, Tensor3, make_icorol, mode::Mode};

impl Spectro {
    pub(crate) fn alphas(
        &self,
        freq: &Dvec,
        wila: &Dmat,
        zmat: &Tensor3,
        f3qcm: &F3qcm,
        modes: &[Mode],
        states: &[State],
        coriolis: &[Coriolis],
    ) -> Dmat {
        if let Rotor::SphericalTop = self.rotor {
            todo!();
        }
        let (ia, ib, ic, iaia, iaib, ibib, ibic) =
            if let Rotor::ProlateSymmTop = self.rotor {
                (2, 1, 0, 5, 3, 0, 1)
            } else {
                (2, 1, 0, 5, 4, 2, 1)
            };
        let icorol = make_icorol(coriolis);
        let (n1dm, n2dm, _) = Mode::count(modes);
        let (i1mode, i2mode, _) = Mode::partition(modes);
        let alpha = self.make_alpha(
            n1dm, &i1mode, ia, freq, wila, iaia, icorol, zmat, f3qcm, n2dm,
            &i2mode, iaib, ib, ibib, ibic, ic,
        );

        // vibrationally averaged rotational constants

        // use nstop again to compute only the fundamentals, not every state
        let nstop = states.len();
        let mut rotnst = Dmat::zeros(nstop, 3);
        // only one axis for linear molecules
        let axes = if self.rotor.is_linear() {
            vec![ib]
        } else {
            vec![ia, ib]
        };
        let StatePartition { i1sts, i2sts, .. } = State::partition(states);
        for x in axes {
            for n in 0..nstop {
                let mut suma = 0.0;
                for ii in 0..n1dm {
                    let i = i1mode[ii];
                    suma += alpha[(i, x)] * (i1sts[n][ii] as f64 + 0.5);
                }

                for ii in 0..n2dm {
                    let i = i2mode[ii].0;
                    suma += alpha[(i, x)] * (i2sts[n][ii].0 as f64 + 1.0);
                }
                rotnst[(n, x)] = self.rotcon[x] + suma;
            }
        }
        rotnst
    }

    /// make the alpha matrix for a symmetric top
    pub(crate) fn make_alpha(
        &self,
        n1dm: usize,
        i1mode: &[usize],
        ia: usize,
        freq: &Dvec,
        wila: &Dmat,
        iaia: usize,
        icorol: HashMap<(usize, usize), usize>,
        zmat: &tensor::Tensor3<f64>,
        f3qcm: &F3qcm,
        n2dm: usize,
        i2mode: &[(usize, usize)],
        iaib: usize,
        ib: usize,
        ibib: usize,
        ibic: usize,
        ic: usize,
    ) -> Dmat {
        let mut alpha = Dmat::zeros(self.nvib, 3);
        // alpha A
        for kk in 0..n1dm {
            let k = i1mode[kk];
            let valu0 = 2.0 * self.rotcon[ia].powi(2) / freq[k];
            let valu1 = 0.75 * wila[(k, iaia)].powi(2) / self.primat[ia];

            let mut valu2 = 0.0;
            let mut valu3 = 0.0;
            for jj in 0..n1dm {
                let j = i1mode[jj];
                let wj32 = freq[j].powf(1.5);
                valu3 += wila[(j, iaia)] * f3qcm[(j, k, k)] * freq[k] / wj32;
                if j == k {
                    continue;
                }
                let wksq = freq[k].powi(2);
                let wjsq = freq[j].powi(2);
                let tmp = icorol.get(&(j, k));
                if tmp.is_some() && *tmp.unwrap() == ia {
                    valu2 -= 0.5
                        * zmat[(j, k, ia)].powi(2)
                        * (freq[k] - freq[j]).powi(2)
                        / (freq[j] * (freq[k] + freq[j]));
                } else {
                    valu2 += zmat[(j, k, ia)].powi(2) * (3.0 * wksq + wjsq)
                        / (wksq - wjsq);
                }
            }
            alpha[(k, ia)] = valu0 * (valu1 + valu2 + ALPHA_CONST * valu3);
        }

        for kk in 0..n2dm {
            let k = i2mode[kk].0;
            let valu0 = 2.0 * self.rotcon[ia].powi(2) / freq[k];
            let valu1 = 0.75 * wila[(k, iaib)].powi(2) / self.primat[ib];

            let mut valu2 = 0.0;
            let mut valu3 = 0.0;
            for jj in 0..n2dm {
                let (j1, j2) = i2mode[jj];
                if jj == kk {
                    continue;
                }
                let wksq = freq[k].powi(2);
                let wjsq = freq[j1].powi(2);
                let tmp = icorol.get(&(j1, k));
                if tmp.is_some() && *tmp.unwrap() == ia {
                    valu2 -= 0.5
                        * zmat[(k, j2, ia)].powi(2)
                        * (freq[k] - freq[j1]).powi(2)
                        / (freq[j1] * (freq[k] + freq[j1]));
                } else {
                    valu2 += zmat[(k, j2, ia)].powi(2) * (3.0 * wksq + wjsq)
                        / (wksq - wjsq);
                }
            }
            for jj in 0..n1dm {
                let j = i1mode[jj];
                let wj32 = freq[j].powf(1.5);
                valu3 += wila[(j, iaia)] * f3qcm[(j, k, k)] * freq[k] / wj32;
            }
            alpha[(k, ia)] = valu0 * (valu1 + valu2 + ALPHA_CONST * valu3);
        }

        // alpha B - sadly different enough that I can't loop
        for kk in 0..n1dm {
            let k = i1mode[kk];
            let valu0 = 2.0 * self.rotcon[ib].powi(2) / freq[k];
            let valu1 = 0.75
                * (wila[(k, ibib)].powi(2) + wila[(k, ibic)].powi(2))
                / self.primat[ib];

            let mut valu2 = 0.0;
            for jj in 0..n2dm {
                let j = i2mode[jj].0;
                let wjsq = freq[j].powi(2);
                let wksq = freq[k].powi(2);
                let valu3 = (3.0 * wksq + wjsq) / (wksq - wjsq);
                let tmp = icorol.get(&(j, k));
                if tmp.is_some() && *tmp.unwrap() == ic {
                    valu2 -= 0.5
                        * zmat[(k, j, ic)].powi(2)
                        * (freq[k] - freq[j]).powi(2)
                        / (freq[j] * (freq[k] + freq[j]));
                } else {
                    valu2 += zmat[(k, j, ic)].powi(2) * valu3;
                }
                if tmp.is_some() && *tmp.unwrap() == ib {
                    valu2 -= 0.5
                        * zmat[(k, j, ib)].powi(2)
                        * (freq[k] - freq[j]).powi(2)
                        / (freq[j] * (freq[k] + freq[j]));
                } else {
                    valu2 += zmat[(k, j, ib)].powi(2) * valu3;
                }
            }

            let mut valu4 = 0.0;
            for jj in 0..n1dm {
                let j = i1mode[jj];
                let wj32 = freq[j].powf(1.5);
                valu4 += f3qcm[(j, k, k)] * wila[(j, ibib)] * freq[k] / wj32;
            }
            alpha[(k, ib)] = valu0 * (valu1 + valu2 + ALPHA_CONST * valu4);
        }

        for kk in 0..n2dm {
            let (k, _) = i2mode[kk];
            let valu0 = 2.0 * self.rotcon[ib].powi(2) / freq[k];
            // TODO do diatomics count as linear?
            let (valu1, valu2) = if let Rotor::Linear = self.rotor {
                (0.0, 0.0)
            } else {
                let valu1 = 0.375
                    * (2.0 * wila[(k, ibib)].powi(2) / self.primat[ib]
                        + wila[(k, iaib)].powi(2) / self.primat[ia]);
                let mut valu2 = 0.0;
                for jj in 0..n2dm {
                    let j = i2mode[jj].0;
                    if jj == kk {
                        continue;
                    }
                    let wjsq = freq[j].powi(2);
                    let wksq = freq[k].powi(2);
                    let valu3 = (3.0 * wksq + wjsq) / (wksq - wjsq);
                    let tmp = icorol.get(&(j, k));
                    if tmp.is_some() && *tmp.unwrap() == ic {
                        valu2 -= 0.5
                            * zmat[(k, j, ic)].powi(2)
                            * (freq[k] - freq[j]).powi(2)
                            / (freq[j] * (freq[k] + freq[j]));
                    } else {
                        valu2 += zmat[(k, j, ic)].powi(2) * valu3;
                    }
                    if tmp.is_some() && *tmp.unwrap() == ib {
                        valu2 -= 0.5
                            * zmat[(k, j, ib)].powi(2)
                            * (freq[k] - freq[j]).powi(2)
                            / (freq[j] * (freq[k] + freq[j]));
                    } else {
                        valu2 += zmat[(k, j, ib)].powi(2) * valu3;
                    }
                }
                (valu1, valu2)
            };

            let mut valu4 = 0.0;
            for jj in 0..n1dm {
                let j = i1mode[jj];
                let wjsq = freq[j].powi(2);
                let wksq = freq[k].powi(2);
                let valu3 = (3.0 * wksq + wjsq) / (wksq - wjsq);
                let tmp = icorol.get(&(j, k));
                if tmp.is_some() && *tmp.unwrap() == ic {
                    valu4 -= 0.5
                        * zmat[(k, j, ic)].powi(2)
                        * (freq[k] - freq[j]).powi(2)
                        / (freq[j] * (freq[k] + freq[j]));
                } else {
                    valu4 += zmat[(k, j, ic)].powi(2) * valu3;
                }
                if tmp.is_some() && *tmp.unwrap() == ib {
                    valu4 -= 0.5
                        * zmat[(k, j, ib)].powi(2)
                        * (freq[k] - freq[j]).powi(2)
                        / (freq[j] * (freq[k] + freq[j]));
                } else {
                    valu4 += zmat[(k, j, ib)].powi(2) * valu3;
                }
            }
            valu4 *= 0.5;

            let mut valu5 = 0.0;
            for jj in 0..n1dm {
                let j = i1mode[jj];
                let wj32 = freq[j].powf(1.5);
                valu5 += f3qcm[(j, k, k)] * wila[(j, ibib)] * freq[k] / wj32;
            }
            alpha[(k, ib)] =
                valu0 * (valu1 + valu2 + valu4 + ALPHA_CONST * valu5);
        }
        // end alpha b
        alpha
    }
}
