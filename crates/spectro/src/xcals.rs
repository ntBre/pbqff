//! calculate anharmonic constants for symmetric tops

type Tensor3 = tensor::tensor3::Tensor3<f64>;

use crate::{
    Dmat, Dvec, Spectro,
    f3qcm::F3qcm,
    f4qcm::F4qcm,
    ifrm1::Ifrm1,
    ifrm2::Ifrm2,
    mode::Mode,
    resonance::{Fermi1, Fermi2},
    utils::make_e0,
};

/// make the second component of E0 for symmetric tops
fn make_e2(
    modes: &[Mode],
    freq: &Dvec,
    f4qcm: &F4qcm,
    f3qcm: &F3qcm,
    ifrm1: &Ifrm1,
) -> f64 {
    let (n1dm, n2dm, _) = Mode::count(modes);
    let (i1mode, i2mode, _) = Mode::partition(modes);
    // e2
    // sss and ssss terms
    let mut f4s = 0.0;
    let mut f3s = 0.0;
    let mut f3kss = 0.0;
    for kk in 0..n2dm {
        let (k, _) = i2mode[kk];
        let wk = freq[k].powi(2);
        f4s += f4qcm[(k, k, k, k)] / 48.0;
        f3s += 11.0 * f3qcm[(k, k, k)].powi(2) / (freq[k] * 144.0);

        // kss and fss terms
        for ll in 0..n1dm {
            let l = i1mode[ll];
            let zval4 = f3qcm[(k, k, l)].powi(2);
            let wl = freq[l].powi(2);
            if ifrm1.check(k, l) {
                let delta4 = 2.0 * (2.0 * freq[k] + freq[l]);
                f3kss += zval4 / (16.0 * delta4);
            } else {
                let delta4 = 4.0 * wk - wl;
                f3kss += zval4 * freq[l] / (16.0 * delta4);
            }
        }
    }
    f4s + f3s + f3kss
}

/// make the third component of E0 for symmetric tops
fn make_e3(
    modes: &[Mode],
    freq: &Dvec,
    f3qcm: &F3qcm,
    ifrm1: &Ifrm1,
    ifrm2: &Ifrm2,
    ifrmchk: &tensor::Tensor3<usize>,
) -> f64 {
    let (n1dm, n2dm, _) = Mode::count(modes);
    let (i1mode, i2mode, _) = Mode::partition(modes);
    let mut f3kst = 0.0;
    let mut f3sst = 0.0;
    let mut f3stu = 0.0;
    // they check that ilin == 0 again here, but this all should only be called
    // if the molecule isn't linear
    for kk in 0..n2dm {
        let (k, _) = i2mode[kk];
        let wk = freq[k].powi(2);
        for ll in 0..n2dm {
            let (l, _) = i2mode[ll];
            if k == l {
                continue;
            }
            let wl = freq[l].powi(2);
            let zval5 = f3qcm[(k, k, l)].powi(2);
            if ifrm1.check(k, l) {
                // very encouraging comment from jan martin here: "this is
                // completely hosed!!!! dimension is cm-1 and should be cm !!!
                // perhaps he's multiplying by 2*FREQ(L) when he should divide?"
                // I'm assuming he fixed this and the code I've translated is
                // right.
                let delta5 =
                    (8.0 * wk + wl) / (2.0 * freq[k] + freq[l]) / (2.0 * wl);
                f3sst += zval5 * delta5 / 16.0;
            } else {
                let delta5 = (8.0 * wk + wl) / (freq[l] * (4.0 * wk - wl));
                f3sst += zval5 * delta5 / 16.0;
            }

            for mm in 0..n1dm {
                let m = i1mode[mm];
                if k <= l {
                    f3kst = 0.0;
                } else {
                    let zval6 = f3qcm[(k, l, m)].powi(2);
                    let xkst = freq[k] * freq[l] * freq[m];
                    let d1 = freq[k] + freq[l] + freq[m];
                    let d2 = freq[k] - freq[l] + freq[m];
                    let d3 = freq[k] + freq[l] - freq[m];
                    let d4 = freq[k] - freq[l] - freq[m];
                    if ifrm2.check((l, k), m) {
                        f3kst -= zval6 * (1.0 / d1 + 1.0 / d2 + 1.0 / d4) / 8.0;
                    } else if ifrm2.check((m, l), k) {
                        f3kst -= zval6 * (1.0 / d1 + 1.0 / d2 - 1.0 / d3) / 8.0;
                    } else if ifrm2.check((m, k), l) {
                        f3kst -= zval6 * (1.0 / d1 - 1.0 / d3 + 1.0 / d4) / 8.0;
                    } else {
                        f3kst -= zval6 * xkst / (2.0 * (d1 * d2 * d3 * d4));
                    }
                }
            }

            // stu term
            for mm in 0..n2dm {
                let (m, _) = i2mode[mm];
                if k <= l || l <= m {
                    f3stu = 0.0;
                } else {
                    let zval7 = f3qcm[(k, l, m)].powi(2);
                    let xstu = freq[k] * freq[l] * freq[m];
                    let d1 = freq[k] + freq[l] + freq[m];
                    let d2 = freq[k] - freq[l] + freq[m];
                    let d3 = freq[k] + freq[l] - freq[m];
                    let d4 = freq[k] - freq[l] - freq[m];
                    // TODO these should be ifrm2 as well, I think, but we
                    // haven't had any issues here yet
                    f3stu -= zval7
                        * if ifrmchk[(k, l, m)] != 0 {
                            (1.0 / d1 + 1.0 / d2 + 1.0 / d4) / 2.0
                        } else if ifrmchk[(l, m, k)] != 0 {
                            (1.0 / d1 + 1.0 / d2 - 1.0 / d3) / 2.0
                        } else if ifrmchk[(k, m, l)] != 0 {
                            (1.0 / d1 - 1.0 / d3 + 1.0 / d4) / 2.0
                        } else {
                            2.0 * xstu / (d1 * d2 * d3 * d4)
                        }
                }
            }
        }
    }
    f3kst + f3sst + f3stu
}

impl Spectro {
    /// degenerate-degenerate interactions for anharmonic constants of
    /// symmetric tops
    fn deg_deg(
        &self,
        i2mode: &Vec<(usize, usize)>,
        f4qcm: &F4qcm,
        freq: &Dvec,
        i1mode: &Vec<usize>,
        f3qcm: &F3qcm,
        ifrm1: &Ifrm1,
        xcnst: &mut Dmat,
        n2dm: usize,
        ifrm2: &Ifrm2,
        ia: usize,
        zmat: &tensor::Tensor3<f64>,
        ib: usize,
        ixyz: usize,
    ) {
        deg_deg1(i2mode, f4qcm, freq, i1mode, f3qcm, ifrm1, xcnst);

        for kk in 1..n2dm {
            let (k, k2) = i2mode[kk];
            // might be -1
            for ll in 0..kk {
                let (l, l2) = i2mode[ll];
                let val1 = (f4qcm[(k, k, l, l)] + f4qcm[(k, k, l2, l2)]) / 8.0;

                let val2: f64 = i1mode
                    .iter()
                    .map(|&m| {
                        -(f3qcm[(k, k, m)] * f3qcm[(l, l, m)] / (4.0 * freq[m]))
                    })
                    .sum();

                let valu: f64 = i1mode
                    .iter()
                    .map(|&m| {
                        let d1 = freq[k] + freq[l] + freq[m];
                        let d2 = freq[k] - freq[l] + freq[m];
                        let d3 = freq[k] + freq[l] - freq[m];
                        let d4 = -freq[k] + freq[l] + freq[m];

                        // how many more of these should be check_either?
                        if ifrm2.check_either((l, k), m) {
                            let delta = 1.0 / d1 + 1.0 / d2 + 1.0 / d4;
                            -(f3qcm[(k, l, m)].powi(2)) * delta / 16.0
                        } else if ifrm2.check((m, l), k) {
                            let delta = 1.0 / d1 + 1.0 / d2 - 1.0 / d3;
                            -(f3qcm[(k, l, m)].powi(2)) * delta / 16.0
                        } else if ifrm2.check((m, k), l) {
                            let delta = 1.0 / d1 + 1.0 / d3 + 1.0 / d4;
                            -(f3qcm[(k, l, m)].powi(2)) * delta / 16.0
                        } else {
                            let delta = -d1 * d2 * d3 * d4;
                            let val3 = freq[m].powi(2)
                                - freq[k].powi(2)
                                - freq[l].powi(2);
                            -0.25 * (f3qcm[(k, l, m)].powi(2)) * freq[m] * val3
                                / delta
                        }
                    })
                    .sum();

                let valus: f64 = i2mode
                    .iter()
                    .map(|&(m, _)| {
                        let d1 = freq[k] + freq[l] + freq[m];
                        let d2 = freq[k] - freq[l] + freq[m];
                        let d3 = freq[k] + freq[l] - freq[m];
                        let d4 = -freq[k] + freq[l] + freq[m];

                        let klm = (k, l, m);
                        if l == m && ifrm2.check((l, m), k) {
                            let delta = 8.0 * (2.0 * freq[l] + freq[k]);
                            -(f3qcm[klm].powi(2)) / delta
                        } else if k == m && ifrm2.check((k, m), l) {
                            let delta = 8.0 * (2.0 * freq[k] + freq[l]);
                            -(f3qcm[klm].powi(2)) / delta
                        } else if ifrm2.check((k, l), m) {
                            let delta = 1.0 / d1 + 1.0 / d2 + 1.0 / d4;
                            -(f3qcm[klm].powi(2)) * delta / 8.0
                        } else if ifrm2.check((l, m), k) {
                            let delta = 1.0 / d1 + 1.0 / d2 + 1.0 / d3;
                            -(f3qcm[klm].powi(2)) * delta / 8.0
                        } else if ifrm2.check((k, m), l) {
                            let delta = 1.0 / d1 + 1.0 / d3 + 1.0 / d4;
                            -(f3qcm[klm].powi(2)) * delta / 8.0
                        } else {
                            let delta = -d1 * d2 * d3 * d4;
                            let val3 = freq[m].powi(2)
                                - freq[k].powi(2)
                                - freq[l].powi(2);
                            -0.5 * (f3qcm[klm].powi(2)) * freq[m] * val3 / delta
                        }
                    })
                    .sum();

                let val5 = freq[k] / freq[l];
                let val6 = freq[l] / freq[k];
                let val7 = 0.5 * self.rotcon[ia] * zmat[(k, l2, 2)].powi(2)
                    + self.rotcon[ib] * zmat[(k, l, ixyz)].powi(2);
                let val8 = (val5 + val6) * val7;
                let value = val1 + val2 + valu + valus + val8;
                xcnst[(k, l)] = value;
                xcnst[(l, k)] = value;
                xcnst[(k2, l2)] = value;
                xcnst[(l2, k2)] = value;
            }
        }
    }

    /// g constants for degenerate modes
    fn make_gcnst(
        &self,
        n2dm: usize,
        i2mode: Vec<(usize, usize)>,
        f4qcm: &F4qcm,
        freq: &Dvec,
        i1mode: Vec<usize>,
        f3qcm: &F3qcm,
        ifrm1: Ifrm1,
        ia: usize,
        zmat: &tensor::Tensor3<f64>,
        ifrm2: &Ifrm2,
        ib: usize,
        ixyz: usize,
    ) -> Dmat {
        let mut gcnst = Dmat::zeros(self.nvib, self.nvib);

        self.gcnst1(
            n2dm, &i2mode, f4qcm, freq, &i1mode, f3qcm, &ifrm1, ia, zmat,
            &mut gcnst,
        );

        for kk in 1..n2dm {
            let (k, k2) = i2mode[kk];
            for ll in 0..kk {
                let (l, l2) = i2mode[ll];

                let valu: f64 = i1mode
                    .iter()
                    .map(|&m| {
                        let d1 = freq[k] + freq[l] + freq[m];
                        let d2 = freq[k] - freq[l] + freq[m];
                        let d3 = freq[k] + freq[l] - freq[m];
                        let d4 = -freq[k] + freq[l] + freq[m];

                        let klm = (k, l, m);
                        // TODO should more be check_either?

                        -(f3qcm[klm].powi(2))
                            * if ifrm2.check_either((l, k), m) {
                                (1.0 / d1 + 1.0 / d2 + 1.0 / d4) / 8.0
                            } else if ifrm2.check((m, l), k) {
                                (1.0 / d1 + 1.0 / d2 - 1.0 / d3) / 8.0
                            } else if ifrm2.check((m, k), l) {
                                (1.0 / d1 + 1.0 / d3 + 1.0 / d4) / 8.0
                            } else {
                                -0.5 * (freq[m] * freq[k] * freq[l])
                                    / (-d1 * d2 * d3 * d4)
                            }
                    })
                    .sum();

                let valus: f64 = i2mode
                    .iter()
                    .map(|&(m, _)| {
                        let d1 = freq[k] + freq[l] + freq[m];
                        let d2 = freq[k] - freq[l] + freq[m];
                        let d3 = freq[k] + freq[l] - freq[m];
                        let d4 = -freq[k] + freq[l] + freq[m];

                        let klm = (k, l, m);
                        -(f3qcm[klm].powi(2))
                            * if ifrm2.check((l, m), k) {
                                1.0 / (8.0 * (2.0 * freq[l] + freq[k]))
                            } else if ifrm2.check((k, m), l) {
                                1.0 / (8.0 * (2.0 * freq[k] + freq[l]))
                            } else if ifrm2.check((l, k), m) {
                                (1.0 / d1 + 1.0 / d2 + 1.0 / d4) / 8.0
                            } else if ifrm2.check((m, l), k) {
                                (1.0 / d1 + 1.0 / d2 - 1.0 / d3) / 8.0
                            } else if ifrm2.check((m, k), l) {
                                (1.0 / d1 + 1.0 / d3 + 1.0 / d4) / 8.0
                            } else {
                                (freq[m] * freq[k] * freq[l])
                                    / (-d1 * d2 * d3 * d4)
                            }
                    })
                    .sum();

                let val7 = -2.0 * self.rotcon[ib] * zmat[(k, l, ixyz)].powi(2)
                    + self.rotcon[ia] * zmat[(k, l2, 2)].powi(2)
                    + 2.0
                        * self.rotcon[ia]
                        * zmat[(k, k2, 2)]
                        * zmat[(l, l2, 2)];
                let value = valu + valus + val7;
                gcnst[(k, l)] = value;
                gcnst[(l, k)] = value;
                gcnst[(k2, l2)] = value;
                gcnst[(l2, k2)] = value;
            }
        }
        gcnst
    }

    /// the diagonal elements of G
    pub(crate) fn gcnst1(
        &self,
        n2dm: usize,
        i2mode: &[(usize, usize)],
        f4qcm: &F4qcm,
        freq: &Dvec,
        i1mode: &[usize],
        f3qcm: &F3qcm,
        ifrm1: &Ifrm1,
        ia: usize,
        zmat: &tensor::Tensor3<f64>,
        gcnst: &mut Dmat,
    ) {
        for kk in 0..n2dm {
            let (k, k2) = i2mode[kk];
            let wk = freq[k].powi(2);

            let valu: f64 = i1mode
                .iter()
                .map(|&l| {
                    let val2 = f3qcm[(k, k, l)].powi(2);
                    if ifrm1.check(k, l) {
                        let val3 = 32.0 * (2.0 * freq[k] + freq[l]);
                        val2 / val3
                    } else {
                        let wl = freq[l] * freq[l];
                        let val3 = val2 * freq[l];
                        let val4 = 16.0 * (4.0 * wk - wl);
                        -val3 / val4
                    }
                })
                .sum();

            let valus: f64 = i2mode
                .iter()
                .map(|&(m, _)| {
                    let val2 = f3qcm[(k, k, m)].powi(2);
                    if ifrm1.check(k, m) {
                        let val3 = 1.0 / (8.0 * freq[m]);
                        let val4 = 1.0 / (32.0 * (2.0 * freq[k] + freq[m]));
                        -val2 * (val3 + val4)
                    } else {
                        let wl = freq[m] * freq[m];
                        let val3 = 8.0 * wk - wl;
                        let val4 = 16.0 * freq[m] * (4.0 * wk - wl);
                        val2 * val3 / val4
                    }
                })
                .sum();

            let val7 = self.rotcon[ia] * zmat[(k, k2, 2)].powi(2);
            let value = (-f4qcm[(k, k, k, k)] / 48.0) + valu + valus + val7;
            gcnst[(k, k)] = value;
            gcnst[(k2, k2)] = value;
        }
    }

    /// nondeg-deg interactions for anharmonic constants of a symmetric top
    pub(crate) fn nondeg_deg(
        &self,
        i1mode: &Vec<usize>,
        i2mode: &Vec<(usize, usize)>,
        f4qcm: &F4qcm,
        f3qcm: &F3qcm,
        freq: &Dvec,
        ifrm2: &Ifrm2,
        ib: usize,
        zmat: &tensor::Tensor3<f64>,
        xcnst: &mut Dmat,
    ) {
        for &k in i1mode {
            for &(l, l2) in i2mode {
                let mut val2 = 0.0;
                for &m in i1mode {
                    val2 -=
                        f3qcm[(k, k, m)] * f3qcm[(l, l, m)] / (4. * freq[m]);
                }

                let valu = i2mode
                    .iter()
                    .map(|&(m, _)| {
                        let klm = (k, l, m);
                        let d1 = freq[k] + freq[l] + freq[m];
                        let d2 = freq[k] - freq[l] + freq[m];
                        let d3 = freq[k] + freq[l] - freq[m];
                        let d4 = -freq[k] + freq[l] + freq[m];

                        // ifrm2a checks mean both elements of the key are the same,
                        // really more of an ifrm1 check
                        -f3qcm[klm].powi(2)
                            * if l == m && ifrm2.check((l, m), k) {
                                1.0 / (8.0 * (2.0 * freq[l] + freq[k]))
                            } else if ifrm2.check((k, l), m) {
                                1.0 * (1.0 / d1 + 1.0 / d2 + 1.0 / d4) / 8.0
                            } else if ifrm2.check((l, m), k) {
                                1.0 * (1.0 / d1 + 1.0 / d2 + 1.0 / d3) / 8.0
                            } else if ifrm2.check((k, m), l) {
                                1.0 * (1.0 / d1 + 1.0 / d3 + 1.0 / d4) / 8.0
                            } else {
                                0.5 * freq[m]
                                    * (freq[m].powi(2)
                                        - freq[k].powi(2)
                                        - freq[l].powi(2))
                                    / (-d1 * d2 * d3 * d4)
                            }
                    })
                    .sum::<f64>();
                let val5 = freq[k] / freq[l];
                let val6 = freq[l] / freq[k];
                let val7 = self.rotcon[ib]
                    * (zmat[(k, l, 0)].powi(2) + zmat[(k, l, 1)].powi(2));
                let val8 = (val5 + val6) * val7;
                let value = (f4qcm[(k, k, l, l)] / 4.0) + val2 + valu + val8;
                xcnst[(k, l)] = value;
                xcnst[(l, k)] = value;
                xcnst[(k, l2)] = value;
                xcnst[(l2, k)] = value;
            }
        }
    }

    /// calculate the nondeg-nondeg anharmonic constants
    pub(crate) fn nondeg_nondeg(
        &self,
        i1mode: &Vec<usize>,
        f4qcm: &F4qcm,
        f3qcm: &F3qcm,
        freq: &Dvec,
        ifrm2: &Ifrm2,
        ia: usize,
        zmat: &tensor::Tensor3<f64>,
        ifrm1: &Ifrm1,
        xcnst: &mut Dmat,
    ) {
        for &k in i1mode {
            let wk = freq[k].powi(2);
            let mut valu = 0.0;
            for &l in i1mode {
                let val2 = f3qcm[(k, k, l)].powi(2);
                if ifrm1.check(k, l) {
                    let val3 = 1.0 / (8.0 * freq[l]);
                    let val4 = 1.0 / (32.0 * (2.0 * freq[k] + freq[l]));
                    valu -= val2 * (val3 + val4);
                } else {
                    let wl = freq[l].powi(2);
                    let val3 = 8.0 * wk - 3.0 * wl;
                    let val4 = 16.0 * freq[l] * (4.0 * wk - wl);
                    valu -= val2 * val3 / val4;
                }
            }
            xcnst[(k, k)] = (f4qcm[(k, k, k, k)] / 16.0) + valu;
        }

        // pretty sure it's safe to ignore this because the next loop should
        // never be entered for a diatomic since there should only be one mode
        if self.rotor.is_diatomic() {
            todo!("goto funds section. return?")
        }

        for kk in 1..i1mode.len() {
            let k = i1mode[kk];
            // might be kk-1 not sure
            for ll in 0..kk {
                let l = i1mode[ll];
                let kkll = (k, k, l, l);
                let val1 = f4qcm[kkll] / 4.0;

                let mut val2 = 0.0;
                for &m in i1mode {
                    val2 -=
                        f3qcm[(k, k, m)] * f3qcm[(l, l, m)] / (4.0 * freq[m]);
                }

                let valu = i1mode
                    .iter()
                    .map(|&m| {
                        let klm = (k, m, l);
                        let d1 = freq[k] + freq[l] + freq[m];
                        let d2 = freq[k] - freq[l] + freq[m];
                        let d3 = freq[k] + freq[l] - freq[m];
                        let d4 = -freq[k] + freq[l] + freq[m];

                        -f3qcm[klm].powi(2)
                            * if ifrm2.check((l, m), k) {
                                1.0 / (8.0 * (2.0 * freq[l] + freq[k]))
                            } else if ifrm2.check((k, m), l) {
                                1.0 / (8.0 * (2.0 * freq[k] + freq[l]))
                            } else if ifrm2.check((k, l), m) {
                                (1.0 / d1 + 1.0 / d2 + 1.0 / d4) / 8.0
                            } else if ifrm2.check((l, m), k) {
                                (1.0 / d1 + 1.0 / d2 + 1.0 / d3) / 8.0
                            } else if ifrm2.check((k, m), l) {
                                (1.0 / d1 + 1.0 / d3 + 1.0 / d4) / 8.0
                            } else {
                                0.5 * freq[m]
                                    * (freq[m].powi(2)
                                        - freq[k].powi(2)
                                        - freq[l].powi(2))
                                    / (-d1 * d2 * d3 * d4)
                            }
                    })
                    .sum::<f64>();

                let val5 = freq[k] / freq[l];
                let val6 = freq[l] / freq[k];
                let val7 = self.rotcon[ia] * zmat[(k, l, 2)].powi(2);
                let val8 = (val5 + val6) * val7;
                let value = val1 + val2 + valu + val8;
                xcnst[(k, l)] = value;
                xcnst[(l, k)] = value;
            }
        }
    }

    /// compute `ia`, `ib`, `n2dm`, `i1mode`, `i2mode`, and `ixyz` for [xcals]
    #[allow(clippy::type_complexity)]
    pub(crate) fn setup_xcals(
        &self,
        modes: &[Mode],
        wila: &Dmat,
    ) -> (usize, usize, usize, Vec<usize>, Vec<(usize, usize)>, usize) {
        let (ia, ib) = (2, 1);
        let (_, n2dm, _) = Mode::count(modes);
        let (i1mode, i2mode, _) = Mode::partition(modes);
        // find out which of a(xz)tb or a(yz)tb are zero
        let ixyz = if !self.rotor.is_linear() {
            const TOL: f64 = 0.000001;
            let mut ixz = 0;
            let mut iyz = 0;
            for &(_, i2) in &i2mode {
                if wila[(i2, 3)].abs() <= TOL {
                    ixz += 1;
                }
                if wila[(i2, 4)].abs() <= TOL {
                    iyz += 1;
                }
            }
            if ixz > 0 {
                1
            } else if iyz > 0 {
                0
            } else {
                eprintln!(
                    "warning: potential problem in xcals {}:{}:{}",
                    file!(),
                    line!(),
                    column!()
                );
                1
            }
        } else {
            1
        };
        (ia, ib, n2dm, i1mode, i2mode, ixyz)
    }

    /// calculate the anharmonic constants and E_0 for a symmetric top
    pub fn xcals(
        &self,
        f4qcm: &F4qcm,
        freq: &Dvec,
        f3qcm: &F3qcm,
        zmat: &Tensor3,
        fermi1: &[Fermi1],
        fermi2: &[Fermi2],
        modes: &[Mode],
        wila: &Dmat,
    ) -> (Dmat, Dmat, f64) {
        let (ia, ib, n2dm, i1mode, i2mode, ixyz) =
            self.setup_xcals(modes, wila);
        // NOTE skipping zeta checks, but they only print stuff

        let (ifrmchk, ifrm1, ifrm2) = self.make_fermi_checks(fermi1, fermi2);

        let e1 = make_e0(modes, f4qcm, f3qcm, freq, &ifrm1, &ifrmchk);
        let e2 = make_e2(modes, freq, f4qcm, f3qcm, &ifrm1);
        let e3 = if self.rotor.is_linear() {
            0.0
        } else {
            make_e3(modes, freq, f3qcm, &ifrm1, &ifrm2, &ifrmchk)
        };
        let e0 = e1 + e2 + e3;

        // start calculating anharmonic constants
        let mut xcnst = Dmat::zeros(self.nvib, self.nvib);

        self.nondeg_nondeg(
            &i1mode, f4qcm, f3qcm, freq, &ifrm2, ia, zmat, &ifrm1, &mut xcnst,
        );

        self.nondeg_deg(
            &i1mode, &i2mode, f4qcm, f3qcm, freq, &ifrm2, ib, zmat, &mut xcnst,
        );

        self.deg_deg(
            &i2mode, f4qcm, freq, &i1mode, f3qcm, &ifrm1, &mut xcnst, n2dm,
            &ifrm2, ia, zmat, ib, ixyz,
        );

        let gcnst = self.make_gcnst(
            n2dm, i2mode, f4qcm, freq, i1mode, f3qcm, ifrm1, ia, zmat, &ifrm2,
            ib, ixyz,
        );

        (xcnst, gcnst, e0)
    }
}

pub(crate) fn deg_deg1(
    i2mode: &Vec<(usize, usize)>,
    f4qcm: &F4qcm,
    freq: &Dvec,
    i1mode: &Vec<usize>,
    f3qcm: &F3qcm,
    ifrm1: &Ifrm1,
    xcnst: &mut Dmat,
) {
    for &(k, k2) in i2mode {
        let val1 = f4qcm[(k, k, k, k)] / 16.0;
        let wk = freq[k].powi(2);

        let mut valu = 0.0;
        for &l in i1mode {
            let val2 = f3qcm[(k, k, l)].powi(2);
            if ifrm1.check(k, l) {
                let val3 = 1.0 / (8.0 * freq[l]);
                let val4 = 1.0 / (32.0 * (2.0 * freq[k] + freq[l]));
                valu -= val2 * (val3 + val4);
            } else {
                let wl = freq[l] * freq[l];
                let val3 = 8.0 * wk - 3.0 * wl;
                let val4 = 16.0 * freq[l] * (4.0 * wk - wl);
                valu -= val2 * val3 / val4;
            }
        }

        let mut valus = 0.0;
        for &(l, _) in i2mode {
            let val2 = f3qcm[(k, k, l)].powi(2);
            if ifrm1.check(k, l) {
                let val3 = 1.0 / (8.0 * freq[l]);
                let val4 = 1.0 / (32.0 * (2.0 * freq[k] + freq[l]));
                valus -= val2 * (val3 + val4);
            } else {
                let wl = freq[l] * freq[l];
                let val3 = 8.0 * wk - 3.0 * wl;
                let val4 = 16.0 * freq[l] * (4.0 * wk - wl);
                valus -= val2 * val3 / val4;
            }
        }

        let value = val1 + valu + valus;
        xcnst[(k, k)] = value;
        xcnst[(k2, k2)] = value;
    }
}
