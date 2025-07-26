use std::{
    cmp::{max, min},
    fmt::Display,
};

use tensor::Tensor4;
type Tensor3 = tensor::tensor3::Tensor3<f64>;

use crate::{
    Dmat, Dvec, Spectro,
    consts::SQLAM,
    f3qcm::F3qcm,
    utils::{ioff, make_tau, princ_cart, tau_prime},
};

/// struct holding the sextic distortion constants. For an asymmetric top, the
/// field names are correct. For a symmetric top, phi->H and sphij->h1,
/// sphijk->h2, and sphik->h3. TODO make this an enum to get rid of this comment
#[derive(
    Clone, Debug, Default, PartialEq, serde::Serialize, serde::Deserialize,
)]
pub struct Sextic {
    // a reduced constants
    pub phij: f64,
    pub phijk: f64,
    pub phikj: f64,
    pub phik: f64,
    pub sphij: f64,
    pub sphijk: f64,
    pub sphik: f64,
    // s reduced constants
    pub hj: f64,
    pub hjk: f64,
    pub hkj: f64,
    pub hk: f64,
    pub h1: f64,
    pub h2: f64,
    pub h3: f64,
    // linear molecules
    pub he: f64,
}

#[cfg(test)]
macro_rules! impl_abs_diff {
    ($lhs:ident, $rhs:ident, $eps:ident, $($field:ident$(,)?),*) => {
	$(
	    $lhs.$field.abs_diff_eq(&$rhs.$field, $eps) &&
	)*
	    // so I can have the trailing and
	    true
    };
}

#[cfg(test)]
impl approx::AbsDiffEq for Sextic {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        impl_abs_diff!(
            self, other, epsilon, phij, phijk, phikj, phik, sphij, sphijk,
            sphik, hj, hjk, hkj, hk, h1, h2, h3, he,
        )
    }
}

fn scc(
    maxcor: usize,
    tau: Tensor4,
    rotcon: &[f64],
    nvib: usize,
    freq: &Dvec,
    cc: Tensor4,
    f3qcm: &F3qcm,
    c: &Dmat,
    spectro: &Spectro,
) -> Tensor3 {
    // some kind of tolerance for messing with certain values
    const TOL: f64 = 1e-4;
    let mut scc = Tensor3::zeros(3, 3, 3);
    let mut rotcon = rotcon.to_vec();
    if freq.len() == 6 {
        // proxy for is_diatomic
        rotcon.reverse();
    }
    for ixyz in 1..=maxcor {
        let x = ixyz - 1;
        let iixyz = ioff(ixyz + 1);
        let s = iixyz - 1;

        let mut val1 = 0.0;
        for y in 0..maxcor {
            val1 += tau[(y, x, x, x)].powi(2) / rotcon[y];
        }
        val1 = 3.0 * val1 / 16.0;

        let mut val2 = 0.0;
        for i in 0..nvib {
            val2 += freq[i] * cc[(i, x, x, x)].powi(2);
        }
        val2 *= 2.0;

        let mut val3 = 0.0;
        for i in 0..nvib {
            for j in 0..nvib {
                for k in 0..nvib {
                    val3 +=
                        f3qcm[(i, j, k)] * c[(i, s)] * c[(j, s)] * c[(k, s)];
                }
            }
        }
        val3 /= 6.0;

        let mut val4 = 0.0;
        if !spectro.rotor.is_linear() {
            for y in 0..maxcor {
                if x != y {
                    let div = rotcon[x] - rotcon[y];
                    if div.abs() > TOL {
                        val4 += tau[(y, x, x, x)].powi(2) / div;
                    }
                }
            }
        };
        val4 *= 0.25;

        let value = val1 - val2 + val3 + val4;
        scc[(x, x, x)] = value;
    }

    if spectro.is_linear() {
        return scc;
    }

    for ixyz in 1..=3 {
        let x = ixyz - 1;
        let s = ioff(ixyz + 1) - 1;
        for jxyz in 1..=3 {
            let y = jxyz - 1;
            if x == y {
                continue;
            }
            let t = ioff(jxyz + 1) - 1;
            let u = ioff(max(ixyz, jxyz)) + min(ixyz, jxyz) - 1;

            let mut vala = 0.0;
            for z in 0..3 {
                let val1 =
                    (tau[(z, y, x, x)] + 2.0 * tau[(z, x, y, x)]).powi(2);
                let val2 = 2.0
                    * tau[(z, x, x, x)]
                    * (tau[(z, x, y, y)] + 2.0 * tau[(z, y, y, x)]);
                vala += (val1 + val2) / rotcon[z];
            }
            vala *= 3.0 / 32.0;

            let mut valb = 0.0;
            for i in 0..nvib {
                let val1 = cc[(i, y, x, x)].powi(2)
                    + 2.0 * cc[(i, x, x, x)] * cc[(i, y, y, x)];
                valb += freq[i] * val1;
            }

            let mut valc = 0.0;
            for i in 0..nvib {
                for j in 0..nvib {
                    let val1 =
                        c[(j, s)] * c[(i, t)] + 4.0 * c[(j, u)] * c[(i, u)];
                    for k in 0..nvib {
                        valc += f3qcm[(i, j, k)] * c[(k, s)] * val1;
                    }
                }
            }
            valc /= 4.0;

            let div = 8.0 * (rotcon[x] - rotcon[y]);
            let mut vald = 0.0;
            if div.abs() > TOL {
                let val1 = 4.0 * tau[(y, y, y, x)] - 3.0 * tau[(y, x, x, x)];
                vald += val1 * tau[(y, x, x, x)] / div;
            }

            let mut vale = 0.0;
            let mut valf = 0.0;
            for kxyz in 1..=3 {
                let z = kxyz - 1;
                if z == y || z == x {
                    continue;
                }
                let div = 4.0 * (rotcon[x] - rotcon[z]).powi(2);
                if div > TOL {
                    let val1 = (rotcon[x] - rotcon[z])
                        * (tau[(z, x, y, y)] + 2.0e0 * tau[(z, y, y, x)]);
                    let val2 = (rotcon[x] - rotcon[y])
                        * (tau[(z, z, z, x)] - tau[(z, x, x, x)]);
                    vale += tau[(z, x, x, x)] * (val1 + val2) / div;
                }

                let div = 8.0 * (rotcon[y] - rotcon[z]).powi(2);
                if div > TOL {
                    let val1 = (rotcon[y] - rotcon[z])
                        * (tau[(z, y, x, x)] + 2.0 * tau[(z, x, y, x)]);
                    let val2 = 2.0
                        * (rotcon[x] - rotcon[y])
                        * (tau[(z, y, y, y)] - tau[(z, z, z, y)]);
                    let val3 = tau[(z, y, x, x)] + 2.0 * tau[(z, x, y, x)];
                    valf += val3 * (val1 + val2) / div;
                }
            }

            let value = vala - valb + valc + vald + vale + valf;
            scc[(y, x, x)] = value;
            scc[(x, y, x)] = value;
            scc[(x, x, y)] = value;
        }
    } // end loop at line 397

    let mut vala = 0.0;
    for x in 0..3 {
        vala += ((2.0
            * (tau[(x, 0, 1, 2)] + tau[(x, 1, 2, 0)] + tau[(x, 2, 0, 1)])
                .powi(2))
            + ((tau[(x, 0, 1, 1)] + 2.0 * tau[(x, 1, 0, 1)])
                * (tau[(x, 0, 2, 2)] + 2.0 * tau[(x, 2, 0, 2)]))
            + ((tau[(x, 1, 2, 2)] + 2.0 * tau[(x, 2, 1, 2)])
                * (tau[(x, 1, 0, 0)] + 2.0 * tau[(x, 0, 1, 0)]))
            + ((tau[(x, 2, 0, 0)] + 2.0 * tau[(x, 0, 2, 0)])
                * (tau[(x, 2, 1, 1)] + 2.0 * tau[(x, 1, 2, 1)])))
            / rotcon[x];
    }
    vala = 3.0 * vala / 16.0;

    let mut valb = 0.0;
    for i in 0..nvib {
        valb += freq[i]
            * (2.0 * cc[(i, 2, 1, 0)].powi(2)
                + cc[(i, 1, 1, 0)] * cc[(i, 2, 2, 0)]
                + cc[(i, 2, 2, 1)] * cc[(i, 1, 0, 0)]
                + cc[(i, 2, 0, 0)] * cc[(i, 2, 1, 1)]);
    }
    valb *= 2.0;

    let mut valc = 0.0;
    for i in 0..nvib {
        for j in 0..nvib {
            for k in 0..nvib {
                valc += f3qcm[(i, j, k)]
                    * (c[(i, 0)] * c[(j, 2)] * c[(k, 5)]
                        + 2.0 * c[(i, 0)] * c[(j, 4)] * c[(k, 4)]
                        + 2.0 * c[(i, 2)] * c[(j, 3)] * c[(k, 3)]
                        + 2.0 * c[(i, 5)] * c[(j, 1)] * c[(k, 1)]
                        + 8.0 * c[(i, 4)] * c[(j, 3)] * c[(k, 1)]);
            }
        }
    }
    valc *= 0.5;

    let mut vald = 0.0;
    let div = 4.0 * (rotcon[1] - rotcon[2]).powi(2);
    if div > TOL {
        vald = (3.0 * (tau[(2, 1, 1, 1)] - tau[(2, 2, 2, 1)]))
            * ((rotcon[2] - rotcon[0]) * tau[(2, 1, 1, 1)]
                + (rotcon[0] - rotcon[1]) * tau[(2, 2, 2, 1)])
            / div;
    }

    let mut vale = 0.0;
    let div = 4.0 * (rotcon[2] - rotcon[0]).powi(2);
    if div > TOL {
        vale = (3.0 * (tau[(2, 2, 2, 0)] - tau[(2, 0, 0, 0)]))
            * ((rotcon[0] - rotcon[1]) * tau[(2, 2, 2, 0)]
                + (rotcon[1] - rotcon[2]) * tau[(2, 0, 0, 0)])
            / div;
    }

    let mut valf = 0.0;
    let div = 4.0 * (rotcon[0] - rotcon[1]).powi(2);
    if div > TOL {
        valf = (3.0 * (tau[(1, 0, 0, 0)] - tau[(1, 1, 1, 0)]))
            * ((rotcon[1] - rotcon[2]) * tau[(1, 0, 0, 0)]
                + (rotcon[2] - rotcon[0]) * tau[(1, 1, 1, 0)])
            / div;
    }

    let value = vala - valb + valc + vald + vale + valf;

    scc[(0, 1, 2)] = value;
    scc[(0, 2, 1)] = value;
    scc[(1, 0, 2)] = value;
    scc[(1, 2, 0)] = value;
    scc[(2, 0, 1)] = value;
    scc[(2, 1, 0)] = value;
    scc
}

fn cc_tensor(
    nvib: usize,
    maxcor: usize,
    freq: &Dvec,
    c: &Dmat,
    zmat: &Tensor3,
    rotcon: &[f64],
) -> Tensor4 {
    let mut cc = Tensor4::zeros(nvib, 3, 3, 3);
    for i in 0..nvib {
        for x in 1..=maxcor {
            let s = ioff(x + 1);

            let mut val = 0.0;
            for j in 0..nvib {
                let val1 = 1.0 / freq[j].powi(2) + 2.0 / freq[i].powi(2);
                let val2 = c[(j, s - 1)]
                    * zmat[(j, i, x - 1)]
                    * rotcon[x - 1]
                    * freq[j].powf(1.5);
                val += val2 * val1;
            }
            val /= freq[i].sqrt();
            cc[(i, x - 1, x - 1, x - 1)] = val;

            for y in 1..=maxcor {
                if x == y {
                    continue;
                }
                let t = ioff(max(x, y)) + min(x, y);

                let mut val = 0.0;
                for j in 0..nvib {
                    let val1 = 1.0 / freq[j].powi(2) + 2.0 / freq[i].powi(2);
                    let val2 = freq[j].powf(1.5)
                        * (c[(j, s - 1)] * zmat[(j, i, y - 1)] * rotcon[y - 1]
                            + 2.0
                                * c[(j, t - 1)]
                                * zmat[(j, i, x - 1)]
                                * rotcon[x - 1]);
                    val += val2 * val1;
                }
                val /= freq[i].sqrt();
                cc[(i, y - 1, x - 1, x - 1)] = val;
                cc[(i, x - 1, y - 1, x - 1)] = val;
                cc[(i, x - 1, x - 1, y - 1)] = val;

                for z in 1..=maxcor {
                    if x == z || y == z {
                        continue;
                    }
                    let u = ioff(max(x, z)) + min(x, z);
                    let v = ioff(max(y, z)) + min(y, z);

                    let mut val = 0.0;
                    for j in 0..nvib {
                        let val1 =
                            1.0 / freq[j].powi(2) + 2.0 / freq[i].powi(2);
                        let val2 = freq[j].powf(1.5)
                            * (c[(j, v - 1)]
                                * zmat[(j, i, x - 1)]
                                * rotcon[x - 1]
                                + c[(j, u - 1)]
                                    * zmat[(j, i, y - 1)]
                                    * rotcon[y - 1]
                                + c[(j, t - 1)]
                                    * zmat[(j, i, z - 1)]
                                    * rotcon[z - 1]);
                        val += val2 * val1;
                    }
                    val /= freq[i].sqrt();
                    cc[(i, z - 1, y - 1, x - 1)] = val;
                }
            }
        }
    }
    cc
}

pub(crate) fn c_mat(
    maxcor: usize,
    nvib: usize,
    freq: &Dvec,
    primat: &[f64],
    wila: &Dmat,
) -> Dmat {
    let mut c = Dmat::zeros(nvib, 6);
    for i in 0..nvib {
        let cnst = (SQLAM * SQLAM * freq[i]).powf(1.5);
        for x in 0..maxcor {
            for y in 0..=x {
                let z = ioff(x + 1) + y;
                let div = 2.0 * primat[x] * primat[y] * cnst;
                c[(i, z)] = wila[(i, z)] / div;
            }
        }
    }
    c
}

impl Sextic {
    pub fn new(
        s: &Spectro,
        wila: &Dmat,
        zmat: &Tensor3,
        freq: &Dvec,
        f3qcm: &F3qcm,
    ) -> Self {
        let mut ret = Self::default();
        let nvib = s.nvib;
        let maxcor = if s.is_linear() { 2 } else { 3 };
        let primat = if s.rotor.is_diatomic() {
            let mut p = s.primat.clone();
            p.reverse();
            p
        } else {
            s.primat.clone()
        };
        let c = c_mat(maxcor, nvib, freq, &primat, wila);
        let cc = cc_tensor(nvib, maxcor, freq, &c, zmat, &s.rotcon);
        let tau = make_tau(maxcor, nvib, freq, &primat, wila);
        let taucpm = tau_prime(maxcor, &tau);
        let scc = scc(maxcor, tau, &s.rotcon, nvib, freq, cc, f3qcm, &c, s);

        if s.rotor.is_diatomic() {
            ret.he = scc[(1, 1, 1)];
            return ret;
        } else if s.rotor.is_linear() {
            ret.he = scc[(0, 0, 0)];
            return ret;
        }
        // says this is the default representation and sets it to 1 if it was
        // originally 0. NOTE my 0 is fortran 1
        let irep = if s.rotor.is_sym_top() { 5 } else { 0 };
        let (ic, id) = princ_cart(irep);
        let mut t = Dmat::zeros(maxcor, maxcor);
        for ixyz in 0..maxcor {
            for jxyz in 0..maxcor {
                t[(ic[ixyz], ic[jxyz])] = taucpm[(ixyz, jxyz)] / 4.0;
            }
        }
        let phi = make_phi(maxcor, ic, scc);
        let t400 = (3.0 * t[(0, 0)] + 3.0 * t[(1, 1)] + 2.0 * t[(0, 1)]) / 8.0;
        let t220 = t[(0, 2)] + t[(1, 2)] - 2.0 * t400;
        let t040 = t[(2, 2)] - t220 - t400;
        let t202 = (t[(0, 0)] - t[(1, 1)]) / 4.0;
        let t022 = (t[(0, 2)] - t[(1, 2)]) / 2.0 - t202;
        let t004 = (t[(0, 0)] + t[(1, 1)] - 2.0 * t[(0, 1)]) / 16.0;
        let b002 = 0.25e0 * (s.rotcon[id[0]] - s.rotcon[id[1]]);
        let phi600 = 5.0 * (phi[(0, 0, 0)] + phi[(1, 1, 1)]) / 16.0
            + (phi[(0, 0, 1)] + phi[(1, 1, 0)]) / 8.0;
        let phi420 = 3.0 * (phi[(0, 0, 2)] + phi[(1, 1, 2)]) / 4.0
            + phi[(0, 1, 2)] / 4.0
            - 3.0 * phi600;
        let phi240 =
            phi[(2, 2, 0)] + phi[(2, 2, 1)] - 2.0 * phi420 - 3.0 * phi600;
        let phi060 = phi[(2, 2, 2)] - phi240 - phi420 - phi600;
        let phi402 = 15.0 * (phi[(0, 0, 0)] - phi[(1, 1, 1)]) / 64.0
            + (phi[(0, 0, 1)] - phi[(1, 1, 0)]) / 32.0;
        let phi222 = (phi[(0, 0, 2)] - phi[(1, 1, 2)]) / 2.0 - 2.0 * phi402;
        let phi042 = (phi[(2, 2, 0)] - phi[(2, 2, 1)]) / 2.0 - phi222 - phi402;
        let phi204 = 3.0 * (phi[(0, 0, 0)] + phi[(1, 1, 1)]) / 32.0
            - (phi[(0, 0, 1)] + phi[(1, 1, 0)]) / 16.0;
        let phi024 =
            (phi[(0, 0, 2)] + phi[(1, 1, 2)] - phi[(0, 1, 2)]) / 8.0 - phi204;
        let phi006 = (phi[(0, 0, 0)] - phi[(1, 1, 1)]) / 64.0
            - (phi[(0, 0, 1)] - phi[(1, 1, 0)]) / 32.0;

        let sigma = if !s.rotor.is_sym_top() {
            // asymmetric top
            let sigma =
                (2.0 * s.rotcon[id[2]] - s.rotcon[id[0]] - s.rotcon[id[1]])
                    / (s.rotcon[id[0]] - s.rotcon[id[1]]);
            ret.phij = phi600 + 2.0 * phi204;
            ret.phijk = phi420 - 12.0 * phi204
                + 2.0 * phi024
                + 16.0 * sigma * phi006
                + 8.0 * t022 * t004 / b002;
            ret.phikj = phi240 + 10.0 * phi420 / 3.0
                - 30.0 * phi204
                - 10.0 * ret.phijk / 3.0;
            ret.phik = phi060 - 7.0 * phi420 / 3.0
                + 28.0 * phi204
                + 7.0 * ret.phijk / 3.0;
            ret.sphij = phi402 + phi006;
            ret.sphijk = phi222 + 4.0 * sigma * phi204 - 10.0 * phi006
                + 2.0 * (t220 - 2.0 * sigma * t202 - 4.0 * t004) * t004 / b002;
            ret.sphik = phi042
                + 4.0 * sigma * phi024 / 3.0
                + (32.0 * sigma * sigma / 3.0 + 9.0) * phi006
                + 4.0
                    * (t040 + sigma * t022 / 3.0
                        - 2.0 * (sigma * sigma - 2.0) * t004)
                    * t004
                    / b002;
            sigma
        } else {
            let irep = if s.rotor.is_prolate() { 0 } else { 5 };
            let (_, id) = princ_cart(irep);
            let rkappa =
                (2.0 * s.rotcon[id[0]] - s.rotcon[id[1]] - s.rotcon[id[2]])
                    / (s.rotcon[id[1]] - s.rotcon[id[2]]);

            if rkappa < 0.0 {
                -999999999999.0
            } else {
                let bo = (rkappa - 1.0) / (rkappa + 3.0);
                1.0 / bo
            }
        };

        let irep = if s.rotor.is_prolate() || s.rotor.is_asymm_top() {
            0
        } else {
            5
        };
        let (_, id) = princ_cart(irep);

        let b200 = 0.5 * (s.rotcon[id[0]] + s.rotcon[id[1]]) - 4.0 * t004;
        let b020 = s.rotcon[id[2]] - b200 + 6.0 * t004;

        let div = 2.0 * sigma * sigma + 27.0 / 16.0;
        let mu = (sigma * phi042 - 9.0 * phi024 / 8.0
            + (-2.0 * sigma * t040 + (sigma * sigma + 3.0) * t022
                - 5.0 * sigma * t004)
                * t022
                / b020)
            / div;
        let nu = 3.0 * mu / (16.0 * sigma)
            + phi024 / (8.0 * sigma)
            + t004 * t022 / b020;
        let lamda = 5.0 * nu / sigma
            + phi222 / (sigma * 2.0)
            + (-t220 / (sigma * 2.0) + t202
                - t022 / (sigma * sigma)
                - 2.0 * t004 / sigma)
                * t022
                / b020;

        ret.hj = phi600 - lamda;
        ret.hjk = phi420 + 6.0 * lamda - 3.0 * mu;
        ret.hkj = phi240 - 5.0 * lamda + 10.0 * mu;
        ret.hk = phi060 - 7.0 * mu;
        ret.h1 = phi402 - nu;
        ret.h2 = phi204 + lamda / 2.0;
        ret.h3 = phi006 + nu;

        // TODO spherical top case

        ret
    }
}

fn make_phi(
    maxcor: usize,
    ic: [usize; 3],
    scc: tensor::Tensor3<f64>,
) -> tensor::Tensor3<f64> {
    let mut phi = Tensor3::zeros(maxcor, maxcor, maxcor);
    for ixyz in 0..maxcor {
        for jxyz in 0..maxcor {
            for kxyz in 0..maxcor {
                phi[(ic[ixyz], ic[jxyz], ic[kxyz])] = scc[(ixyz, jxyz, kxyz)];
            }
        }
    }
    phi
}

#[cfg(test)]
impl Display for Sextic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", serde_json::to_string_pretty(self).unwrap())
    }
}

#[cfg(not(test))]
impl Display for Sextic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Watson A Reduction")?;
        writeln!(f, "Phi  J: {:20.12e}", self.phij)?;
        writeln!(f, "Phi  K: {:20.12e}", self.phik)?;
        writeln!(f, "Phi JK: {:20.12e}", self.phijk)?;
        writeln!(f, "Phi KJ: {:20.12e}", self.phikj)?;
        writeln!(f, "phi  J: {:20.12e}", self.sphij)?;
        writeln!(f, "phi JK: {:20.12e}", self.sphijk)?;
        writeln!(f, "phi  k: {:20.12e}", self.sphik)?;

        writeln!(f)?;

        writeln!(f, "Watson S Reduction")?;
        writeln!(f, "H  J: {:20.12e}", self.hj)?;
        writeln!(f, "H  K: {:20.12e}", self.hk)?;
        writeln!(f, "H JK: {:20.12e}", self.hjk)?;
        writeln!(f, "H KJ: {:20.12e}", self.hkj)?;
        writeln!(f, "h  1: {:20.12e}", self.h1)?;
        writeln!(f, "h  2: {:20.12e}", self.h2)?;
        writeln!(f, "h  3: {:20.12e}", self.h3)?;

        writeln!(f)?;

        writeln!(f, "Linear Molecule")?;
        writeln!(f, "He  : {:20.12e}", self.he)?;

        Ok(())
    }
}
