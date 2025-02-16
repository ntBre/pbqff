use std::{
    cmp::{max, min},
    collections::HashSet,
    fmt::Display,
};

use serde::{Deserialize, Serialize};

use crate::{
    f3qcm::F3qcm,
    mode::Mode,
    rotor::{coriolis, darling, fermi1, fermi2},
    state::State,
    utils::{close, ioff},
    Dvec, Spectro, Tensor3,
};

#[derive(
    Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize,
)]
pub enum Axis {
    #[default]
    A = 0,
    B = 1,
    C = 2,
}

impl Display for Axis {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.pad(match self {
            Axis::A => "A",
            Axis::B => "B",
            Axis::C => "C",
        })
    }
}

/// coriolis resonance wᵢ = wⱼ
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Coriolis {
    pub i: usize,
    pub j: usize,
    pub axis: Axis,
}

impl From<usize> for Axis {
    fn from(i: usize) -> Self {
        match i {
            0 => Axis::A,
            1 => Axis::B,
            2 => Axis::C,
            _ => panic!("unmatched axis"),
        }
    }
}

impl From<i32> for Axis {
    fn from(i: i32) -> Self {
        match i {
            0 => Axis::A,
            1 => Axis::B,
            2 => Axis::C,
            _ => panic!("unmatched axis"),
        }
    }
}

impl Coriolis {
    pub fn new(i: usize, j: usize, axis: impl Into<Axis>) -> Self {
        Self {
            i,
            j,
            axis: axis.into(),
        }
    }
}

/// type 1 Fermi resonance 2wᵢ = wⱼ
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Fermi1 {
    pub i: usize,
    pub j: usize,
}

impl Fermi1 {
    pub fn new(i: usize, j: usize) -> Self {
        Self { i, j }
    }
}

/// type 2 Fermi resonance wₖ = wⱼ + wᵢ
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Fermi2 {
    pub i: usize,
    pub j: usize,
    pub k: usize,
}

impl Fermi2 {
    pub fn new(i: usize, j: usize, k: usize) -> Self {
        Self { i, j, k }
    }
}

impl Display for Fermi2 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:5}{:5}{:5}", self.i, self.j, self.k)
    }
}

/// Darling-Dennison resonance 2wᵢ = 2wⱼ
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct Darling {
    pub i: usize,
    pub j: usize,
}

impl Darling {
    pub fn new(i: usize, j: usize) -> Self {
        Self { i, j }
    }
}

#[derive(Clone, Debug, Default, Eq, PartialEq, Serialize, Deserialize)]
pub struct Restst {
    pub coriolis: Vec<Coriolis>,
    pub fermi1: Vec<Fermi1>,
    pub fermi2: Vec<Fermi2>,
    pub darling: Vec<Darling>,
    pub states: Vec<State>,
    pub modes: Vec<Mode>,
    pub ifunda: Vec<usize>,
    pub iovrtn: Vec<usize>,
    pub icombn: Vec<usize>,
}

impl Restst {
    // should return all of the resonances, as well as the states (I think)
    pub fn new(
        spectro: &Spectro,
        zmat: &Tensor3,
        f3qcm: &F3qcm,
        freq: &Dvec,
    ) -> Self {
        let modes = deg_modes(spectro, freq);

        let r = &spectro.rotor;
        let coriolis = coriolis(r, &modes, freq, zmat);
        let fermi1 = fermi1(r, &modes, freq, f3qcm);
        let fermi2 = fermi2(r, &modes, freq, f3qcm);
        let darling = darling(r, &modes, freq);

        let (n1dm, n2dm, n3dm) = Mode::count(&modes);
        let (i1mode, i2mode, i3mode) = Mode::partition(&modes);
        if n3dm > 0 {
            eprintln!("triply degenerate modes are untested");
        }

        let mut states = Vec::new();
        use crate::state::State::*;

        let mut ifunda = vec![0; spectro.nvib];
        let mut iovrtn = vec![0; spectro.nvib];
        let mut icombn = vec![0; spectro.i2vib];
        let mut ist = 1;
        // ground state, all zeros
        states.push(I1st(vec![0; spectro.nvib]));

        // fundamentals, single excitations
        for ii in 0..n1dm {
            let mut tmp = vec![0; spectro.nvib];
            tmp[ii] = 1;
            states.push(I1st(tmp));
            let i = i1mode[ii];
            ifunda[i] = ist;
            ist += 1;
        }
        for ii in 0..n2dm {
            let mut tmp = vec![(0, 0); spectro.nvib];
            tmp[ii] = (1, 1);
            states.push(I2st(tmp));
            let (i1, i2) = i2mode[ii];
            ifunda[i1] = ist;
            ifunda[i2] = ist;
            ist += 1
        }
        for ii in 0..n3dm {
            let mut tmp = vec![0; spectro.nvib];
            tmp[ii] = 1;
            states.push(I3st(tmp));
            let (i1, i2, i3) = i3mode[ii];
            ifunda[i1] = ist;
            ifunda[i2] = ist;
            ifunda[i3] = ist;
            ist += 1
        }
        // NOTE: triply-degenerate modes are only handled for fundamentals. also
        // their resonances are not handled at all

        // overtones, double excitations in one mode
        for ii in 0..n1dm {
            let mut tmp = vec![0; spectro.nvib];
            tmp[ii] = 2;
            states.push(I1st(tmp));
            let i = i1mode[ii];
            iovrtn[i] = ist;
            ist += 1;
        }
        for ii in 0..n2dm {
            let mut tmp = vec![(0, 0); spectro.nvib];
            tmp[ii] = (2, 0);
            states.push(I2st(tmp));
            let (i1, i2) = i2mode[ii];
            iovrtn[i1] = ist;
            iovrtn[i2] = ist;
            ist += 1;
        }

        // combination bands
        for ii in 0..n1dm {
            let i = i1mode[ii];
            // nondeg - nondeg combination
            for jj in ii + 1..n1dm {
                let mut tmp = vec![0; spectro.nvib];
                tmp[ii] = 1;
                tmp[jj] = 1;
                states.push(I1st(tmp));
                let j = i1mode[jj];
                let ij = ioff(max(i, j) + 1) + min(i, j);
                icombn[ij] = ist;
                ist += 1;
            }

            // nondeg - deg combination
            for jj in 0..n2dm {
                // I guess you push this to states as well. might have to change
                // how I index these because they set states(ii, ist) to this,
                // which I guess syncs it with ist in i2sts. we'll see how it's
                // accessed. I'm not that interested in combination bands anyway
                let mut i1st = vec![0; spectro.nvib];
                i1st[ii] = 1;
                let mut i2st = vec![(0, 0); spectro.nvib];
                i2st[jj] = (1, 1);
                states.push(I12st { i1st, i2st });
                let (j1, j2) = i2mode[jj];
                let ij1 = ioff(max(i, j1) + 1) + min(i, j1);
                let ij2 = ioff(max(i, j2) + 1) + min(i, j2);
                icombn[ij1] = ist;
                icombn[ij2] = ist;
                ist += 1;
            }
        }

        // deg-deg combination
        for ii in 0..n2dm {
            let (i1, i2) = i2mode[ii];
            for jj in ii + 1..n2dm {
                let (j1, j2) = i2mode[jj];
                let mut tmp = vec![(0, 0); spectro.nvib];
                tmp[ii] = (1, 1);
                tmp[jj] = (1, 1);
                states.push(I2st(tmp));
                let i1j1 = ioff(max(i1, j1) + 1) + min(i1, j1);
                let i1j2 = ioff(max(i1, j2) + 1) + min(i1, j2);
                let i2j1 = ioff(max(i2, j1) + 1) + min(i2, j1);
                let i2j2 = ioff(max(i2, j2) + 1) + min(i2, j2);
                icombn[i1j1] = ist;
                icombn[i1j2] = ist;
                icombn[i2j1] = ist;
                icombn[i2j2] = ist;
                ist += 1;
            }
        }

        Self {
            coriolis,
            fermi1,
            fermi2,
            darling,
            states,
            modes,
            ifunda,
            iovrtn,
            icombn,
        }
    }
}

pub(crate) fn deg_modes(spectro: &Spectro, freq: &Dvec) -> Vec<Mode> {
    let mut modes = Vec::new();
    // probably I could get rid of this if, but I guess it protects against
    // accidental degmodes in asymmetric tops. is an accident still
    // degenerate?
    use Mode::*;
    if spectro.rotor.is_sym_top() {
        const DEG_TOL: f64 = 0.1;
        let mut triples = HashSet::new();
        for i in 0..spectro.nvib {
            let fi = freq[i];
            for j in i + 1..spectro.nvib {
                let fj = freq[j];
                for k in j + 1..spectro.nvib {
                    let fk = freq[k];
                    if close(fi, fj, DEG_TOL)
                        && close(fi, fk, DEG_TOL)
                        && [fi, fj, fk].iter().all(|f| f.abs() > DEG_TOL)
                    {
                        triples.insert(i);
                        triples.insert(j);
                        triples.insert(k);
                        modes.push(I3(i, j, k));
                    }
                }
            }
        }
        // run these in two passes to avoid counting something as
        // triply-degenerate AND doubly-degenerate
        let mut doubles = HashSet::new();
        for i in 0..spectro.nvib {
            let fi = freq[i];
            for j in i + 1..spectro.nvib {
                let fj = freq[j];
                if close(fi, fj, DEG_TOL)
                    && [fi, fj].iter().all(|f| f.abs() > DEG_TOL)
                    && !triples.contains(&i)
                    && !triples.contains(&j)
                {
                    doubles.insert(i);
                    doubles.insert(j);
                    modes.push(I2(i, j));
                }
            }
        }
        for i in 0..spectro.nvib {
            if !triples.contains(&i) && !doubles.contains(&i) {
                modes.push(I1(i));
            }
        }
    } else {
        for i in 0..spectro.nvib {
            modes.push(I1(i));
        }
    };
    modes
}
