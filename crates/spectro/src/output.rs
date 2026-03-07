use std::fmt::Display;

use serde::{Deserialize, Serialize};
use symm::{Irrep, Molecule};

use crate::{
    F3qcm, F4qcm,
    quartic::Quartic,
    resonance::{Coriolis, Fermi1, Fermi2, Restst},
    rot::Rot,
    sextic::Sextic,
};

/// contains all of the output data from running Spectro
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Output {
    /// harmonic frequencies
    pub harms: Vec<f64>,

    /// partially resonance-corrected anharmonic frequencies
    pub funds: Vec<f64>,

    /// fully resonance-corrected anharmonic frequencies
    pub corrs: Vec<f64>,

    /// vibrationally averaged rotational constants
    pub rots: Vec<Rot>,

    /// equilibrium rotational constants
    pub rot_equil: Vec<f64>,

    pub irreps: Vec<Irrep>,

    /// quartic distortion coefficients
    pub quartic: Quartic,

    /// sextic distortion coefficients
    pub sextic: Sextic,

    /// zero-point vibrational energy
    pub zpt: f64,

    pub geom: Molecule,

    pub lxm: Vec<Vec<f64>>,

    #[serde(default)]
    pub lx: Vec<Vec<f64>>,

    pub linear: bool,

    pub resonances: Restst,

    #[serde(default)]
    pub f3qcm: Option<F3qcm>,

    #[serde(default)]
    pub f4qcm: Option<F4qcm>,
}

impl Display for Output {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let Output {
            harms,
            funds,
            corrs,
            rots,
            rot_equil,
            irreps,
            quartic,
            sextic,
            zpt,
            geom,
            lxm: _,
            lx: _,
            linear,
            resonances:
                Restst {
                    coriolis,
                    fermi1,
                    fermi2,
                    darling: _,
                    states: _,
                    modes: _,
                    ifunda: _,
                    iovrtn: _,
                    icombn: _,
                },
            f3qcm,
            f4qcm,
        } = self;
        writeln!(f, "Geometry: {geom:.8}")?;
        writeln!(
            f,
            "Vibrational Frequencies (cm-1):\n{:>5}{:>8}{:>8}{:>8}{:>8}",
            "Mode", "Symm", "Harm", "Fund", "Corr"
        )?;
        for i in 0..harms.len() {
            writeln!(
                f,
                "{:5}{:>8}{:8.1}{:8.1}{:8.1}",
                i + 1,
                irreps[i],
                harms[i],
                funds[i],
                corrs[i]
            )?;
        }

        writeln!(f, "\nZPT = {zpt:8.1}")?;

        writeln!(f, "\nRotational Constants (cm-1):")?;
        if !rots.is_empty() {
            let width = 5 * rots[0].state.len();
            if !linear {
                writeln!(
                    f,
                    "{:^width$}{:^12}{:^12}{:^12}",
                    "State", "A", "B", "C"
                )?;
                for rot in rots {
                    writeln!(f, "{rot}")?;
                }
            } else {
                writeln!(f, "{:^width$}{:^12}", "State", "B")?;
                for rot in rots {
                    writeln!(f, "{}{:12.7}", rot.state, rot.b + rot_equil[1])?;
                }
            }
        }

        writeln!(f, "\nQuartic Distortion Constants (cm-1):\n{quartic}")?;

        writeln!(f, "Sextic Distortion Constants (cm-1):\n{sextic}")?;

        writeln!(f, "\nCoriolis Resonances:")?;
        for Coriolis { i, j, axis } in coriolis {
            writeln!(f, "{i:5}{j:5}{axis:>5}")?;
        }

        writeln!(f, "\nType 1 Fermi Resonances:")?;
        for Fermi1 { i, j } in fermi1 {
            writeln!(f, "{i:5}{j:5}")?;
        }

        writeln!(f, "\nType 2 Fermi Resonances:")?;
        for Fermi2 { i, j, k } in fermi2 {
            writeln!(f, "{i:5}{j:5}{k:5}")?;
        }

        let nvib = harms.len();

        // Threshold for considering force constants to be zero.
        const FC_THRESH: f64 = 1e-8;

        writeln!(f, "\nHarmonic force constants (cm-1)")?;
        for (i, value) in harms.iter().enumerate() {
            if value.abs() < FC_THRESH {
                continue;
            }
            let mode = i + 1;
            writeln!(f, "{mode:5}{mode:5} {value:15.8}")?;
        }

        if let Some(f3qcm) = f3qcm {
            writeln!(f, "\nNon-zero cubic force constants (cm-1)")?;
            for i in 0..nvib {
                for j in 0..=i {
                    for k in 0..=j {
                        let value = f3qcm[(i, j, k)];
                        if value.abs() < FC_THRESH {
                            continue;
                        }
                        let (i, j, k) = (i + 1, j + 1, k + 1);
                        writeln!(f, "{i:5}{j:5}{k:5} {value:15.8}")?;
                    }
                }
            }
        }

        if let Some(f4qcm) = f4qcm {
            writeln!(f, "\nNon-zero quartic force constants (cm-1)")?;
            for i in 0..nvib {
                for j in 0..=i {
                    for k in 0..=j {
                        for l in 0..=k {
                            let value = f4qcm[(i, j, k, l)];
                            if value.abs() < FC_THRESH {
                                continue;
                            }
                            let (i, j, k, l) = (i + 1, j + 1, k + 1, l + 1);
                            writeln!(f, "{i:5}{j:5}{k:5}{l:5}{value:15.8}")?;
                        }
                    }
                }
            }
        }

        Ok(())
    }
}

impl Output {
    pub fn to_json(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string(&self)
    }

    pub fn to_json_pretty(&self) -> Result<String, serde_json::Error> {
        serde_json::to_string_pretty(&self)
    }
}
