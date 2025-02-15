use std::path::{Path, PathBuf};

use approx::assert_abs_diff_eq;

use crate::{
    consts::FACT2,
    load_fc2, load_fc3, load_fc4,
    resonance::Restst,
    utils::{linalg::symm_eigen_decomp, to_wavenumbers},
    Dmat, Spectro,
};

use super::{force3, force4};

#[derive(Clone)]
pub(crate) struct Test {
    infile: PathBuf,
    fort15: PathBuf,
    fort30: PathBuf,
    fort40: PathBuf,
    want: Dmat,
    eps: f64,
}

impl Test {
    fn new(dir: &'static str, eps: f64) -> Self {
        let start = Path::new("testfiles");
        Self {
            infile: start.join(dir).join("spectro.in"),
            fort15: start.join(dir).join("fort.15"),
            fort30: start.join(dir).join("fort.30"),
            fort40: start.join(dir).join("fort.40"),
            want: load_fc2(start.join(dir).join("polyad"), 7),
            eps,
        }
    }
}

#[test]
fn c2h4_polyad() {
    let test = Test::new("c2h4", 1e-3);
    let s = Spectro::load(&test.infile);
    let fc2 = load_fc2(&test.fort15, s.n3n);
    let fc2 = s.rot2nd(&fc2);
    let fc2 = FACT2 * fc2;
    let w = s.geom.weights();
    let sqm: Vec<_> = w.iter().map(|w| 1.0 / w.sqrt()).collect();
    let fxm = s.form_sec(fc2, &sqm);
    let (harms, mut lxm) = symm_eigen_decomp(fxm, true);
    let freq = to_wavenumbers(&harms);
    let mut lx = s.make_lx(&sqm, &lxm);
    if s.rotor.is_sym_top() {
        s.bdegnl(&freq, &mut lxm, &w, &mut lx);
    }
    let f3x = load_fc3(&test.fort30, s.n3n);
    let f3x = s.rot3rd(f3x);
    let f3qcm = force3(s.n3n, f3x, &lx, s.nvib, &freq);
    let f4x = load_fc4(&test.fort40, s.n3n);
    let f4x = s.rot4th(f4x);
    let f4qcm = force4(s.n3n, &f4x, &lx, s.nvib, &freq);
    let (zmat, wila) = s.zeta(&lxm, &w);
    let restst = Restst::new(&s, &zmat, &f3qcm, &freq);
    let Restst {
        fermi1,
        fermi2,
        modes,
        states,
        ..
    } = &restst;
    let (xcnst, _, e0) = if s.rotor.is_sym_top() {
        let (x, g, e) =
            s.xcals(&f4qcm, &freq, &f3qcm, &zmat, fermi1, fermi2, modes, &wila);
        (x, Some(g), e)
    } else {
        let (x, e) =
            s.xcalc(&f4qcm, &freq, &f3qcm, &zmat, modes, fermi1, fermi2);
        (x, None, e)
    };

    let mut eng = vec![0.0; states.len()];
    let got = crate::polyads::resona(
        &zmat, &f3qcm, &f4qcm, e0, modes, &freq, &s.rotcon, &xcnst, fermi1,
        fermi2, &mut eng,
    );

    assert_abs_diff_eq!(
        got.unwrap().abs(),
        test.want.abs().transpose(),
        epsilon = test.eps
    );
}
