use std::path::PathBuf;

use approx::{abs_diff_ne, assert_abs_diff_eq};
use nalgebra::dmatrix;

use crate::{
    consts::FACT2, resonance::Restst, utils::linalg::symm_eigen_decomp, *,
};

use super::load_dmat;

#[derive(Clone)]
struct Test {
    infile: PathBuf,
    fort15: PathBuf,
    fort30: PathBuf,
    want: Dmat,
}

impl Test {
    /// take want as a str so it can be used for both tests herein
    fn new(dir: &'static str, want: &str, size: usize) -> Self {
        let start = Path::new("testfiles");
        Self {
            infile: start.join(dir).join("spectro.in"),
            fort15: start.join(dir).join("fort.15"),
            fort30: start.join(dir).join("fort.30"),
            want: load_dmat(start.join(dir).join(want), size, 3),
        }
    }
}

#[test]
fn alpha() {
    let tests = [
        Test::new("h2o", "alpha", 3),
        Test::new("h2co", "alpha", 6),
        Test::new("c3h2", "alpha", 9),
        Test::new("c3hf", "alpha", 9),
        Test::new("c3hcn", "alpha", 12),
        Test::new("c3hcn010", "alpha", 12),
    ];
    for test in Vec::from(&tests[..]) {
        let s = Spectro::load(test.infile);
        let fc2 = load_fc2(test.fort15, s.n3n);
        let fc2 = s.rot2nd(&fc2);
        let fc2 = FACT2 * fc2;
        let w = s.geom.weights();
        let sqm: Vec<_> = w.iter().map(|w| 1.0 / w.sqrt()).collect();
        let fxm = s.form_sec(fc2, &sqm);
        let (harms, lxm) = symm_eigen_decomp(fxm, true);
        let freq = to_wavenumbers(&harms);
        let lx = s.make_lx(&sqm, &lxm);
        let (zmat, wila) = s.zeta(&lxm, &w);
        let f3x = load_fc3(test.fort30, s.n3n);
        let f3x = s.rot3rd(f3x);
        let f3qcm = force3(s.n3n, f3x, &lx, s.nvib, &freq);
        let r = Restst::new(&s, &zmat, &f3qcm, &freq);
        let got = s.alpha(&freq, &wila, &zmat, &f3qcm, &r.coriolis);
        if abs_diff_ne!(got, test.want, epsilon = 3e-6) {
            println!("got\n{got:.8}");
            println!("want\n{:.8}", test.want);
            println!("diff\n{:.8}", got.clone() - test.want.clone());
            println!(
                "max diff = {:.2e}",
                (got.clone() - test.want.clone()).abs().max()
            );
        }
    }
}

#[test]
fn test_alphaa() {
    let s = Spectro::load("testfiles/h2o/spectro.in");
    let fc2 = load_fc2("testfiles/fort.15", s.n3n);
    let fc2 = s.rot2nd(&fc2);
    let fc2 = FACT2 * fc2;
    let w = s.geom.weights();
    let sqm: Vec<_> = w.iter().map(|w| 1.0 / w.sqrt()).collect();
    let fxm = s.form_sec(fc2, &sqm);
    let (harms, lxm) = symm_eigen_decomp(fxm, true);
    let freq = to_wavenumbers(&harms);
    let lx = s.make_lx(&sqm, &lxm);
    let (zmat, wila) = s.zeta(&lxm, &w);
    let f3x = load_fc3("testfiles/fort.30", s.n3n);
    let f3x = s.rot3rd(f3x);
    let f3qcm = force3(s.n3n, f3x, &lx, s.nvib, &freq);
    let Restst {
        coriolis,
        fermi1: _,
        fermi2: _,
        darling: _,
        states,
        modes,
        ifunda: _,
        iovrtn: _,
        icombn: _,
    } = Restst::new(&s, &zmat, &f3qcm, &freq);
    let got = s.alphaa(&freq, &wila, &zmat, &f3qcm, &modes, &states, &coriolis);
    let want = dmatrix![
    27.657417987118755, 14.498766626639174, 9.267_303_844_958_324;
     26.500953159400968, 14.400078799009306, 9.124_329_013_234_416;
     26.968_408_648_401_32, 14.279451668708212, 9.090_265_821_548_55;
     30.255750931866725, 14.660626578959441, 9.120_911_467_809_906;
    ];
    assert_abs_diff_eq!(got, want, epsilon = 2e-5);
}
