use super::{check_mat, load_dmat};
use crate::{
    Dmat, Spectro,
    consts::FACT2,
    load_fc2, load_fc3,
    resonance::Restst,
    tests::force3,
    utils::{linalg::symm_eigen_decomp, to_wavenumbers},
};
use crate::{Tensor3, f3qcm::F3qcm, utils::load_vec};
use nalgebra::dvector;
use std::collections::HashMap;
use std::path::Path;

#[test]
#[allow(clippy::excessive_precision)]
fn make_alpha() {
    let s = Spectro {
        rotcon: vec![
            4.5034193376213194,
            4.5034192814538905,
            3.9601029298042021,
        ],
        primat: vec![
            3.7432957395278987,
            3.7432957862149348,
            4.2568667326682146,
        ],
        nvib: 6,
        ..Spectro::default()
    };
    let n1dm = 2;
    let n2dm = 2;
    let i1mode = vec![2, 5];
    let i2mode = vec![(0, 1), (3, 4)];
    let ia = 2;
    let ib = 1;
    let ic = 0;
    let iaia = 5;
    let iaib = 4;
    let ibib = 2;
    let ibic = 1;
    let freq = dvector![
        2437.0024382429601,
        2437.0020300169222,
        2428.1546959531433,
        1145.7614713948005,
        1145.7567008646774,
        1015.4828482214278,
        0.046516934814435321,
        0.025278773879589642,
        0.015292949295790762,
        0.0077266799336841666,
        0.011089653435072953,
        0.014380890514910207
    ];
    let zmat = Tensor3::load("testfiles/ph3/zmat");
    let wila = load_dmat("testfiles/ph3/wila", 6, 6);
    let f3qcm = F3qcm::new(load_vec("testfiles/ph3/f3qcm"));
    let icorol = HashMap::from([((3, 5), 1), ((5, 3), 1)]);
    let got = s.make_alpha(
        n1dm, &i1mode, ia, &freq, &wila, iaia, icorol, &zmat, &f3qcm, n2dm,
        &i2mode, iaib, ib, ibib, ibic, ic,
    );
    let want = load_dmat("testfiles/ph3/alpha", 6, 3);

    check_mat(&got, &want, 1e-10, "ph3");
}

#[derive(Clone)]
struct Test {
    infile: String,
    fort15: String,
    fort30: String,
    want: Dmat,
}

impl Test {
    fn new(dir: &'static str, states: usize) -> Self {
        let start = Path::new("testfiles");
        Self {
            infile: String::from(
                start.join(dir).join("spectro.in").to_str().unwrap(),
            ),
            fort15: String::from(
                start.join(dir).join("fort.15").to_str().unwrap(),
            ),
            fort30: String::from(
                start.join(dir).join("fort.30").to_str().unwrap(),
            ),
            want: load_dmat(start.join(dir).join("rotnst"), states, 3),
        }
    }
}

/// take the values from alphas.f:792
#[test]
fn rotnst() {
    let tests = [
        //
        Test::new("bipy", 66),
    ];
    for test in Vec::from(&tests[..]) {
        let s = Spectro::load(&test.infile);
        let fc2 = load_fc2(test.fort15, s.n3n);
        let fc2 = s.rot2nd(&fc2);
        let fc2 = FACT2 * fc2;
        let w = s.geom.weights();
        let sqm: Vec<_> = w.iter().map(|w| 1.0 / w.sqrt()).collect();
        let fxm = s.form_sec(fc2, &sqm);
        let (harms, mut lxm) = symm_eigen_decomp(fxm, true);
        let freq = to_wavenumbers(&harms);
        let mut lx = s.make_lx(&sqm, &lxm);
        s.bdegnl(&freq, &mut lxm, &w, &mut lx);
        let f3x = load_fc3(test.fort30, s.n3n);
        let f3x = s.rot3rd(f3x);
        let f3qcm = force3(s.n3n, f3x, &lx, s.nvib, &freq);

        let (zmat, wila) = s.zeta(&lxm, &w);

        let r = Restst::new(&s, &zmat, &f3qcm, &freq);

        let got = s.alphas(
            &freq,
            &wila,
            &zmat,
            &f3qcm,
            &r.modes,
            &r.states,
            &r.coriolis,
        );

        check_mat!(&got, &test.want, 5e-13, test.infile);
    }
}
