use crate::{consts::FACT2, utils::linalg::symm_eigen_decomp, *};

use super::check_vec;

#[derive(Clone)]
struct Test {
    infile: String,
    fort15: String,
    fort30: String,
    want: Vec<f64>,
    eps: f64,
}

impl Test {
    fn new(dir: &'static str, eps: f64) -> Self {
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
            want: load_vec(start.join(dir).join("f3qcm")),
            eps,
        }
    }
}

fn inner(tests: &[Test]) {
    for test in tests {
        let s = Spectro::load(&test.infile);
        let fc2 = load_fc2(&test.fort15, s.n3n);
        let fc2 = s.rot2nd(&fc2);
        let fc2 = FACT2 * fc2;
        let w = s.geom.weights();
        let sqm: Vec<_> = w.iter().map(|w| 1.0 / w.sqrt()).collect();
        let fxm = s.form_sec(fc2, &sqm);
        let (harms, lxm) = symm_eigen_decomp(fxm, true);
        let freq = to_wavenumbers(&harms);
        let lx = s.make_lx(&sqm, &lxm);
        let f3x = load_fc3(&test.fort30, s.n3n);
        let f3x = s.rot3rd(f3x);
        let got = force3(s.n3n, f3x, &lx, s.nvib, &freq);
        let got = Dvec::from(got).abs();
        let want = Dvec::from(test.want.clone()).abs();
        check_vec!(got, want, test.eps, &test.infile);
    }
}

#[test]
fn asym() {
    let tests = [
        Test::new("h2o", 4e-6),
        Test::new("h2co", 3e-6),
        Test::new("c3h2", 3e-5),
        Test::new("c3hf", 5e-6),
        Test::new("c3hcn", 5e-6),
    ];
    inner(&tests);
}

#[test]
fn sym() {
    let tests = [
        Test::new("nh3", 5e-6),
        Test::new("ph3", 3e-6),
        Test::new("bipy", 3e-6),
        Test::new("c2h-", 3e-6),
        Test::new("hmgnc", 3e-6),
    ];
    for test in tests {
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
        s.bdegnl(&freq, &mut lxm, &w, &mut lx);
        let f3x = load_fc3(&test.fort30, s.n3n);
        let f3x = s.rot3rd(f3x);
        let got = force3(s.n3n, f3x, &lx, s.nvib, &freq);
        let got = Dvec::from(got).abs();
        let want = Dvec::from(test.want.clone()).abs();
        check_vec!(got, want, test.eps, &test.infile);
    }
}
