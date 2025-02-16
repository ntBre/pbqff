use crate::{consts::FACT2, utils::linalg::symm_eigen_decomp, *};

use super::check_vec;

#[derive(Clone)]
struct Test {
    infile: String,
    fort15: String,
    fort40: String,
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
            fort40: String::from(
                start.join(dir).join("fort.40").to_str().unwrap(),
            ),
            want: load_vec(start.join(dir).join("f4qcm")),
            eps,
        }
    }
}

#[test]
fn asym() {
    let tests = [
        Test::new("h2o", 2.2e-6),
        Test::new("h2co", 1.7e-6),
        Test::new("c3h2", 7e-6),
        Test::new("c3hf", 3.1e-6),
        Test::new("c3hcn", 3.2e-6),
    ];
    for test in tests {
        let s = Spectro::load(&test.infile);
        let fc2 = load_fc2(test.fort15, s.n3n);
        let fc2 = s.rot2nd(&fc2);
        let fc2 = FACT2 * fc2;
        let w = s.geom.weights();
        let sqm: Vec<_> = w.iter().map(|w| 1.0 / w.sqrt()).collect();
        let fxm = s.form_sec(fc2, &sqm);
        let (harms, lxm) = symm_eigen_decomp(fxm, true);
        let freq = to_wavenumbers(&harms);
        let lx = s.make_lx(&sqm, &lxm);
        let f4x = load_fc4(test.fort40, s.n3n);
        let f4x = s.rot4th(f4x);
        let got = force4(s.n3n, &f4x, &lx, s.nvib, &freq);
        let got = Dvec::from(got).abs();
        let want = Dvec::from(test.want).abs();
        check_vec!(got, want, test.eps, &test.infile);
    }
}

#[test]
fn sym() {
    let tests = [
        Test::new("nh3", 2.2e-6),
        Test::new("ph3", 2.2e-6),
        Test::new("bipy", 2.2e-6),
        Test::new("hmgnc", 2.2e-6),
    ];
    for test in tests {
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
        let f4x = load_fc4(test.fort40, s.n3n);
        let f4x = s.rot4th(f4x);
        let got = force4(s.n3n, &f4x, &lx, s.nvib, &freq);
        let got = Dvec::from(got).abs();
        let want = Dvec::from(test.want).abs();
        check_vec!(got, want, test.eps, &test.infile);
    }
}
