use crate::{sextic::Sextic, *};

use super::*;

#[derive(Clone)]
struct Test {
    infile: String,
    fort15: String,
    fort30: String,
    want: Sextic,
}

impl Test {
    fn new(dir: &'static str) -> Self {
        let start = Path::new("testfiles");
        let w = File::open(start.join(dir).join("phis.json")).unwrap();
        let want: Sextic = serde_json::from_reader(w).unwrap();
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
            want,
        }
    }
}

impl Sub<Sextic> for Sextic {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self {
            phij: self.phij - rhs.phij,
            phijk: self.phijk - rhs.phijk,
            phikj: self.phikj - rhs.phikj,
            phik: self.phik - rhs.phik,
            sphij: self.sphij - rhs.sphij,
            sphijk: self.sphijk - rhs.sphijk,
            sphik: self.sphik - rhs.sphik,
            hj: self.hj - rhs.hj,
            hjk: self.hjk - rhs.hjk,
            hkj: self.hkj - rhs.hkj,
            hk: self.hk - rhs.hk,
            h1: self.h1 - rhs.h1,
            h2: self.h2 - rhs.h2,
            h3: self.h3 - rhs.h3,
            he: self.he - rhs.he,
        }
    }
}

macro_rules! check {
    ($got:expr, $test:ident) => {
        if abs_diff_ne!($got, $test.want, epsilon = 2e-10) {
            println!("\n{}", $test.infile);
            println!("got\n{}", $got);
            println!("want\n{}", $test.want);
            println!("diff\n{}", $got - $test.want);
            panic!("failed, {}", $test.infile);
        }
    };
}

#[test]
fn asym() {
    let tests = [
        Test::new("h2o"),
        Test::new("h2co"),
        Test::new("c3h2"),
        Test::new("c3hf"),
        Test::new("c3hcn"),
    ];
    for test in Vec::from(&tests[..]) {
        let s = Spectro::load(&test.infile);
        let fc2 = load_fc2(test.fort15, s.n3n);
        let fc2 = s.rot2nd(&fc2);
        let fc2 = FACT2 * fc2;
        let w = s.geom.weights();
        let sqm: Vec<_> = w.iter().map(|w| 1.0 / w.sqrt()).collect();
        let fxm = s.form_sec(fc2, &sqm);
        let (harms, lxm) = utils::linalg::symm_eigen_decomp(fxm, true);
        let freq = to_wavenumbers(&harms);
        let lx = s.make_lx(&sqm, &lxm);
        let (zmat, wila) = s.zeta(&lxm, &w);
        let f3x = load_fc3(test.fort30, s.n3n);
        let f3x = s.rot3rd(f3x);
        let f3qcm = force3(s.n3n, f3x, &lx, s.nvib, &freq);
        let got = Sextic::new(&s, &wila, &zmat, &freq, &f3qcm);

        check!(got, test);
    }
}

#[test]
fn sym() {
    let tests = [
        Test::new("nh3"),
        Test::new("ph3"),
        Test::new("bipy"),
        // linear
        Test::new("c2h-"),
        // diatomic
        Test::new("alh"),
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
        if s.rotor.is_sym_top() {
            s.bdegnl(&freq, &mut lxm, &w, &mut lx);
        }
        let f3x = load_fc3(test.fort30, s.n3n);
        let f3x = s.rot3rd(f3x);
        let f3qcm = force3(s.n3n, f3x, &lx, s.nvib, &freq);
        let (zmat, wila) = s.zeta(&lxm, &w);

        let got = Sextic::new(&s, &wila, &zmat, &freq, &f3qcm);

        check!(got, test);
    }
}
