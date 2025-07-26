use std::{fs::read_to_string, path::Path};

use approx::{abs_diff_ne, assert_abs_diff_eq};

use crate::{
    Spectro,
    consts::FACT2,
    load_fc2,
    quartic::Quartic,
    utils::{linalg::symm_eigen_decomp, to_wavenumbers},
};

#[derive(Clone)]
struct Test {
    infile: String,
    fort15: String,
    want: Quartic,
}

impl Test {
    fn new(dir: &'static str) -> Self {
        let start = Path::new("testfiles");
        let data =
            read_to_string(start.join(dir).join("quartic.json")).unwrap();
        let want: Quartic = serde_json::from_str(&data).unwrap();
        Self {
            infile: String::from(
                start.join(dir).join("spectro.in").to_str().unwrap(),
            ),
            fort15: String::from(
                start.join(dir).join("fort.15").to_str().unwrap(),
            ),
            want,
        }
    }
}

/// generate these tests with the pq function in the old-spectro/gdb/help.py
/// directory with a break point at qcent.f:443
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
        let sqm: Vec<_> = w.iter().map(|w: &f64| 1.0 / w.sqrt()).collect();
        let fxm = s.form_sec(fc2, &sqm);
        let (harms, lxm) = symm_eigen_decomp(fxm, true);
        let freq = to_wavenumbers(&harms);
        let (_zmat, wila) = s.zeta(&lxm, &w);
        let got = Quartic::new(&s, &freq, &wila);

        // accept this size of epsilon because this is about how good the
        // rotational constant agreement is and b[xyz][as] are the largest
        // differences
        if abs_diff_ne!(got, test.want, epsilon = 2e-5) {
            println!("got\n{got}");
            println!("want\n{}", test.want);
            println!("diff\n{}", got - test.want.clone());
            panic!("{} failed", test.infile);
        }
    }
}

/// this only checks bxs, bys, and bzs, the values that are used in rots
#[test]
fn sym() {
    let tests = [
        Test::new("nh3"),
        Test::new("ph3"),
        Test::new("bipy"),
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
        s.bdegnl(&freq, &mut lxm, &w, &mut lx);
        let (_zmat, wila) = s.zeta(&lxm, &w);
        let got = Quartic::new(&s, &freq, &wila);

        // sigma is so crazy that its epsilon is much higher. this would be a
        // good use for the approx relative difference
        assert_abs_diff_eq!(got.sigma, test.want.sigma, epsilon = 1.0);
        let got = Quartic {
            sigma: test.want.sigma,
            ..got.clone()
        };
        if abs_diff_ne!(got, test.want, epsilon = 2e-5) {
            println!("got\n{got}");
            println!("want\n{}", test.want);
            println!("diff\n{}", got - test.want.clone());
            panic!("{} failed", test.infile);
        }
    }
}
