use super::*;

#[derive(Clone)]
struct Test {
    infile: String,
    fort15: String,
    wila: Dmat,
    zmat: Tensor3,
    zmat_eps: f64,
    wila_eps: f64,
}

impl Test {
    fn new(
        dir: &'static str,
        rows: usize,
        cols: usize,
        zmat_eps: f64,
        wila_eps: f64,
    ) -> Self {
        let start = Path::new("testfiles");
        Self {
            infile: String::from(
                start.join(dir).join("spectro.in").to_str().unwrap(),
            ),
            fort15: String::from(
                start.join(dir).join("fort.15").to_str().unwrap(),
            ),
            wila: load_dmat(
                start.join(dir).join("wila").to_str().unwrap(),
                rows,
                cols,
            ),
            zmat: Tensor3::load(start.join(dir).join("zmat").to_str().unwrap()),
            zmat_eps,
            wila_eps,
        }
    }
}

#[test]
fn asym() {
    let tests = [
        // eps increases with mass, which I guess is from mass dependence of lxm
        // and also wila itself
        Test::new("h2o", 3, 6, 1.53e-10, 7.6e-7),
        Test::new("h2co", 6, 6, 1.41e-9, 1.8e-6),
        Test::new("c3h2", 9, 6, 5.78e-9, 2.8e-6),
        Test::new("c3hf", 9, 6, 9.85e-10, 4.2e-6),
        Test::new("c3hcn", 12, 6, 8.39e-10, 6.6e-6),
        Test::new("c3hcn010", 12, 6, 8.39e-10, 6.6e-6),
    ];
    for test in Vec::from(&tests[..]) {
        let s = Spectro::load(&test.infile);
        let fc2 = load_fc2(&test.fort15, s.n3n);
        let fc2 = s.rot2nd(&fc2);
        let fc2 = FACT2 * fc2;
        let w = s.geom.weights();
        let sqm: Vec<_> = w.iter().map(|w: &f64| 1.0 / w.sqrt()).collect();
        let fxm = s.form_sec(fc2, &sqm);
        let (_harms, lxm) = symm_eigen_decomp(fxm, true);

        let (zmat, wila) = s.zeta(&lxm, &w);

        check_tens!(
            &zmat.abs(),
            &test.zmat.abs(),
            test.zmat_eps,
            "zmat",
            &test.infile
        );

        check_mat!(&wila.abs(), &test.wila.abs(), test.wila_eps, &test.infile);
    }
}

#[test]
fn sym() {
    let tests = [
        Test::new("nh3", 6, 6, 3.52e-10, 8.73e-7),
        Test::new("ph3", 6, 6, 3.52e-10, 8.73e-7),
        Test::new("bipy", 15, 6, 1.02e-9, 8.73e-7),
    ];
    for test in Vec::from(&tests[..]) {
        let s = Spectro::load(&test.infile);
        let fc2 = load_fc2(&test.fort15, s.n3n);
        let fc2 = s.rot2nd(&fc2);
        let fc2 = FACT2 * fc2;
        let w = s.geom.weights();
        let sqm: Vec<_> = w.iter().map(|w: &f64| 1.0 / w.sqrt()).collect();
        let fxm = s.form_sec(fc2, &sqm);
        let (harms, mut lxm) = symm_eigen_decomp(fxm, true);
        let freq = to_wavenumbers(&harms);
        let mut lx = s.make_lx(&sqm, &lxm);
        s.bdegnl(&freq, &mut lxm, &w, &mut lx);

        let (zmat, wila) = s.zeta(&lxm, &w);

        check_tens!(
            &zmat.abs(),
            &test.zmat.abs(),
            test.zmat_eps,
            "zmat",
            &test.infile
        );

        check_mat!(&wila.abs(), &test.wila.abs(), test.wila_eps, &test.infile);
    }
}
