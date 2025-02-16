use approx::{abs_diff_eq, assert_abs_diff_eq};

use crate::{consts::FACT2, resonance::Restst, *};

use super::load_dmat;

#[derive(Clone)]
struct Test {
    infile: String,
    fort15: String,
    fort30: String,
    fort40: String,
    xcnst: Dmat,
    gcnst: Dmat,

    /// xcnst after only the nondeg-nondeg part
    xcnst_nn: Dmat,

    /// xcnst after the nondeg-deg part
    xcnst_nd: Dmat,

    /// xcnst after the first deg-deg part
    xcnst_dd: Dmat,

    /// gcnst after the first part
    gcnst1: Dmat,

    e0: f64,
    xcnst_eps: f64,
    gcnst_eps: f64,
}

impl Test {
    fn new(
        dir: &'static str,
        nvib: usize,
        e0: f64,
        xcnst_eps: f64,
        gcnst_eps: f64,
    ) -> Self {
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
            fort40: String::from(
                start.join(dir).join("fort.40").to_str().unwrap(),
            ),
            xcnst: load_dmat(start.join(dir).join("xcnst"), nvib, nvib),
            xcnst_nn: load_dmat(start.join(dir).join("xcnst_nn"), nvib, nvib),
            xcnst_nd: load_dmat(start.join(dir).join("xcnst_nd"), nvib, nvib),
            xcnst_dd: load_dmat(start.join(dir).join("xcnst_dd"), nvib, nvib),
            gcnst: load_dmat(start.join(dir).join("gcnst"), nvib, nvib),
            gcnst1: load_dmat(start.join(dir).join("gcnst1"), nvib, nvib),
            e0,
            xcnst_eps,
            gcnst_eps,
        }
    }
}

macro_rules! check {
    ($got: expr, $want: expr, $eps: expr, $label: expr, $infile: expr) => {
        if !abs_diff_eq!($got, $want, epsilon = $eps) {
            println!("got\n{:.6}", $got);
            println!("want\n{:.6}", $want);
            println!("diff\n={:.6}", $got.clone() - $want.clone());
            println!(
                "max diff = {:.2e}",
                ($got.clone() - $want.clone()).abs().max()
            );
            assert!(false, "{} differs on {}", $label, $infile);
        }
    };
}

/// setup for calling xcals
#[allow(clippy::type_complexity)]
fn setup(
    test: &Test,
    s: &Spectro,
) -> (
    Dvec,
    tensor::Tensor3<f64>,
    Dmat,
    F3qcm,
    F4qcm,
    Vec<Fermi1>,
    Vec<Fermi2>,
    Vec<Mode>,
) {
    let fc2 = load_fc2(&test.fort15, s.n3n);
    let fc2 = s.rot2nd(&fc2);
    let fc2 = FACT2 * fc2;
    let w = s.geom.weights();
    let sqm: Vec<_> = w.iter().map(|w| 1.0 / w.sqrt()).collect();
    let fxm = s.form_sec(fc2, &sqm);
    let (harms, mut lxm) = utils::linalg::symm_eigen_decomp(fxm, true);
    let freq = to_wavenumbers(&harms);
    let mut lx = s.make_lx(&sqm, &lxm);
    s.bdegnl(&freq, &mut lxm, &w, &mut lx);
    let (zmat, wila) = s.zeta(&lxm, &w);
    let f3x = load_fc3(&test.fort30, s.n3n);
    let f3x = s.rot3rd(f3x);
    let f3qcm = force3(s.n3n, f3x, &lx, s.nvib, &freq);
    let f4x = load_fc4(&test.fort40, s.n3n);
    let f4x = s.rot4th(f4x);
    let f4qcm = force4(s.n3n, &f4x, &lx, s.nvib, &freq);
    let Restst {
        coriolis: _,
        fermi1,
        fermi2,
        darling: _,
        states: _,
        modes,
        ifunda: _,
        iovrtn: _,
        icombn: _,
    } = Restst::new(s, &zmat, &f3qcm, &freq);
    (freq, zmat, wila, f3qcm, f4qcm, fermi1, fermi2, modes)
}

#[test]
fn nondeg_nondeg() {
    let tests = [
        Test::new("nh3", 6, 24.716378286389887, 1e-11, 0.0),
        Test::new("ph3", 6, 20.748849036017717, 1e-11, 0.0),
        Test::new("bipy", 15, 32.906_770_783_666_87, 1e-11, 0.0),
        Test::new("c2h-", 4, -1.0534319575869713, 1e-11, 6e-12),
        Test::new("hmgnc", 7, -0.12246977241439683, 1e-11, 6e-12),
    ];
    for test in Vec::from(&tests[..]) {
        let s = Spectro::load(&test.infile);
        let (freq, zmat, wila, f3qcm, f4qcm, fermi1, fermi2, modes) =
            setup(&test, &s);
        let (ia, _ib, _n2dmm, i1mode, _i2mode, _ixyz) =
            s.setup_xcals(&modes, &wila);
        let (_ifrmchk, ifrm1, ifrm2) = s.make_fermi_checks(&fermi1, &fermi2);
        let mut got = Dmat::zeros(s.nvib, s.nvib);

        s.nondeg_nondeg(
            &i1mode, &f4qcm, &f3qcm, &freq, &ifrm2, ia, &zmat, &ifrm1, &mut got,
        );
        check!(
            &got,
            &test.xcnst_nn,
            test.xcnst_eps,
            "nondeg-nondeg",
            &test.infile
        );
    }
}

#[test]
fn nondeg_deg() {
    let tests = [
        Test::new("nh3", 6, 24.716378286389887, 4e-10, 0.0),
        Test::new("ph3", 6, 20.748849036017717, 1e-11, 0.0),
        Test::new("bipy", 15, 32.906_770_783_666_87, 1.52e-11, 0.0),
        Test::new("hmgnc", 7, -0.12246977241439683, 1e-11, 6e-12),
    ];
    for test in Vec::from(&tests[..]) {
        let s = Spectro::load(&test.infile);
        let (freq, zmat, wila, f3qcm, f4qcm, fermi1, fermi2, modes) =
            setup(&test, &s);
        let (ia, ib, _n2dmm, i1mode, i2mode, _ixyz) =
            s.setup_xcals(&modes, &wila);
        let (_ifrmchk, ifrm1, ifrm2) = s.make_fermi_checks(&fermi1, &fermi2);
        let mut got = Dmat::zeros(s.nvib, s.nvib);

        s.nondeg_nondeg(
            &i1mode, &f4qcm, &f3qcm, &freq, &ifrm2, ia, &zmat, &ifrm1, &mut got,
        );

        s.nondeg_deg(
            &i1mode, &i2mode, &f4qcm, &f3qcm, &freq, &ifrm2, ib, &zmat,
            &mut got,
        );

        check!(
            &got,
            &test.xcnst_nd,
            test.xcnst_eps,
            "nondeg-deg",
            &test.infile
        );
    }
}

#[test]
fn deg_deg1() {
    let tests = [
        Test::new("nh3", 6, 24.716378286389887, 5e-10, 0.0),
        Test::new("ph3", 6, 20.748849036017717, 1e-11, 0.0),
        Test::new("bipy", 15, 32.906_770_783_666_87, 1.52e-11, 0.0),
        Test::new("hmgnc", 7, -0.12246977241439683, 1e-11, 6e-12),
    ];
    for test in Vec::from(&tests[..]) {
        let s = Spectro::load(&test.infile);
        let (freq, zmat, wila, f3qcm, f4qcm, fermi1, fermi2, modes) =
            setup(&test, &s);
        let (ia, ib, _n2dmm, i1mode, i2mode, _ixyz) =
            s.setup_xcals(&modes, &wila);
        let (_ifrmchk, ifrm1, ifrm2) = s.make_fermi_checks(&fermi1, &fermi2);
        let mut got = Dmat::zeros(s.nvib, s.nvib);

        s.nondeg_nondeg(
            &i1mode, &f4qcm, &f3qcm, &freq, &ifrm2, ia, &zmat, &ifrm1, &mut got,
        );

        s.nondeg_deg(
            &i1mode, &i2mode, &f4qcm, &f3qcm, &freq, &ifrm2, ib, &zmat,
            &mut got,
        );

        xcals::deg_deg1(
            &i2mode, &f4qcm, &freq, &i1mode, &f3qcm, &ifrm1, &mut got,
        );

        check!(
            &got,
            &test.xcnst_dd,
            test.xcnst_eps,
            "deg-deg",
            &test.infile
        );
    }
}

#[test]
fn gcnst1() {
    let tests = [
        Test::new("nh3", 6, 24.716378286389887, 5e-10, 3e-9),
        Test::new("ph3", 6, 20.748849036017717, 1e-11, 9.32e-12),
        Test::new("bipy", 15, 32.906_770_783_666_87, 1e-11, 6e-12),
        Test::new("c2h-", 4, -1.0534319575869713, 1e-11, 6e-12),
        Test::new("hmgnc", 7, -1.0534319575869713, 1e-11, 6e-12),
    ];
    for test in Vec::from(&tests[..]) {
        let s = Spectro::load(&test.infile);
        let (freq, zmat, wila, f3qcm, f4qcm, fermi1, fermi2, modes) =
            setup(&test, &s);
        let (ia, _ib, n2dm, i1mode, i2mode, _ixyz) =
            s.setup_xcals(&modes, &wila);
        let (_ifrmchk, ifrm1, _ifrm2) = s.make_fermi_checks(&fermi1, &fermi2);

        let mut got = Dmat::zeros(s.nvib, s.nvib);

        s.gcnst1(
            n2dm, &i2mode, &f4qcm, &freq, &i1mode, &f3qcm, &ifrm1, ia, &zmat,
            &mut got,
        );

        check!(&got, &test.gcnst1, test.gcnst_eps, "gcnst", &test.infile);
    }
}

macro_rules! warn {
    ($msg: expr) => {
        eprintln!(
            "\nwarning: {} --> {}:{}:{}",
            $msg,
            file!(),
            line!(),
            column!()
        );
    };
}

#[test]
fn sym() {
    let tests = [
        Test::new("nh3", 6, 24.716378286389887, 5e-10, 3e-9),
        Test::new("ph3", 6, 20.748849036017717, 1e-11, 9.32e-12),
        {
            warn!("high bipy gcnst eps");
            Test::new("bipy", 15, 32.906_770_783_666_87, 1.52e-11, 1.01)
        },
        Test::new("c2h-", 4, -1.0534319575869713, 1e-11, 6e-12),
        Test::new("hmgnc", 7, -0.12246977241439683, 1e-11, 6e-12),
    ];
    for test in Vec::from(&tests[..]) {
        let s = Spectro::load(&test.infile);
        let (freq, zmat, wila, f3qcm, f4qcm, fermi1, fermi2, modes) =
            setup(&test, &s);
        let (xcnst, gcnst, e0) = s.xcals(
            &f4qcm, &freq, &f3qcm, &zmat, &fermi1, &fermi2, &modes, &wila,
        );
        assert_abs_diff_eq!(e0, test.e0, epsilon = 2e-9);
        check!(&xcnst, &test.xcnst, test.xcnst_eps, "xcnst", &test.infile);
        check!(&gcnst, &test.gcnst, test.gcnst_eps, "gcnst", &test.infile);
    }
}
