use approx::assert_abs_diff_eq;

use super::*;

#[derive(Clone)]
struct Test {
    infile: String,
    fort15: String,
    lxm: Dmat,
    lx: Dmat,
    harm: Vec<f64>,
    eps: f64,
}

impl Test {
    /// create a new [Test] with its directory in testfiles and the dimensions
    /// of its `lxm` matrix
    fn new(dir: &'static str, lxm: (usize, usize), eps: f64) -> Self {
        let start = Path::new("testfiles");
        Self {
            infile: String::from(
                start.join(dir).join("spectro.in").to_str().unwrap(),
            ),
            fort15: String::from(
                start.join(dir).join("fort.15").to_str().unwrap(),
            ),
            lxm: load_dmat(
                start.join(dir).join("lxm").to_str().unwrap(),
                lxm.0,
                lxm.1,
            ),
            lx: load_dmat(
                start.join(dir).join("lx").to_str().unwrap(),
                lxm.0,
                lxm.1,
            ),
            harm: load_vec(start.join(dir).join("harm").to_str().unwrap()),
            eps,
        }
    }
}

#[macro_export]
macro_rules! check_eigen {
    ($got: expr, $want: expr, $eps: expr, $label: expr, $infile: expr) => {
        $crate::tests::lxm::check_eigen(
            $got,
            $want,
            $eps,
            $label,
            &format!("'{}', {}:{}:{}", $infile, file!(), line!(), column!(),),
        )
    };
}

/// check that two matrices, `got` and `want`, differ by at most `eps`, but
/// allowing the sign of whole columns to differ in line with non-unique
/// eigendecompositions. Really, columns could be allowed to differ by any
/// factor, not just -1, but I haven't encountered that situation yet.
pub(crate) fn check_eigen<'a, S, R: nalgebra::Dim, C: nalgebra::Dim>(
    got: &'a Matrix<f64, R, C, S>,
    want: &'a Matrix<f64, R, C, S>,
    eps: f64,
    label: &str,
    infile: &str,
) where
    &'a Matrix<f64, R, C, S>: Sub<Output = Matrix<f64, R, C, S>>,
    S: Storage<f64, R, C>,
    DefaultAllocator: Allocator<R, C>,
    DefaultAllocator: nalgebra::allocator::Allocator<R>,
{
    assert_eq!(got.shape(), want.shape());
    let (_, cols) = got.shape();
    for col in 0..cols {
        let g = got.column(col);
        let w = want.column(col);
        if abs_diff_ne!(g, w, epsilon = eps) {
            let diff1 = (g - w).abs().max();
            let g = -1.0 * g;
            // stupid nalgebra typing issue
            let w = 1.0 * w;
            if abs_diff_ne!(g, w, epsilon = eps) {
                println!("got\n{got:.8}");
                println!("want\n{want:.8}");
                let diff = got - want;
                let diff2 = (g - w).abs().max();
                println!("diff\n{diff:.8}");
                println!("max diff = {:.2e}", diff1.min(diff2));
                panic!("{label} failed on column {col} of {infile}");
            }
        }
    }
}

#[test]
fn asym() {
    let tests = [
        Test::new("h2o", (9, 9), 5e-9),
        Test::new("h2co", (12, 12), 5e-9),
        Test::new("c3h2", (15, 15), 5e-9),
        Test::new("c3hf", (15, 15), 5e-9),
        Test::new("c3hcn", (18, 18), 5e-9),
        Test::new("c3hcn010", (18, 18), 5e-9),
    ];

    for test in Vec::from(&tests[..]) {
        let s = Spectro::load(&test.infile);
        let fc2 = load_fc2(&test.fort15, s.n3n);
        let fc2 = s.rot2nd(&fc2);
        let fc2 = FACT2 * fc2;
        let w = s.geom.weights();
        let sqm: Vec<_> = w.iter().map(|w| 1.0 / w.sqrt()).collect();
        let fxm = s.form_sec(fc2, &sqm);

        let (harms, lxm) = symm_eigen_decomp(fxm, true);
        let lx = s.make_lx(&sqm, &lxm);

        check_vec!(
            to_wavenumbers(&harms),
            Dvec::from(test.harm),
            6e-6,
            test.infile
        );

        // only really care about the part with frequencies. there is more noise
        // in the rotations and translations, so this allows tightening epsilon
        let got = lxm.view((0, 0), (s.n3n, s.nvib));
        let want = test.lxm.view((0, 0), (s.n3n, s.nvib));

        // println!("{:.2e}", (got.clone() - want.clone()).max());
        check_mat!(
            &Dmat::from(got).abs(),
            &Dmat::from(want).abs(),
            test.eps,
            &test.infile
        );
        // assert_abs_diff_eq!(got, want, epsilon = 2e-9);

        let got = lx.view((0, 0), (s.n3n, s.nvib)).abs();
        let want = test.lx.view((0, 0), (s.n3n, s.nvib)).abs();
        // println!("{:.2e}", (got.clone() - want.clone()).max());
        // a little looser, but I guess that's from mass differences since these
        // are multiplied by 1/√w
        assert_abs_diff_eq!(got, want, epsilon = test.eps);
    }
}

/// testing symm_eigen_decomp LXM against LXM printed in spectro2.out, using the
/// spectro2.out FXM as input
#[test]
fn c3hf_lxm() {
    let fxm = load_dmat("testfiles/c3hf/fort_fxm", 15, 15);
    let want = load_dmat("testfiles/c3hf/fort_lxm", 15, 15);
    let (_, got) = symm_eigen_decomp(fxm, true);

    let got = got.view((0, 0), (15, 9));
    let want = want.view((0, 0), (15, 9));
    // all the precision I got from spectro2.out
    check_eigen!(&Dmat::from(got), &Dmat::from(want), 2.4e-7, "lxm", "c3hf");
}

#[test]
fn c3hcn_lxm() {
    let fxm = load_lower_triangle("testfiles/c3hcn/fort_fxm", 18);
    let want = load_dmat("testfiles/c3hcn/fort_lxm", 18, 18);
    let (_, got) = symm_eigen_decomp(fxm, true);

    let got = got.view((0, 0), (18, 12));
    let want = want.view((0, 0), (18, 12));
    // all the precision I got from spectro2.out
    check_eigen!(&Dmat::from(got), &Dmat::from(want), 1.32e-6, "lxm", "c3hcn");
}

#[test]
fn sym() {
    let tests = [
        Test::new("nh3", (12, 6), 2e-9),
        Test::new("ph3", (12, 6), 2e-9),
        Test::new("bipy", (21, 21), 2e-9),
        Test::new("c2h-", (9, 9), 2e-9),
        Test::new("hmgnc", (12, 12), 2e-9),
    ];

    for test in Vec::from(&tests[..]) {
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

        // had to loosen this for one of the imaginary frequencies in nh3 and
        // then a little more after using fortran symm eigen
        check_vec!(
            to_wavenumbers(&harms),
            Dvec::from(test.harm),
            5.9e-5,
            &test.infile
        );

        // only really care about the part with frequencies. there is more noise
        // in the rotations and translations, so this allows tightening epsilon
        let got = lxm.view((0, 0), (s.n3n, s.nvib));
        let want = test.lxm.view((0, 0), (s.n3n, s.nvib));

        check_mat!(&got.abs(), &want.abs(), test.eps, &test.infile);

        let got = lx.view((0, 0), (s.n3n, s.nvib)).abs();
        let want = test.lx.view((0, 0), (s.n3n, s.nvib)).abs();
        check_mat!(&got, &want, test.eps, &test.infile);
    }
}

/// test that the lxm and lx matrices going into bdegnl are correct, at least as
/// correct as the other lxm tests, to within signs. this part is really tested
/// by every asym test, but I'm trying to debug ph3
#[test]
fn pre_bdegnl() {
    let s = Spectro::load("testfiles/ph3/spectro.in");
    let fc2 = load_fc2("testfiles/ph3/fort.15", s.n3n);
    let fc2 = s.rot2nd(&fc2);
    let fc2 = FACT2 * fc2;
    let w = s.geom.weights();
    let sqm: Vec<_> = w.iter().map(|w| 1.0 / w.sqrt()).collect();
    let fxm = s.form_sec(fc2, &sqm);
    let (_, lxm) = symm_eigen_decomp(fxm, true);
    let lx = s.make_lx(&sqm, &lxm);

    let got = lxm.view((0, 0), (s.n3n, s.nvib));
    let want = load_dmat("testfiles/ph3/pre_bdegnl_lxm", 12, 12);
    let want = want.view((0, 0), (s.n3n, s.nvib));

    // println!("{:.2e}", (got.clone() - want.clone()).max());
    check_mat!(&got.abs(), &want.abs(), 2e-9, "ph3");

    let got = lx.view((0, 0), (s.n3n, s.nvib)).abs();
    let want = load_dmat("testfiles/ph3/pre_bdegnl_lx", 12, 12);
    let want = want.view((0, 0), (s.n3n, s.nvib));
    check_mat!(&got.abs(), &want.abs(), 2e-9, "ph3");
}

/// this is testing that `bdegnl` is working properly with values input from the
/// Fortran version
#[test]
fn bdegnl() {
    // ph3 values from vibfx.f:169
    let mut lxm = load_dmat("testfiles/ph3/pre_bdegnl_lxm", 12, 12);
    let mut lx = load_dmat("testfiles/ph3/pre_bdegnl_lx", 12, 12);
    let freq = na::dvector![
        2_437.002_438_242_96,
        2_437.002_030_016_922,
        2428.1546959531433,
        1145.7614713948005,
        1145.7567008646774,
        1015.4828482214278,
        0.046_516_934_814_435_32,
        0.025_278_773_879_589_64,
        0.015292949295790762,
        0.007_726_679_933_684_167,
        0.011089653435072953,
        0.014380890514910207
    ];
    use std::str::FromStr;
    let s = Spectro {
        nvib: 6,
        axis_order: 3,
        geom: Molecule::from_str(
            "
H  1.18656578  0.00000000  0.69757310
P  0.00000000 -0.00000000 -0.06809295
H -0.59328292  1.02759614  0.69757310
H -0.59328292 -1.02759614  0.69757310
",
        )
        .unwrap(),
        ..Spectro::default()
    };
    let w = s.geom.weights();

    s.bdegnl(&freq, &mut lxm, &w, &mut lx);

    let got = lxm.view((0, 0), (12, 6)).abs();
    let want = load_dmat("testfiles/ph3/lxm", 12, 6).abs();
    // println!("{:.2e}", (got.clone() - want.clone()).max());
    assert_abs_diff_eq!(got, want, epsilon = 1e-10);

    let got = lx.view((0, 0), (12, 6)).abs();
    let want = load_dmat("testfiles/ph3/lx", 12, 6).abs();
    // println!("{:.2e}", (got.clone() - want.clone()).max());
    assert_abs_diff_eq!(got, want, epsilon = 1e-10);
}

/// check if I can get the same FXM without absolute value using the same
/// geometry as the fortran code.
#[test]
fn ph3_fxm() {
    use std::str::FromStr;
    let s = Spectro {
        n3n: 12,
        geom: Molecule::from_str(
            "
H  1.18656578  0.00000000  0.69757313
P  0.00000000 -0.00000000 -0.06809295
H -0.59328292  1.02759614  0.69757310
H -0.59328292 -1.02759614  0.69757310
",
        )
        .unwrap(),
        axes: nalgebra::matrix![
        0.999999986417,0.000000000000,0.000164823265;
        -0.000164823265,0.000000000000,0.999999986417;
        0.000000000000,1.000000000000,0.000000000000;
        ],
        ..Spectro::default()
    };
    let fort15 = "testfiles/ph3/fort.15";
    let fc2 = load_fc2(fort15, s.n3n);
    let fc2 = s.rot2nd(&fc2);
    let fc2 = FACT2 * fc2;
    let w = s.geom.weights();
    let sqm: Vec<_> = w.iter().map(|w| 1.0 / w.sqrt()).collect();
    let got = s.form_sec(fc2, &sqm);

    let want = load_lower_triangle("testfiles/ph3/fxm", 12);

    // println!("{:.2e}", (got.clone() - want.clone()).max());
    check_mat(&got, &want, 1e-7, "ph3");
}

/// all of these values are from spectro2.out
#[test]
fn c3hcn_fxm() {
    use std::str::FromStr;
    let s = Spectro {
        n3n: 18,
        geom: Molecule::from_str(
            "
H      2.0345490     -1.5884146      0.0000000
C      1.5163524     -0.6441862      0.0000000
C      1.5721454      0.7755306      0.0000000
C      0.3627860      0.0129312      0.0000000
C     -1.0464358      0.0002536      0.0000000
N     -2.2072758     -0.0095340      0.0000000
",
        )
        .unwrap(),
        axes: nalgebra::matrix![
        0.00000000,0.00000000,1.00000000;
        0.00005289,-1.00000000,0.00000000;
        1.00000000,0.00005289,0.00000000;
                ],
        ..Spectro::default()
    };
    let fort15 = "testfiles/c3hcn/fort.15";
    let fc2 = load_fc2(fort15, s.n3n);
    let fc2 = s.rot2nd(&fc2);
    let fc2 = FACT2 * fc2;
    let w = s.geom.weights();
    let sqm: Vec<_> = w.iter().map(|w| 1.0 / w.sqrt()).collect();
    let got = s.form_sec(fc2, &sqm);

    let want = load_lower_triangle("testfiles/c3hcn/fort_fxm", 18);

    // println!("{:.2e}", (got.clone() - want.clone()).max());
    check_mat(&got, &want, 7e-8, "c3hcn");
}

/// build on c3hcn_fxm to see if I can take *that* fxm forward to get a good lxm
#[test]
fn c3hcn_fxm_lxm() {
    use std::str::FromStr;
    let s = Spectro {
        n3n: 18,
        geom: Molecule::from_str(
            "
H      2.0345490     -1.5884146      0.0000000
C      1.5163524     -0.6441862      0.0000000
C      1.5721454      0.7755306      0.0000000
C      0.3627860      0.0129312      0.0000000
C     -1.0464358      0.0002536      0.0000000
N     -2.2072758     -0.0095340      0.0000000
",
        )
        .unwrap(),
        axes: nalgebra::matrix![
        0.00000000,0.00000000,1.00000000;
        0.00005289,-1.00000000,0.00000000;
        1.00000000,0.00005289,0.00000000;
                ],
        ..Spectro::default()
    };
    let fort15 = "testfiles/c3hcn/fort.15";
    let fc2 = load_fc2(fort15, s.n3n);
    let fc2 = s.rot2nd(&fc2);
    let fc2 = FACT2 * fc2;
    let w = s.geom.weights();
    let sqm: Vec<_> = w.iter().map(|w| 1.0 / w.sqrt()).collect();
    let fxm = s.form_sec(fc2, &sqm);

    let want = load_dmat("testfiles/c3hcn/fort_lxm", 18, 18);
    let (_, got) = symm_eigen_decomp(fxm, true);

    let got = got.view((0, 0), (18, 12));
    let want = want.view((0, 0), (18, 12));
    // all the precision I got from spectro2.out
    check_eigen!(&Dmat::from(got), &Dmat::from(want), 5e-8, "lxm", "c3hcn");
}
/// check if I can get the same LXM without absolute value using the same
/// geometry as the fortran code. this is a copy-paste of fxm but with the
/// additional line of the eigendecomposition. this is basically a test of
/// `symm_eigen_decomp`
#[test]
fn lxm() {
    use std::str::FromStr;
    let s = Spectro {
        n3n: 12,
        geom: Molecule::from_str(
            "
H  1.18656578  0.00000000  0.69757313
P  0.00000000 -0.00000000 -0.06809295
H -0.59328292  1.02759614  0.69757310
H -0.59328292 -1.02759614  0.69757310
",
        )
        .unwrap(),
        axes: nalgebra::matrix![
        0.999999986417,0.000000000000,0.000164823265;
        -0.000164823265,0.000000000000,0.999999986417;
        0.000000000000,1.000000000000,0.000000000000;
        ],
        ..Spectro::default()
    };
    let fort15 = "testfiles/ph3/fort.15";
    let fc2 = load_fc2(fort15, s.n3n);
    let fc2 = s.rot2nd(&fc2);
    let fc2 = FACT2 * fc2;
    let w = s.geom.weights();
    let sqm: Vec<_> = w.iter().map(|w| 1.0 / w.sqrt()).collect();
    let fxm = s.form_sec(fc2, &sqm);

    let (_harms, lxm) = utils::linalg::symm_eigen_decomp(fxm, true);

    let want = load_dmat("testfiles/ph3/pre_bdegnl_lxm", 12, 12);
    let got = lxm.view((0, 0), (12, 6));
    let want = want.view((0, 0), (12, 6));

    // println!("{:.2e}", (got.clone() - want.clone()).max());
    check_eigen(&Dmat::from(got), &Dmat::from(want), 1e-7, "lxm", "ph3");
}
