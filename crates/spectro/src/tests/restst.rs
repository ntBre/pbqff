use crate::{
    consts::FACT2,
    resonance::{Darling, Restst},
    state::{I1states, I2states, StatePartition},
    *,
};

struct Test {
    infile: String,
    fort15: String,
    fort30: String,
    want: Restst,
}

impl Test {
    fn new(dir: &'static str, want: Restst) -> Self {
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
            want,
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
        let (harms, mut lxm) = utils::linalg::symm_eigen_decomp(fxm, true);
        let freq = to_wavenumbers(&harms);
        let mut lx = s.make_lx(&sqm, &lxm);
        if s.rotor.is_sym_top() {
            s.bdegnl(&freq, &mut lxm, &w, &mut lx);
        }
        let (zmat, _wila) = s.zeta(&lxm, &w);
        let f3x = load_fc3(&test.fort30, s.n3n);
        let f3x = s.rot3rd(f3x);
        let f3qcm = force3(s.n3n, f3x, &lx, s.nvib, &freq);
        let got = Restst::new(&s, &zmat, &f3qcm, &freq);
        assert_eq!(got.coriolis, test.want.coriolis, "{}", test.infile);
        assert_eq!(got.fermi1, test.want.fermi1, "{}", test.infile);
        // for i in 0..got.fermi2.len() {
        //     println!("{} vs {}", got.fermi2[i], test.want.fermi2[i]);
        // }
        assert_eq!(got.fermi2.len(), test.want.fermi2.len());
        assert_eq!(got.fermi2, test.want.fermi2, "{}", test.infile);
        assert_eq!(got.darling, test.want.darling, "{}", test.infile);
        assert_eq!(got.states.len(), test.want.states.len(), "{}", test.infile);
        let StatePartition {
            i1sts,
            i2sts,
            i3sts,
        } = State::partition(&got.states);
        let StatePartition {
            i1sts: want_i1,
            i2sts: want_i2,
            i3sts: want_i3,
        } = State::partition(&test.want.states);
        assert_eq!(i1sts.len(), want_i1.len());
        assert_eq!(
            i1sts,
            want_i1,
            "got\n{}\n\nwant{}",
            I1states(i1sts.clone()),
            I1states(want_i1.clone())
        );
        assert_eq!(i2sts.len(), want_i2.len());
        for i in 0..i2sts.len() {
            if i2sts[i] != want_i2[i] {
                println!("error on {i}");
            }
        }
        assert_eq!(
            i2sts,
            want_i2,
            "got\n{}\n\nwant{}",
            I2states(i2sts.clone()),
            I2states(want_i2.clone())
        );
        assert_eq!(i3sts, want_i3, "{}", test.infile);
        assert_eq!(got.states, test.want.states, "{}", test.infile);
        assert_eq!(got.modes, test.want.modes, "{}", test.infile);
        assert_eq!(got.ifunda, test.want.ifunda, "{}", test.infile);
        assert_eq!(got.iovrtn, test.want.iovrtn, "{}", test.infile);
        assert_eq!(got.icombn, test.want.icombn, "{}", test.infile);
        assert_eq!(got, test.want, "{}", test.infile);
    }
}

#[test]
fn asym() {
    use resonance::Axis::*;
    use state::State::*;
    use Mode::*;
    let tests = [
        Test::new(
            "h2o",
            Restst {
                coriolis: vec![],
                fermi1: vec![],
                fermi2: vec![],
                darling: vec![Darling::new(1, 0)],
                states: vec![
                    I1st(vec![0, 0, 0]),
                    I1st(vec![1, 0, 0]),
                    I1st(vec![0, 1, 0]),
                    I1st(vec![0, 0, 1]),
                    I1st(vec![2, 0, 0]),
                    I1st(vec![0, 2, 0]),
                    I1st(vec![0, 0, 2]),
                    I1st(vec![1, 1, 0]),
                    I1st(vec![1, 0, 1]),
                    I1st(vec![0, 1, 1]),
                ],
                modes: vec![I1(0), I1(1), I1(2)],
                ifunda: vec![1, 2, 3],
                iovrtn: vec![4, 5, 6],
                icombn: vec![0, 7, 0, 8, 9, 0],
            },
        ),
        Test::new(
            "h2co",
            Restst {
                coriolis: vec![Coriolis::new(5, 4, A)],
                fermi1: vec![Fermi1::new(3, 1)],
                fermi2: vec![Fermi2::new(4, 2, 0)],
                darling: vec![Darling::new(1, 0), Darling::new(5, 4)],
                states: vec![
                    I1st(vec![0, 0, 0, 0, 0, 0]),
                    I1st(vec![1, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 1, 0, 0, 0, 0]),
                    I1st(vec![0, 0, 1, 0, 0, 0]),
                    I1st(vec![0, 0, 0, 1, 0, 0]),
                    I1st(vec![0, 0, 0, 0, 1, 0]),
                    I1st(vec![0, 0, 0, 0, 0, 1]),
                    I1st(vec![2, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 2, 0, 0, 0, 0]),
                    I1st(vec![0, 0, 2, 0, 0, 0]),
                    I1st(vec![0, 0, 0, 2, 0, 0]),
                    I1st(vec![0, 0, 0, 0, 2, 0]),
                    I1st(vec![0, 0, 0, 0, 0, 2]),
                    I1st(vec![1, 1, 0, 0, 0, 0]),
                    I1st(vec![1, 0, 1, 0, 0, 0]),
                    I1st(vec![1, 0, 0, 1, 0, 0]),
                    I1st(vec![1, 0, 0, 0, 1, 0]),
                    I1st(vec![1, 0, 0, 0, 0, 1]),
                    I1st(vec![0, 1, 1, 0, 0, 0]),
                    I1st(vec![0, 1, 0, 1, 0, 0]),
                    I1st(vec![0, 1, 0, 0, 1, 0]),
                    I1st(vec![0, 1, 0, 0, 0, 1]),
                    I1st(vec![0, 0, 1, 1, 0, 0]),
                    I1st(vec![0, 0, 1, 0, 1, 0]),
                    I1st(vec![0, 0, 1, 0, 0, 1]),
                    I1st(vec![0, 0, 0, 1, 1, 0]),
                    I1st(vec![0, 0, 0, 1, 0, 1]),
                    I1st(vec![0, 0, 0, 0, 1, 1]),
                ],
                modes: vec![I1(0), I1(1), I1(2), I1(3), I1(4), I1(5)],
                ifunda: vec![1, 2, 3, 4, 5, 6],
                iovrtn: vec![7, 8, 9, 10, 11, 12],
                icombn: vec![
                    0, 13, 0, 14, 18, 0, 15, 19, 22, 0, 16, 20, 23, 25, 0, 17,
                    21, 24, 26, 27, 0,
                ],
            },
        ),
        Test::new(
            "c3h2",
            Restst {
                coriolis: vec![
                    Coriolis::new(5, 4, A),
                    Coriolis::new(6, 5, A),
                    Coriolis::new(7, 4, C),
                    Coriolis::new(7, 5, B),
                    Coriolis::new(7, 6, C),
                    Coriolis::new(8, 6, B),
                    Coriolis::new(8, 7, A),
                ],
                fermi1: vec![
                    Fermi1::new(2, 0),
                    Fermi1::new(6, 2),
                    Fermi1::new(7, 2),
                    Fermi1::new(8, 2),
                ],
                fermi2: vec![],
                darling: vec![
                    Darling::new(1, 0),
                    Darling::new(5, 4),
                    Darling::new(6, 5),
                    Darling::new(7, 5),
                    Darling::new(7, 6),
                    Darling::new(8, 6),
                    Darling::new(8, 7),
                ],
                states: vec![
                    I1st(vec![0, 0, 0, 0, 0, 0, 0, 0, 0]),
                    I1st(vec![1, 0, 0, 0, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 1, 0, 0, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 0, 1, 0, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 0, 0, 1, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 0, 0, 0, 1, 0, 0, 0, 0]),
                    I1st(vec![0, 0, 0, 0, 0, 1, 0, 0, 0]),
                    I1st(vec![0, 0, 0, 0, 0, 0, 1, 0, 0]),
                    I1st(vec![0, 0, 0, 0, 0, 0, 0, 1, 0]),
                    I1st(vec![0, 0, 0, 0, 0, 0, 0, 0, 1]),
                    I1st(vec![2, 0, 0, 0, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 2, 0, 0, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 0, 2, 0, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 0, 0, 2, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 0, 0, 0, 2, 0, 0, 0, 0]),
                    I1st(vec![0, 0, 0, 0, 0, 2, 0, 0, 0]),
                    I1st(vec![0, 0, 0, 0, 0, 0, 2, 0, 0]),
                    I1st(vec![0, 0, 0, 0, 0, 0, 0, 2, 0]),
                    I1st(vec![0, 0, 0, 0, 0, 0, 0, 0, 2]),
                    I1st(vec![1, 1, 0, 0, 0, 0, 0, 0, 0]),
                    I1st(vec![1, 0, 1, 0, 0, 0, 0, 0, 0]),
                    I1st(vec![1, 0, 0, 1, 0, 0, 0, 0, 0]),
                    I1st(vec![1, 0, 0, 0, 1, 0, 0, 0, 0]),
                    I1st(vec![1, 0, 0, 0, 0, 1, 0, 0, 0]),
                    I1st(vec![1, 0, 0, 0, 0, 0, 1, 0, 0]),
                    I1st(vec![1, 0, 0, 0, 0, 0, 0, 1, 0]),
                    I1st(vec![1, 0, 0, 0, 0, 0, 0, 0, 1]),
                    I1st(vec![0, 1, 1, 0, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 1, 0, 1, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 1, 0, 0, 1, 0, 0, 0, 0]),
                    I1st(vec![0, 1, 0, 0, 0, 1, 0, 0, 0]),
                    I1st(vec![0, 1, 0, 0, 0, 0, 1, 0, 0]),
                    I1st(vec![0, 1, 0, 0, 0, 0, 0, 1, 0]),
                    I1st(vec![0, 1, 0, 0, 0, 0, 0, 0, 1]),
                    I1st(vec![0, 0, 1, 1, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 0, 1, 0, 1, 0, 0, 0, 0]),
                    I1st(vec![0, 0, 1, 0, 0, 1, 0, 0, 0]),
                    I1st(vec![0, 0, 1, 0, 0, 0, 1, 0, 0]),
                    I1st(vec![0, 0, 1, 0, 0, 0, 0, 1, 0]),
                    I1st(vec![0, 0, 1, 0, 0, 0, 0, 0, 1]),
                    I1st(vec![0, 0, 0, 1, 1, 0, 0, 0, 0]),
                    I1st(vec![0, 0, 0, 1, 0, 1, 0, 0, 0]),
                    I1st(vec![0, 0, 0, 1, 0, 0, 1, 0, 0]),
                    I1st(vec![0, 0, 0, 1, 0, 0, 0, 1, 0]),
                    I1st(vec![0, 0, 0, 1, 0, 0, 0, 0, 1]),
                    I1st(vec![0, 0, 0, 0, 1, 1, 0, 0, 0]),
                    I1st(vec![0, 0, 0, 0, 1, 0, 1, 0, 0]),
                    I1st(vec![0, 0, 0, 0, 1, 0, 0, 1, 0]),
                    I1st(vec![0, 0, 0, 0, 1, 0, 0, 0, 1]),
                    I1st(vec![0, 0, 0, 0, 0, 1, 1, 0, 0]),
                    I1st(vec![0, 0, 0, 0, 0, 1, 0, 1, 0]),
                    I1st(vec![0, 0, 0, 0, 0, 1, 0, 0, 1]),
                    I1st(vec![0, 0, 0, 0, 0, 0, 1, 1, 0]),
                    I1st(vec![0, 0, 0, 0, 0, 0, 1, 0, 1]),
                    I1st(vec![0, 0, 0, 0, 0, 0, 0, 1, 1]),
                ],
                modes: vec![
                    I1(0),
                    I1(1),
                    I1(2),
                    I1(3),
                    I1(4),
                    I1(5),
                    I1(6),
                    I1(7),
                    I1(8),
                ],
                ifunda: vec![1, 2, 3, 4, 5, 6, 7, 8, 9],
                iovrtn: vec![10, 11, 12, 13, 14, 15, 16, 17, 18],
                icombn: vec![
                    0, 19, 0, 20, 27, 0, 21, 28, 34, 0, 22, 29, 35, 40, 0, 23,
                    30, 36, 41, 45, 0, 24, 31, 37, 42, 46, 49, 0, 25, 32, 38,
                    43, 47, 50, 52, 0, 26, 33, 39, 44, 48, 51, 53, 54, 0,
                ],
            },
        ),
        Test::new(
            "c3hf",
            Restst {
                coriolis: vec![
                    Coriolis::new(5, 4, A),
                    Coriolis::new(5, 4, B),
                    Coriolis::new(8, 7, A),
                ],
                fermi1: vec![
                    Fermi1::new(4, 1),
                    Fermi1::new(5, 1),
                    Fermi1::new(6, 2),
                    Fermi1::new(7, 6),
                    Fermi1::new(8, 4),
                    Fermi1::new(8, 6),
                ],
                fermi2: vec![
                    Fermi2::new(2, 1, 0),
                    Fermi2::new(6, 3, 1),
                    Fermi2::new(6, 4, 1),
                    Fermi2::new(7, 2, 1),
                    Fermi2::new(7, 3, 1),
                    Fermi2::new(7, 4, 2),
                    Fermi2::new(7, 6, 2),
                    Fermi2::new(8, 7, 5),
                ],
                darling: vec![
                    Darling::new(5, 4),
                    Darling::new(6, 5),
                    Darling::new(8, 7),
                ],
                states: vec![
                    I1st(vec![0, 0, 0, 0, 0, 0, 0, 0, 0]),
                    I1st(vec![1, 0, 0, 0, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 1, 0, 0, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 0, 1, 0, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 0, 0, 1, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 0, 0, 0, 1, 0, 0, 0, 0]),
                    I1st(vec![0, 0, 0, 0, 0, 1, 0, 0, 0]),
                    I1st(vec![0, 0, 0, 0, 0, 0, 1, 0, 0]),
                    I1st(vec![0, 0, 0, 0, 0, 0, 0, 1, 0]),
                    I1st(vec![0, 0, 0, 0, 0, 0, 0, 0, 1]),
                    I1st(vec![2, 0, 0, 0, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 2, 0, 0, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 0, 2, 0, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 0, 0, 2, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 0, 0, 0, 2, 0, 0, 0, 0]),
                    I1st(vec![0, 0, 0, 0, 0, 2, 0, 0, 0]),
                    I1st(vec![0, 0, 0, 0, 0, 0, 2, 0, 0]),
                    I1st(vec![0, 0, 0, 0, 0, 0, 0, 2, 0]),
                    I1st(vec![0, 0, 0, 0, 0, 0, 0, 0, 2]),
                    I1st(vec![1, 1, 0, 0, 0, 0, 0, 0, 0]),
                    I1st(vec![1, 0, 1, 0, 0, 0, 0, 0, 0]),
                    I1st(vec![1, 0, 0, 1, 0, 0, 0, 0, 0]),
                    I1st(vec![1, 0, 0, 0, 1, 0, 0, 0, 0]),
                    I1st(vec![1, 0, 0, 0, 0, 1, 0, 0, 0]),
                    I1st(vec![1, 0, 0, 0, 0, 0, 1, 0, 0]),
                    I1st(vec![1, 0, 0, 0, 0, 0, 0, 1, 0]),
                    I1st(vec![1, 0, 0, 0, 0, 0, 0, 0, 1]),
                    I1st(vec![0, 1, 1, 0, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 1, 0, 1, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 1, 0, 0, 1, 0, 0, 0, 0]),
                    I1st(vec![0, 1, 0, 0, 0, 1, 0, 0, 0]),
                    I1st(vec![0, 1, 0, 0, 0, 0, 1, 0, 0]),
                    I1st(vec![0, 1, 0, 0, 0, 0, 0, 1, 0]),
                    I1st(vec![0, 1, 0, 0, 0, 0, 0, 0, 1]),
                    I1st(vec![0, 0, 1, 1, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 0, 1, 0, 1, 0, 0, 0, 0]),
                    I1st(vec![0, 0, 1, 0, 0, 1, 0, 0, 0]),
                    I1st(vec![0, 0, 1, 0, 0, 0, 1, 0, 0]),
                    I1st(vec![0, 0, 1, 0, 0, 0, 0, 1, 0]),
                    I1st(vec![0, 0, 1, 0, 0, 0, 0, 0, 1]),
                    I1st(vec![0, 0, 0, 1, 1, 0, 0, 0, 0]),
                    I1st(vec![0, 0, 0, 1, 0, 1, 0, 0, 0]),
                    I1st(vec![0, 0, 0, 1, 0, 0, 1, 0, 0]),
                    I1st(vec![0, 0, 0, 1, 0, 0, 0, 1, 0]),
                    I1st(vec![0, 0, 0, 1, 0, 0, 0, 0, 1]),
                    I1st(vec![0, 0, 0, 0, 1, 1, 0, 0, 0]),
                    I1st(vec![0, 0, 0, 0, 1, 0, 1, 0, 0]),
                    I1st(vec![0, 0, 0, 0, 1, 0, 0, 1, 0]),
                    I1st(vec![0, 0, 0, 0, 1, 0, 0, 0, 1]),
                    I1st(vec![0, 0, 0, 0, 0, 1, 1, 0, 0]),
                    I1st(vec![0, 0, 0, 0, 0, 1, 0, 1, 0]),
                    I1st(vec![0, 0, 0, 0, 0, 1, 0, 0, 1]),
                    I1st(vec![0, 0, 0, 0, 0, 0, 1, 1, 0]),
                    I1st(vec![0, 0, 0, 0, 0, 0, 1, 0, 1]),
                    I1st(vec![0, 0, 0, 0, 0, 0, 0, 1, 1]),
                ],
                modes: vec![
                    I1(0),
                    I1(1),
                    I1(2),
                    I1(3),
                    I1(4),
                    I1(5),
                    I1(6),
                    I1(7),
                    I1(8),
                ],
                ifunda: vec![1, 2, 3, 4, 5, 6, 7, 8, 9],
                iovrtn: vec![10, 11, 12, 13, 14, 15, 16, 17, 18],
                icombn: vec![
                    0, 19, 0, 20, 27, 0, 21, 28, 34, 0, 22, 29, 35, 40, 0, 23,
                    30, 36, 41, 45, 0, 24, 31, 37, 42, 46, 49, 0, 25, 32, 38,
                    43, 47, 50, 52, 0, 26, 33, 39, 44, 48, 51, 53, 54, 0,
                ],
            },
        ),
    ];
    inner(&tests);
}

#[test]
fn sym() {
    use mode::Mode::*;
    use resonance::Axis::*;
    use state::State::*;
    // NOTE just pasted the states in for now
    let tests = [
        Test::new(
            "nh3",
            Restst {
                coriolis: vec![],
                fermi1: vec![Fermi1::new(3, 2)],
                fermi2: vec![],
                darling: vec![],
                states: vec![
                    // ground state
                    I1st(vec![0, 0, 0, 0, 0, 0]),
                    // non-deg funds
                    I1st(vec![1, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 1, 0, 0, 0, 0]),
                    // deg funds
                    I2st(vec![(1, 1), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0)]),
                    I2st(vec![(0, 0), (1, 1), (0, 0), (0, 0), (0, 0), (0, 0)]),
                    // non-deg overtones
                    I1st(vec![2, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 2, 0, 0, 0, 0]),
                    // deg overtones
                    I2st(vec![(2, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0)]),
                    I2st(vec![(0, 0), (2, 0), (0, 0), (0, 0), (0, 0), (0, 0)]),
                    // nondeg-nondeg combination
                    I1st(vec![1, 1, 0, 0, 0, 0]),
                    // nondeg-deg combinations
                    I12st {
                        i1st: vec![1, 0, 0, 0, 0, 0],
                        i2st: vec![
                            (1, 1),
                            (0, 0),
                            (0, 0),
                            (0, 0),
                            (0, 0),
                            (0, 0),
                        ],
                    },
                    I12st {
                        i1st: vec![1, 0, 0, 0, 0, 0],
                        i2st: vec![
                            (0, 0),
                            (1, 1),
                            (0, 0),
                            (0, 0),
                            (0, 0),
                            (0, 0),
                        ],
                    },
                    I12st {
                        i1st: vec![0, 1, 0, 0, 0, 0],
                        i2st: vec![
                            (1, 1),
                            (0, 0),
                            (0, 0),
                            (0, 0),
                            (0, 0),
                            (0, 0),
                        ],
                    },
                    I12st {
                        i1st: vec![0, 1, 0, 0, 0, 0],
                        i2st: vec![
                            (0, 0),
                            (1, 1),
                            (0, 0),
                            (0, 0),
                            (0, 0),
                            (0, 0),
                        ],
                    },
                    // deg-deg combination
                    I2st(vec![(1, 1), (1, 1), (0, 0), (0, 0), (0, 0), (0, 0)]),
                ],
                modes: vec![I2(0, 1), I2(3, 4), I1(2), I1(5)],
                ifunda: vec![3, 3, 1, 4, 4, 2],
                iovrtn: vec![7, 7, 5, 8, 8, 6],
                icombn: vec![
                    0, 0, 0, 10, 10, 0, 14, 14, 11, 0, 14, 14, 11, 0, 0, 12,
                    12, 9, 13, 13, 0,
                ],
            },
        ),
        Test::new(
            "ph3",
            Restst {
                coriolis: vec![Coriolis::new(5, 3, B), Coriolis::new(3, 5, B)],
                fermi1: vec![Fermi1::new(3, 2), Fermi1::new(3, 0)],
                // you get it twice because we don't know about degeneracies at
                // this point
                fermi2: vec![Fermi2::new(3, 0, 3), Fermi2::new(3, 0, 3)],
                darling: vec![Darling::new(5, 3)],
                states: vec![
                    // ground state
                    I1st(vec![0, 0, 0, 0, 0, 0]),
                    // non-deg funds
                    I1st(vec![1, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 1, 0, 0, 0, 0]),
                    // deg funds
                    I2st(vec![(1, 1), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0)]),
                    I2st(vec![(0, 0), (1, 1), (0, 0), (0, 0), (0, 0), (0, 0)]),
                    // non-deg overtones
                    I1st(vec![2, 0, 0, 0, 0, 0]),
                    I1st(vec![0, 2, 0, 0, 0, 0]),
                    // deg overtones
                    I2st(vec![(2, 0), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0)]),
                    I2st(vec![(0, 0), (2, 0), (0, 0), (0, 0), (0, 0), (0, 0)]),
                    // nondeg-nondeg combination
                    I1st(vec![1, 1, 0, 0, 0, 0]),
                    // nondeg-deg combinations
                    I12st {
                        i1st: vec![1, 0, 0, 0, 0, 0],
                        i2st: vec![
                            (1, 1),
                            (0, 0),
                            (0, 0),
                            (0, 0),
                            (0, 0),
                            (0, 0),
                        ],
                    },
                    I12st {
                        i1st: vec![1, 0, 0, 0, 0, 0],
                        i2st: vec![
                            (0, 0),
                            (1, 1),
                            (0, 0),
                            (0, 0),
                            (0, 0),
                            (0, 0),
                        ],
                    },
                    I12st {
                        i1st: vec![0, 1, 0, 0, 0, 0],
                        i2st: vec![
                            (1, 1),
                            (0, 0),
                            (0, 0),
                            (0, 0),
                            (0, 0),
                            (0, 0),
                        ],
                    },
                    I12st {
                        i1st: vec![0, 1, 0, 0, 0, 0],
                        i2st: vec![
                            (0, 0),
                            (1, 1),
                            (0, 0),
                            (0, 0),
                            (0, 0),
                            (0, 0),
                        ],
                    },
                    // deg-deg combination
                    I2st(vec![(1, 1), (1, 1), (0, 0), (0, 0), (0, 0), (0, 0)]),
                ],
                modes: vec![I2(0, 1), I2(3, 4), I1(2), I1(5)],
                ifunda: vec![3, 3, 1, 4, 4, 2],
                iovrtn: vec![7, 7, 5, 8, 8, 6],
                icombn: vec![
                    0, 0, 0, 10, 10, 0, 14, 14, 11, 0, 14, 14, 11, 0, 0, 12,
                    12, 9, 13, 13, 0,
                ],
            },
        ),
        Test::new(
            "bipy",
            Restst {
                coriolis: vec![
                    Coriolis::new(8, 11, B),
                    Coriolis::new(9, 11, B),
                    Coriolis::new(11, 8, B),
                    Coriolis::new(11, 9, B),
                ],
                fermi1: vec![
                    Fermi1::new(13, 8),
                    Fermi1::new(13, 4),
                    Fermi1::new(13, 9),
                ],
                fermi2: vec![
                    Fermi2::new(2, 9, 13),
                    Fermi2::new(2, 13, 9),
                    Fermi2::new(8, 13, 13),
                    Fermi2::new(13, 4, 13),
                    Fermi2::new(13, 9, 2),
                    Fermi2::new(13, 9, 13),
                    Fermi2::new(13, 4, 13),
                    Fermi2::new(13, 9, 2),
                    Fermi2::new(13, 9, 13),
                ],
                darling: vec![
                    Darling::new(1, 0),
                    Darling::new(3, 2),
                    Darling::new(8, 9),
                    Darling::new(8, 11),
                    Darling::new(6, 4),
                    Darling::new(11, 9),
                ],
                states: include!("../../testfiles/bipy/states"),
                // I just copied these from the output
                modes: vec![
                    I2(4, 5),
                    I2(6, 7),
                    I2(9, 10),
                    I2(11, 12),
                    I2(13, 14),
                    I1(0),
                    I1(1),
                    I1(2),
                    I1(3),
                    I1(8),
                ],
                ifunda: vec![1, 2, 3, 4, 6, 6, 7, 7, 5, 8, 8, 9, 9, 10, 10],
                iovrtn: vec![
                    11, 12, 13, 14, 16, 16, 17, 17, 15, 18, 18, 19, 19, 20, 20,
                ],
                icombn: vec![
                    0, 21, 0, 22, 30, 0, 23, 31, 38, 0, 25, 33, 40, 46, 0, 25,
                    33, 40, 46, 0, 0, 26, 34, 41, 47, 56, 56, 0, 26, 34, 41,
                    47, 56, 56, 0, 0, 24, 32, 39, 45, 51, 51, 52, 52, 0, 27,
                    35, 42, 48, 57, 57, 60, 60, 53, 0, 27, 35, 42, 48, 57, 57,
                    60, 60, 53, 0, 0, 28, 36, 43, 49, 58, 58, 61, 61, 54, 63,
                    63, 0, 28, 36, 43, 49, 58, 58, 61, 61, 54, 63, 63, 0, 0,
                    29, 37, 44, 50, 59, 59, 62, 62, 55, 64, 64, 65, 65, 0, 29,
                    37, 44, 50, 59, 59, 62, 62, 55, 64, 64, 65, 65, 0, 0,
                ],
            },
        ),
        Test::new(
            "hmgnc",
            Restst {
                coriolis: vec![],
                fermi1: vec![Fermi1::new(3, 2)],
                fermi2: vec![
                    Fermi2::new(2, 3, 3),
                    Fermi2::new(2, 3, 5),
                    Fermi2::new(2, 5, 3),
                    Fermi2::new(5, 3, 2),
                    Fermi2::new(5, 3, 2),
                ],
                darling: vec![],
                states: include!("../../testfiles/hmgnc/states.rs"),
                // I just copied these from the output
                modes: vec![I2(3, 4), I2(5, 6), I1(0), I1(1), I1(2)],
                ifunda: vec![1, 2, 3, 4, 4, 5, 5],
                iovrtn: vec![6, 7, 8, 9, 9, 10, 10],
                icombn: vec![
                    0, 11, 0, 12, 15, 0, 13, 16, 18, 0, 13, 16, 18, 0, 0, 14,
                    17, 19, 20, 20, 0, 14, 17, 19, 20, 20, 0, 0,
                ],
            },
        ),
    ];
    inner(&tests[..]);
}
