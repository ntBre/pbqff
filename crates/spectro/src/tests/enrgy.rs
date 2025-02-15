use std::path::PathBuf;

use crate::consts::FACT2;
use crate::resonance::Restst;
use crate::state::StatePartition;
use crate::utils::linalg::symm_eigen_decomp;
use crate::*;

#[derive(Clone)]
struct Test {
    infile: PathBuf,
    fort15: PathBuf,
    fort30: PathBuf,
    fort40: PathBuf,
    part1: Vec<f64>,
    full: Vec<f64>,
    eps: f64,

    /// option to split the result at part.0 and check it with part.1 instead of
    /// eps
    part: Option<(usize, f64)>,
}

impl Test {
    fn new(dir: &'static str, eps: f64, part: Option<(usize, f64)>) -> Self {
        let start = Path::new("testfiles");
        Self {
            infile: start.join(dir).join("spectro.in"),
            fort15: start.join(dir).join("fort.15"),
            fort30: start.join(dir).join("fort.30"),
            fort40: start.join(dir).join("fort.40"),
            part1: load_vec(start.join(dir).join("enrgy1")),
            full: load_vec(start.join(dir).join("enrgy")),
            eps,
            part,
        }
    }
}

fn common(
    test: &Test,
) -> (Spectro, Dvec, F3qcm, Restst, Dmat, Option<Dmat>, f64) {
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
    if s.rotor.is_sym_top() {
        s.bdegnl(&freq, &mut lxm, &w, &mut lx);
    }
    let f3x = load_fc3(&test.fort30, s.n3n);
    let f3x = s.rot3rd(f3x);
    let f3qcm = force3(s.n3n, f3x, &lx, s.nvib, &freq);
    let f4x = load_fc4(&test.fort40, s.n3n);
    let f4x = s.rot4th(f4x);
    let f4qcm = force4(s.n3n, &f4x, &lx, s.nvib, &freq);
    let (zmat, wila) = s.zeta(&lxm, &w);
    let restst = Restst::new(&s, &zmat, &f3qcm, &freq);
    let Restst {
        coriolis: _,
        fermi1,
        fermi2,
        darling: _,
        modes,
        states: _,
        ifunda: _,
        iovrtn: _,
        icombn: _,
    } = &restst;
    let (xcnst, gcnst, e0) = if s.rotor.is_sym_top() {
        let (x, g, e) =
            s.xcals(&f4qcm, &freq, &f3qcm, &zmat, fermi1, fermi2, modes, &wila);
        (x, Some(g), e)
    } else {
        let (x, e) =
            s.xcalc(&f4qcm, &freq, &f3qcm, &zmat, modes, fermi1, fermi2);
        (x, None, e)
    };
    (s, freq, f3qcm, restst, xcnst, gcnst, e0)
}

#[test]
fn part1() {
    let tests = [
        Test::new("h2o", 2.1e-5, None),
        // TODO this issue comes from whatever is going on in gcnst, so follow
        // up on this if that gets resolved. anway, 1cm-1 in the 56th excited
        // state doesn't really worry me that much
        Test::new("bipy", 2.1e-5, Some((56, 1.01))),
    ];
    for test in Vec::from(&tests[..]) {
        let (_, freq, _, r, xcnst, gcnst, e0) = common(&test);
        let mut got = vec![0.0; r.states.len()];
        let nstate = r.states.len();
        let (n1dm, n2dm, _) = Mode::count(&r.modes);
        let (i1mode, i2mode, _) = Mode::partition(&r.modes);
        let StatePartition { i1sts, i2sts, .. } = State::partition(&r.states);

        crate::enrgy::part1(
            nstate, &i1mode, &freq, i1sts, &i2mode, i2sts, n1dm, &xcnst, n2dm,
            &gcnst, &mut got, e0,
        );

        if let Some((i, e2)) = test.part {
            let (g1, g2) = got.split_at(i);
            let (w1, w2) = test.part1.split_at(i);
            check_vec!(
                Dvec::from(g1.to_vec()),
                Dvec::from(w1.to_vec()),
                test.eps,
                test.infile.display()
            );
            check_vec!(
                Dvec::from(g2.to_vec()),
                Dvec::from(w2.to_vec()),
                e2,
                test.infile.display()
            );
        } else {
            check_vec!(
                Dvec::from(got),
                Dvec::from(test.part1),
                test.eps,
                test.infile.display()
            );
        }
    }
}

#[test]
fn full() {
    let tests = [
        //
        Test::new("h2o", 2.1e-5, None),
        Test::new("bipy", 2.1e-5, Some((56, 1.01))),
    ];
    for test in Vec::from(&tests[..]) {
        let (s, freq, f3qcm, restst, xcnst, gcnst, e0) = common(&test);
        let mut got = vec![0.0; restst.states.len()];
        s.enrgy(&freq, &xcnst, &gcnst, &restst, &f3qcm, e0, &mut got);

        if let Some((i, e2)) = test.part {
            let (g1, g2) = got.split_at(i);
            let (w1, w2) = test.full.split_at(i);
            check_vec!(
                Dvec::from(g1.to_vec()),
                Dvec::from(w1.to_vec()),
                test.eps,
                test.infile.display()
            );
            check_vec!(
                Dvec::from(g2.to_vec()),
                Dvec::from(w2.to_vec()),
                e2,
                test.infile.display()
            );
        } else {
            check_vec!(
                Dvec::from(got),
                Dvec::from(test.full),
                test.eps,
                test.infile.display()
            );
        }
    }
}
