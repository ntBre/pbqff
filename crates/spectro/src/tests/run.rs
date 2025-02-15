use std::{fs::read_to_string, path::PathBuf};

use approx::abs_diff_eq;
use serde::Deserialize;

use crate::*;

use na::DVector;
use nalgebra as na;

#[derive(Clone, Deserialize, Debug)]
struct Want {
    harms: Vec<f64>,
    funds: Vec<f64>,
    corrs: Vec<f64>,
    rots: Vec<Vec<f64>>,
}

fn load_want(filename: PathBuf, sym: bool) -> Output {
    let data = read_to_string(filename).unwrap();
    let want: Want = serde_json::from_str(&data).unwrap();
    let rots = if sym {
        rots_sym(&want)
    } else {
        rots_asym(&want)
    };
    Output {
        harms: want.harms,
        funds: want.funds,
        corrs: want.corrs,
        rots,
        // these are all untested here
        ..Output::default()
    }
}

fn rots_sym(want: &Want) -> Vec<Rot> {
    let mut rots = vec![Rot::new(
        State::I1st(vec![0; want.harms.len()]),
        want.rots[0][1],
        want.rots[0][0],
        want.rots[0][2],
    )];
    // TODO fix the states
    for i in 0..want.harms.len() {
        let mut tmp = vec![0; want.harms.len()];
        tmp[i] = 1;
        rots.push(Rot::new(
            State::I1st(tmp),
            want.rots[i + 1][1],
            want.rots[i + 1][0],
            want.rots[i + 1][2],
        ));
    }
    rots
}

fn rots_asym(want: &Want) -> Vec<Rot> {
    let mut rots = vec![Rot::new(
        State::I1st(vec![0; want.harms.len()]),
        want.rots[0][2],
        want.rots[0][0],
        want.rots[0][1],
    )];
    for i in 0..want.harms.len() {
        let mut tmp = vec![0; want.harms.len()];
        tmp[i] = 1;
        rots.push(Rot::new(
            State::I1st(tmp),
            want.rots[i + 1][2],
            want.rots[i + 1][0],
            want.rots[i + 1][1],
        ));
    }
    rots
}

#[derive(Clone)]
struct Test {
    infile: PathBuf,
    fort15: PathBuf,
    fort30: PathBuf,
    fort40: PathBuf,
    want: Output,
}

impl Test {
    fn new(dir: &'static str, sym: bool) -> Self {
        let start = Path::new("testfiles");
        Self {
            infile: start.join(dir).join("spectro.in"),
            fort15: start.join(dir).join("fort.15"),
            fort30: start.join(dir).join("fort.30"),
            fort40: start.join(dir).join("fort.40"),
            want: load_want(start.join(dir).join("summary.json"), sym),
        }
    }
}

fn check(got: Vec<f64>, want: Vec<f64>, msg: &str, label: &str) {
    assert_eq!(
        got.len(),
        want.len(),
        "got len={}, want len={} in {} at {}",
        got.len(),
        want.len(),
        msg,
        label
    );
    if !abs_diff_eq!(
        Dvec::from(got.clone()),
        Dvec::from(want.clone()),
        epsilon = 0.1
    ) {
        println!("\n{:>5}{:>8}{:>8}{:>8}", "Mode", "Got", "Want", "Diff");
        for i in 0..got.len() {
            println!(
                "{:5}{:8.1}{:8.1}{:8.1}",
                i + 1,
                got[i],
                want[i],
                got[i] - want[i],
            );
        }
        panic!("{msg} differ at {label}");
    }
}

macro_rules! check {
    ($got: expr, $want: expr, $msg: expr, $label: expr) => {
        check(
            $got,
            $want,
            $msg,
            &format!("'{}', {}:{}:{}", $label, file!(), line!(), column!(),),
        )
    };
}

fn check_rots(got: Vec<Rot>, want: Vec<Rot>, infile: &str) {
    assert_eq!(got.len(), want.len());
    if !abs_diff_eq!(
        DVector::from(got.clone()),
        DVector::from(want.clone()),
        epsilon = 2e-5
    ) {
        println!("got");
        for g in &got {
            println!("{g}");
        }
        println!("want");
        for g in &want {
            println!("{g}");
        }
        println!("diff");
        let mut max = 0.0;
        for g in 0..got.len() {
            let a = got[g].a - want[g].a;
            let b = got[g].b - want[g].b;
            let c = got[g].c - want[g].c;
            for i in [a, b, c] {
                if i.abs() > max {
                    max = i.abs();
                }
            }
            println!("{}", Rot::new(got[g].state.clone(), a, b, c));
        }
        panic!("rots differ by {max:.2e} at {infile}");
    }
}

macro_rules! inner {
    ($tests: expr) => {
        for test in Vec::from(&$tests[..]) {
            let infile = test.infile.to_str().unwrap();

            let spectro = Spectro::load(infile);
            let got = spectro.run_files(test.fort15, test.fort30, test.fort40);

            check!(got.harms, test.want.harms, "harms", infile);
            check!(got.funds, test.want.funds, "funds", infile);
            check!(got.corrs, test.want.corrs, "corrs", infile);
            check_rots(got.rots, test.want.rots, infile);
        }
    };
}

#[test]
fn diatomic() {
    let tests = [
        //
        Test::new("alh", false),
    ];
    inner!(&tests);
}

#[test]
fn asym() {
    let tests = [
        Test::new("h2o", false),
        Test::new("h2co", false),
        Test::new("c3h2", false),
        Test::new("c2h4", false),
        Test::new("c4h3+", false),
        Test::new("allyl", false),
        Test::new("h2o_sic", false),
        Test::new("c3hf", false),
        Test::new("c3hcl", false),
        Test::new("c3hcn", false),
        Test::new("c3hcn010", false),
        Test::new("c-hoco", false),
        Test::new("hno", false),
        Test::new("hocn", false),
        Test::new("hoco+", false),
        Test::new("hpsi", false),
        Test::new("hso", false),
        Test::new("hss", false),
        Test::new("nh2-", false),
        Test::new("nnoh+", false),
        Test::new("sic2", false),
        Test::new("t-hoco", false),
        Test::new("hoof", false),
        Test::new("hosh", false),
        Test::new("hssh", false),
        Test::new("nosym_tdrane", false),
        Test::new("h2sc", false),
        Test::new("d2sc", false),
    ];
    inner!(&tests);
}

#[test]
fn sym() {
    use State::*;
    // manual workaround for the states not being in the summarize output
    let states = [
        vec![
            I1st(vec![0, 0, 0, 0, 0, 0]),
            I1st(vec![1, 0, 0, 0, 0, 0]),
            I1st(vec![0, 1, 0, 0, 0, 0]),
            I2st(vec![(1, 1), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0)]),
            I2st(vec![(0, 0), (1, 1), (0, 0), (0, 0), (0, 0), (0, 0)]),
        ],
        vec![
            I1st(vec![0, 0, 0, 0, 0, 0]),
            I1st(vec![1, 0, 0, 0, 0, 0]),
            I1st(vec![0, 1, 0, 0, 0, 0]),
            I2st(vec![(1, 1), (0, 0), (0, 0), (0, 0), (0, 0), (0, 0)]),
            I2st(vec![(0, 0), (1, 1), (0, 0), (0, 0), (0, 0), (0, 0)]),
        ],
        vec![
            I1st(vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
            I1st(vec![1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
            I1st(vec![0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
            I1st(vec![0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
            I1st(vec![0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
            I1st(vec![0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]),
            I2st(vec![
                (1, 1),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
            ]),
            I2st(vec![
                (0, 0),
                (1, 1),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
            ]),
            I2st(vec![
                (0, 0),
                (0, 0),
                (1, 1),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
            ]),
            I2st(vec![
                (0, 0),
                (0, 0),
                (0, 0),
                (1, 1),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
            ]),
            I2st(vec![
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (1, 1),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
            ]),
            I2st(vec![
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (1, 1),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 0),
            ]),
        ],
        // include!("../../testfiles/c3h3+/states.rs"),
    ];
    let mut tests = [
        Test::new("nh3", true),
        Test::new("ph3", true),
        Test::new("bipy", true),
        // Test::new("c3h3+", true),
    ];
    for (i, test) in tests.iter_mut().enumerate() {
        for j in 0..test.want.rots.len() {
            test.want.rots[j].state = states[i][j].clone();
        }
    }
    inner!(&tests[..]);
}

#[test]
fn lin() {
    use State::*;
    let lin_tri = vec![
        I1st(vec![0, 0, 0, 0]),
        I1st(vec![1, 0, 0, 0]),
        I1st(vec![0, 1, 0, 0]),
        I2st(vec![(1, 1), (0, 0), (0, 0), (0, 0)]),
    ];
    let states = [
        // triatomics
        lin_tri.clone(),
        lin_tri.clone(),
        lin_tri.clone(),
        lin_tri.clone(),
        lin_tri,
        // tetraatomics
        include!("../../testfiles/hmgnc/states.rs"),
        include!("../../testfiles/hmgnc/states.rs"),
    ];
    let mut tests = [
        Test::new("c2h-", true),
        Test::new("hcn", true),
        Test::new("hco+", true),
        Test::new("hnc", true),
        Test::new("h2o+", true),
        Test::new("hmgnc", true),
        Test::new("alnco", true),
    ];
    for (i, test) in tests.iter_mut().enumerate() {
        for j in 0..test.want.rots.len() {
            test.want.rots[j].state = states[i][j].clone();
        }
    }
    inner!(&tests);
}

#[test]
#[ignore]
fn sphere() {
    use State::*;
    let states = [include!("../../testfiles/tdrane/states.rs")];
    // can't be detected as a spherical top here or in fortran version.
    // currently failing with "big problem in xcals" because spherical tops
    // should actually run through wcals. but I can't decide to use wcals
    // without detecting a true spherical top. in the fortran version it fails
    // at restst.f:1244 because the size of i3[xyz]st wasn't set properly at
    // mains.f:429 because isphtp isn't set from dist.f, so the root issue is
    // the same as I'm having. see if methane has the same issue, and then
    // potentially lower the threshold for the third moi. if two of them are the
    // same to the rigorous threshold and one differs by a slightly larger
    // amount, it might still be worth using the spherical top routines. in
    // fortran I can just set n3sz to nvib even for a symmetric top and then
    // trust the user's degmode input, but that should then hit the wcals vs
    // xcals issue I'm having here.

    // I can already tell ch4 is going to have the same issue. molpro optimized
    // the geometry in Cs symmetry and already shows that it has two different
    // rotational constants. How then do you ever get a true symmetric top? I
    // guess you must have to use a larger tolerance on the third rotational
    // constant
    let mut tests = [Test::new("tdrane", true)];
    for (i, test) in tests.iter_mut().enumerate() {
        for j in 0..test.want.rots.len() {
            test.want.rots[j].state = states[i][j].clone();
        }
    }
    inner!(&tests);
}
