use std::io::BufRead;
use std::io::BufReader;

use approx::assert_abs_diff_eq;
use nalgebra as na;

use crate::fc::Fc;
use crate::Anpass;
use crate::Bias;
use crate::StatKind;

type Dmat = na::DMatrix<f64>;
type Dvec = na::DVector<f64>;

#[test]
fn test_load() {
    let anpass = Anpass::load_file("testfiles/anpass.in");
    let want = Anpass {
        #[rustfmt::skip]
            disps: Dmat::from_row_slice(
                69,
                3,
                &vec![
                    -0.00500000, -0.00500000, -0.01000000,
                    -0.00500000, -0.00500000, 0.00000000,
                    -0.00500000, -0.00500000, 0.01000000,
                    -0.00500000, -0.01000000, 0.00000000,
                    -0.00500000, -0.01500000, 0.00000000,
                    -0.00500000, 0.00000000, -0.01000000,
                    -0.00500000, 0.00000000, 0.00000000,
                    -0.00500000, 0.00000000, 0.01000000,
                    -0.00500000, 0.00500000, -0.01000000,
                    -0.00500000, 0.00500000, 0.00000000,
                    -0.00500000, 0.00500000, 0.01000000,
                    -0.00500000, 0.01000000, 0.00000000,
                    -0.00500000, 0.01500000, 0.00000000,
                    -0.01000000, -0.00500000, 0.00000000,
                    -0.01000000, -0.01000000, 0.00000000,
                    -0.01000000, 0.00000000, -0.01000000,
                    -0.01000000, 0.00000000, 0.00000000,
                    -0.01000000, 0.00000000, 0.01000000,
                    -0.01000000, 0.00500000, 0.00000000,
                    -0.01000000, 0.01000000, 0.00000000,
                    -0.01500000, -0.00500000, 0.00000000,
                    -0.01500000, 0.00000000, 0.00000000,
                    -0.01500000, 0.00500000, 0.00000000,
                    -0.02000000, 0.00000000, 0.00000000,
                    0.00000000, -0.00500000, -0.01000000,
                    0.00000000, -0.00500000, 0.00000000,
                    0.00000000, -0.00500000, 0.01000000,
                    0.00000000, -0.01000000, -0.01000000,
                    0.00000000, -0.01000000, 0.00000000,
                    0.00000000, -0.01000000, 0.01000000,
                    0.00000000, -0.01500000, 0.00000000,
                    0.00000000, -0.02000000, 0.00000000,
                    0.00000000, 0.00000000, -0.01000000,
                    0.00000000, 0.00000000, -0.02000000,
                    0.00000000, 0.00000000, 0.00000000,
                    0.00000000, 0.00000000, 0.01000000,
                    0.00000000, 0.00000000, 0.02000000,
                    0.00000000, 0.00500000, -0.01000000,
                    0.00000000, 0.00500000, 0.00000000,
                    0.00000000, 0.00500000, 0.01000000,
                    0.00000000, 0.01000000, -0.01000000,
                    0.00000000, 0.01000000, 0.00000000,
                    0.00000000, 0.01000000, 0.01000000,
                    0.00000000, 0.01500000, 0.00000000,
                    0.00000000, 0.02000000, 0.00000000,
                    0.00500000, -0.00500000, -0.01000000,
                    0.00500000, -0.00500000, 0.00000000,
                    0.00500000, -0.00500000, 0.01000000,
                    0.00500000, -0.01000000, 0.00000000,
                    0.00500000, -0.01500000, 0.00000000,
                    0.00500000, 0.00000000, -0.01000000,
                    0.00500000, 0.00000000, 0.00000000,
                    0.00500000, 0.00000000, 0.01000000,
                    0.00500000, 0.00500000, -0.01000000,
                    0.00500000, 0.00500000, 0.00000000,
                    0.00500000, 0.00500000, 0.01000000,
                    0.00500000, 0.01000000, 0.00000000,
                    0.00500000, 0.01500000, 0.00000000,
                    0.01000000, -0.00500000, 0.00000000,
                    0.01000000, -0.01000000, 0.00000000,
                    0.01000000, 0.00000000, -0.01000000,
                    0.01000000, 0.00000000, 0.00000000,
                    0.01000000, 0.00000000, 0.01000000,
                    0.01000000, 0.00500000, 0.00000000,
                    0.01000000, 0.01000000, 0.00000000,
                    0.01500000, -0.00500000, 0.00000000,
                    0.01500000, 0.00000000, 0.00000000,
                    0.01500000, 0.00500000, 0.00000000,
                    0.02000000, 0.00000000, 0.00000000,
                ],
            ),
        #[rustfmt::skip]
            energies: Dvec::from(vec![
                0.000128387078, 0.000027809414, 0.000128387078,
                0.000035977201, 0.000048243883, 0.000124321064,
                0.000023720402, 0.000124321065, 0.000124313373,
                0.000023689948, 0.000124313373, 0.000027697745,
                0.000035723392, 0.000102791171, 0.000113093098,
                0.000199639109, 0.000096581025, 0.000199639109,
                0.000094442297, 0.000096354531, 0.000228163468,
                0.000219814727, 0.000215550318, 0.000394681651,
                0.000100159437, 0.000001985383, 0.000100159437,
                0.000106187756, 0.000008036587, 0.000106187756,
                0.000018173585, 0.000032416257, 0.000098196697,
                0.000392997365, 0.000000000000, 0.000098196697,
                0.000392997364, 0.000100279477, 0.000002060371,
                0.000100279477, 0.000106387616, 0.000008146336,
                0.000106387616, 0.000018237641, 0.000032313930,
                0.000119935606, 0.000024112936, 0.000119935606,
                0.000028065156, 0.000036090120, 0.000120058596,
                0.000024213636, 0.000120058597, 0.000124214356,
                0.000028347337, 0.000124214356, 0.000036494030,
                0.000048633604, 0.000093011998, 0.000094882871,
                0.000188725453, 0.000095181193, 0.000188725453,
                0.000101370691, 0.000111560627, 0.000207527972,
                0.000211748039, 0.000219975758, 0.000372784451,
            ]),
        exponents: na::DMatrix::from_row_slice(
            3,
            22,
            &vec![
                0, 1, 0, 2, 1, 0, 0, 3, 2, 1, 0, 1, 0, 4, 3, 2, 1, 0, 2, 1, 0,
                0, 0, 0, 1, 0, 1, 2, 0, 0, 1, 2, 3, 0, 1, 0, 1, 2, 3, 4, 0, 1,
                2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 2,
                2, 2, 4,
            ],
        ),
        bias: None,
    };
    assert_abs_diff_eq!(anpass.disps, want.disps);
    assert_eq!(anpass.energies.len(), want.energies.len());
    assert_abs_diff_eq!(anpass.energies, want.energies);
    assert_eq!(anpass.exponents, want.exponents);
    assert_eq!(anpass.bias, want.bias);

    let got2 = Anpass::load_file("testfiles/anpass2.in");
    let want2 = Anpass {
        bias: Some(Bias {
            disp: na::dvector![
                -0.000045311426,
                -0.000027076533,
                0.000000000000
            ],
            energy: -0.000000002131,
        }),
        ..want
    };
    assert_abs_diff_eq!(got2.disps, want2.disps);
    assert_eq!(got2.energies.len(), want2.energies.len());
    assert_abs_diff_eq!(got2.energies, want2.energies);
    assert_eq!(got2.exponents, want2.exponents);
    assert_eq!(got2.bias, want2.bias);
}

#[test]
fn test_fit() {
    let anpass = Anpass::load_file("testfiles/anpass.in");
    let (got, _) = anpass.fit();
    let want = na::dvector![
        0.000000000002,
        0.000089167279,
        0.000008169225,
        0.958637185499,
        0.083537964315,
        0.080915238785,
        0.981791444232,
        -1.591452827847,
        -0.070076217196,
        -0.051299430149,
        -0.026818359356,
        -4.756703979348,
        0.045054940572,
        1.738582446151,
        -0.011114861654,
        0.021333149445,
        0.039547883446,
        -0.006274696681,
        10.448008327803,
        -0.141251401874,
        -0.047081672237,
        1.754904394438
    ];
    assert_abs_diff_eq!(got, want, epsilon = 1e-9);
}

#[test]
fn test_newton() {
    let anpass = Anpass::load_file("testfiles/c3h2.in");
    let (coeffs, _) = anpass.fit();
    let (got, kind) = anpass.newton(&coeffs).unwrap();
    let want = na::dvector![
        -0.000124209618,
        0.000083980449,
        -0.000036821098,
        -0.000117696241,
        0.000000000000,
        0.000000000000,
        0.000000000000,
        0.000000000000,
        0.000000000000,
        -0.000000022736
    ];
    assert_abs_diff_eq!(got, want, epsilon = 1e-12);
    assert_eq!(kind, StatKind::Min);
}

#[test]
fn test_eval() {
    let anpass = Anpass::load_file("testfiles/c3h2.in");
    let (coeffs, _) = anpass.fit();
    let (x, _) = anpass.newton(&coeffs).unwrap();
    let got = anpass.eval(&x, &coeffs);
    let want = -0.000000022736;
    assert!((got - want).abs() < 1e-12);
}

fn load9903(filename: &str) -> Vec<Fc> {
    let f = std::fs::File::open(filename).unwrap();
    let lines = BufReader::new(f).lines().map_while(Result::ok);
    let mut ret = Vec::new();
    for line in lines {
        ret.push(line.parse::<Fc>().unwrap());
    }
    ret
}

#[test]
fn test_bias() {
    let anpass = Anpass {
        disps: Dmat::from_row_slice(
            3,
            4,
            &[
                0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009,
                0.010, 0.011, 0.012,
            ],
        ),
        energies: na::dvector![10., 20., 30.],
        ..Anpass::load_file("testfiles/anpass.in")
    };
    let got = anpass.bias(&Bias {
        disp: na::dvector![0.001, 0.002, 0.003, 0.004],
        energy: 5.0,
    });
    let want_disps = Dmat::from_row_slice(
        3,
        4,
        &[
            0.000, 0.000, 0.000, 0.000, 0.004, 0.004, 0.004, 0.004, 0.008,
            0.008, 0.008, 0.008,
        ],
    );
    let want_energies = na::dvector![5., 15., 25.];
    assert_abs_diff_eq!(got.energies, want_energies);
    assert_abs_diff_eq!(got.disps, want_disps);
}

struct FullTest<'a> {
    infile: &'a str,
    want_file: &'a str,
    eps: f64,
}

impl<'a> FullTest<'a> {
    fn new(infile: &'a str, want_file: &'a str, eps: f64) -> Self {
        Self {
            infile,
            want_file,
            eps,
        }
    }
}

fn full_test(tests: &[FullTest]) {
    for test in tests {
        let now = std::time::Instant::now();
        println!("\nstarting {}", test.infile);

        let anpass = Anpass::load_file(test.infile);
        // initial fitting
        let (coeffs, _) = anpass.fit();
        // find stationary point
        let (x, _) = anpass.newton(&coeffs).unwrap();
        // determine energy at stationary point
        let e = anpass.eval(&x, &coeffs);
        // bias the displacements and energies to the new stationary point
        let anpass = anpass.bias(&Bias { disp: x, energy: e });
        // perform the refitting
        let (coeffs, _) = anpass.fit();
        let got = anpass.make9903(&coeffs);
        let want = load9903(test.want_file);
        assert_abs_diff_eq!(got[..], want, epsilon = test.eps);

        println!(
            "finishing after {:.3} sec",
            now.elapsed().as_millis() as f64 / 1000.
        );
    }
}

#[test]
fn test_full() {
    let tests = [
        FullTest::new("testfiles/h2o.in", "testfiles/h2o.9903", 9e-9),
        FullTest::new("testfiles/c3h2.in", "testfiles/c3h2.9903", 8.4e-8),
        FullTest::new("testfiles/c3h2_010.in", "testfiles/c3h2_010.9903", 7e-7),
        FullTest::new("testfiles/hoof.in", "testfiles/hoof.9903", 2e-7),
        FullTest::new("testfiles/hcf.in", "testfiles/hcf.9903", 2.8e-8),
    ];
    full_test(&tests);
}

#[test]
#[ignore]
fn test_full_long() {
    let tests = [
        FullTest::new("testfiles/c4h3.in", "testfiles/c4h3.9903", 3e-3),
        FullTest::new("testfiles/c5h2.in", "testfiles/c5h2.9903", 8e-4),
    ];
    full_test(&tests);
}
