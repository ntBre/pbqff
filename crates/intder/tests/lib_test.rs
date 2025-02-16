use std::f64::consts::PI;
use std::io::{BufRead, BufReader, Read};

use approx::{abs_diff_ne, assert_abs_diff_eq};

use intder::geom::*;
use intder::*;
use nalgebra as na;

const S: f64 = std::f64::consts::SQRT_2 / 2.;

#[test]
fn test_load_pts() {
    let got = Intder::load_file("testfiles/intder.in");
    let want = Intder {
        input_options: vec![3, 3, 3, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 14],
        simple_internals: vec![
            Siic::Stretch(0, 1),
            Siic::Stretch(1, 2),
            Siic::Bend(0, 1, 2),
        ],
        symmetry_internals: vec![
            vec![S, S, 0.],
            vec![0., 0., 1.],
            vec![S, -S, 0.],
        ],
        geom: Geom(vec![
            na::Vector3::new(0.000000000000, 1.431390244079, 0.986041163966),
            na::Vector3::new(0.000000000000, 0.000000000000, -0.124238450265),
            na::Vector3::new(0.000000000000, -1.431390244079, 0.986041163966),
        ]),
        disps: vec![
            vec![0.005, 0.0, 0.0],
            vec![0.0, 0.005, 0.0],
            vec![0.0, 0.0, 0.005],
            vec![-0.005, -0.005, -0.01],
            vec![-0.005, -0.005, 0.0],
            vec![-0.005, -0.005, 0.010],
            vec![-0.005, -0.010, 0.0],
            vec![-0.005, -0.015, 0.0],
            vec![0.0, 0.0, 0.0],
        ],
        atoms: vec![],
        fc2: vec![],
        fc3: vec![],
        fc4: vec![],
    };
    assert_eq!(got, want);
}

#[test]
fn load_full_xyz() {
    let got = Intder::load_file("testfiles/full_xyz.in");
    let want = Intder {
        input_options: vec![3, 3, 3, 0, 0, 3, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 14],
        simple_internals: vec![
            Siic::Stretch(0, 1),
            Siic::Stretch(1, 2),
            Siic::Bend(0, 1, 2),
        ],
        symmetry_internals: vec![
            vec![S, S, 0.],
            vec![0., 0., 1.],
            vec![S, -S, 0.],
        ],
        geom: Geom(vec![
            na::Vector3::new(0.000000000000, 1.431390244079, 0.986041163966),
            na::Vector3::new(0.000000000000, 0.000000000000, -0.124238450265),
            na::Vector3::new(0.000000000000, -1.431390244079, 0.986041163966),
        ]),
        disps: vec![
            vec![0.005, 0.0, 0.0],
            vec![0.0, 0.005, 0.0],
            vec![0.0, 0.0, 0.005],
            vec![-0.005, -0.005, -0.01],
            vec![-0.005, -0.005, 0.0],
            vec![-0.005, -0.005, 0.010],
            vec![-0.005, -0.010, 0.0],
            vec![-0.005, -0.015, 0.0],
            vec![0.0, 0.0, 0.0],
        ],
        atoms: vec![
            Atom {
                label: "H".to_string(),
                weight: 1,
            },
            Atom {
                label: "O".to_string(),
                weight: 16,
            },
            Atom {
                label: "H".to_string(),
                weight: 1,
            },
        ],
        fc2: vec![],
        fc3: vec![],
        fc4: vec![],
    };
    assert_eq!(got, want);
}

#[test]
fn test_load_freqs() {
    let got = Intder::load_file("testfiles/h2o.freq.in");
    let want = Intder {
        input_options: vec![3, 3, 3, 4, 0, 3, 2, 0, 0, 1, 3, 0, 0, 0, 0, 0, 14],
        simple_internals: vec![
            Siic::Stretch(0, 1),
            Siic::Stretch(1, 2),
            Siic::Bend(0, 1, 2),
        ],
        symmetry_internals: vec![
            vec![S, S, 0.],
            vec![0., 0., 1.],
            vec![S, -S, 0.],
        ],
        geom: Geom(vec![
            na::Vector3::new(0.0000000000, 1.4313273344, 0.9860352735),
            na::Vector3::new(0.0000000000, 0.0000000000, -0.1242266321),
            na::Vector3::new(0.0000000000, -1.4313273344, 0.9860352735),
        ]),
        disps: vec![],
        atoms: vec![
            Atom {
                label: "H".to_string(),
                weight: 1,
            },
            Atom {
                label: "O".to_string(),
                weight: 16,
            },
            Atom {
                label: "H".to_string(),
                weight: 1,
            },
        ],
        #[rustfmt::skip]
        fc2: vec![
	    8.360863692425, 0.364250381719, 0.0,
	    0.0, 0.705590041836, 0.0,
	    0.0, 0.0, 8.562725561924,
	],
        #[rustfmt::skip]
        fc3: vec![
            -41.638868371821, -0.611029974345, -0.447356783193,
            -0.701565547377, 0.0, 0.0,
            0.0, -41.484904978222, 0.392943158463,
        ],
        #[rustfmt::skip]
        fc4: vec![
            181.917347386091, -0.292134838642, 0.372522604963,
            1.034069182728, -0.6558318021, 0.0,
            0.0, 0.0, 0.0,
            182.20617803073, -1.233550191665, -0.820459303061,
            0.0, 0.0, 183.621273961541,
        ],
    };
    assert_eq!(got.input_options, want.input_options);
    assert_eq!(got.simple_internals, want.simple_internals);
    assert_eq!(got.symmetry_internals, want.symmetry_internals);
    assert_eq!(got.geom, want.geom);
    assert_eq!(got.disps, want.disps);
    assert_eq!(got.atoms, want.atoms);
    assert_eq!(got.fc2, want.fc2);
    assert_eq!(got.fc3, want.fc3);
    assert_eq!(got.fc4, want.fc4);
    assert_eq!(got, want);
}

#[test]
fn test_initial_values_simple() {
    let tests = vec![
        (
            "testfiles/intder.in",
            vec![0.9586143145, 0.9586143145, 1.8221415968],
        ),
        (
            "testfiles/c7h2.in",
            vec![
                1.4260535407,
                1.4260535407,
                1.3992766813,
                1.3992766813,
                2.6090029486,
                2.6090029486,
                3.6728481977,
                3.6728481977,
                2.5991099760,
                2.5991099760,
                2.5961248359,
                2.5961248359,
                2.5945738184,
                2.5945738184,
                1.0819561376,
                PI,
                PI,
                PI,
                PI,
                PI,
                PI,
            ],
        ),
    ];
    for test in tests {
        let intder = Intder::load_file(test.0);
        let got = intder.simple_values(&intder.geom);
        assert_abs_diff_eq!(
            DVec::from(got),
            DVec::from(test.1),
            epsilon = 3e-7
        );
    }
}

#[test]
fn test_initial_values_symmetry() {
    let intder = Intder::load_file("testfiles/intder.in");
    let got = intder.symmetry_values(&intder.geom);
    let got = got.as_slice();
    let want = vec![1.3556853647, 1.8221415968, 0.0000000000];
    let want = want.as_slice();
    assert_abs_diff_eq!(got, want, epsilon = 1e-7);
}

fn load_vec(filename: &str) -> Vec<f64> {
    let mut f = std::fs::File::open(filename).unwrap();
    let mut buf = String::new();
    f.read_to_string(&mut buf).unwrap();
    buf.split_whitespace()
        .map(|x| x.parse::<f64>().unwrap())
        .collect::<Vec<_>>()
}

struct MatTest<'a> {
    infile: &'a str,
    rows: usize,
    cols: usize,
    vecfile: &'a str,
    eps: f64,
}

#[test]
fn test_b_matrix() {
    let tests = vec![
        MatTest {
            infile: "testfiles/intder.in",
            rows: 3,
            cols: 9,
            vecfile: "testfiles/h2o.bmat",
            eps: 2e-7,
        },
        MatTest {
            infile: "testfiles/c7h2.in",
            rows: 21,
            cols: 27,
            vecfile: "testfiles/c7h2.bmat",
            eps: 2.2e-7,
        },
    ];
    for test in tests {
        let intder = Intder::load_file(test.infile);
        let want =
            DMat::from_row_slice(test.rows, test.cols, &load_vec(test.vecfile));
        let got = intder.b_matrix(&intder.geom);
        assert_abs_diff_eq!(got, want, epsilon = test.eps);
    }
}

#[test]
fn test_sym_b() {
    let tests = vec![
        MatTest {
            infile: "testfiles/intder.in",
            rows: 3,
            cols: 9,
            vecfile: "testfiles/h2o.bsmat",
            eps: 2e-7,
        },
        MatTest {
            infile: "testfiles/c7h2.in",
            rows: 21,
            cols: 27,
            vecfile: "testfiles/c7h2.bsmat",
            eps: 2.2e-7,
        },
    ];
    for test in tests {
        let intder = Intder::load_file(test.infile);
        let got = intder.sym_b_matrix(&intder.geom);
        let want =
            DMat::from_row_slice(test.rows, test.cols, &load_vec(test.vecfile));
        assert_abs_diff_eq!(got, want, epsilon = test.eps);
    }
}

#[allow(dead_code)]
fn dbg_mat(a: &DMat, b: &DMat, eps: f64) {
    let a = a.as_slice();
    let b = b.as_slice();
    assert!(a.len() == b.len());
    println!();
    for i in 0..a.len() {
        if (a[i] - b[i]).abs() > eps {
            println!(
                "{:5}{:>15.8}{:>15.8}{:>15.8e}",
                i,
                a[i],
                b[i],
                a[i] - b[i]
            );
        }
    }
}

#[test]
fn test_a_matrix() {
    let tests = vec![
        MatTest {
            infile: "testfiles/intder.in",
            rows: 9,
            cols: 3,
            vecfile: "testfiles/h2o.amat",
            eps: 3e-8,
        },
        // low precision from intder.out
        MatTest {
            infile: "testfiles/c3h2.in",
            rows: 15,
            cols: 9,
            vecfile: "testfiles/c3h2.amat",
            eps: 1e-6,
        },
    ];
    for test in tests {
        let intder = Intder::load_file(test.infile);
        let load = load_vec(test.vecfile);
        let want = DMat::from_row_slice(test.rows, test.cols, &load);
        let got = Intder::a_matrix(&intder.sym_b_matrix(&intder.geom));
        assert_abs_diff_eq!(got, want, epsilon = test.eps);
    }
}

/// load a file where each line is a DVec
fn load_geoms(filename: &str) -> Vec<DVec> {
    let f = std::fs::File::open(filename).unwrap();
    let lines = BufReader::new(f).lines().map_while(Result::ok);
    let mut ret = Vec::new();
    for line in lines {
        if !line.is_empty() {
            ret.push(DVec::from(
                line.split_whitespace()
                    .map(|x| x.parse().unwrap())
                    .collect::<Vec<_>>(),
            ));
        }
    }
    ret
}

#[test]
fn test_convert_disps() {
    struct Test<'a> {
        infile: &'a str,
        wantfile: &'a str,
    }
    let tests = vec![
        Test {
            infile: "testfiles/h2o.in",
            wantfile: "testfiles/h2o.small.07",
        },
        Test {
            infile: "testfiles/thoco.in",
            wantfile: "testfiles/thoco.07",
        },
        Test {
            infile: "testfiles/c3h2.in",
            wantfile: "testfiles/c3h2.07",
        },
        Test {
            infile: "testfiles/c7h2.in",
            wantfile: "testfiles/c7h2.small.07",
        },
        Test {
            infile: "testfiles/c2h4.in",
            wantfile: "testfiles/c2h4.07",
        },
        Test {
            infile: "testfiles/c2h2.in",
            wantfile: "testfiles/c2h2.07",
        },
        Test {
            infile: "testfiles/h2co.in",
            wantfile: "testfiles/h2co.07",
        },
        Test {
            infile: "testfiles/halnh.in",
            wantfile: "testfiles/halnh.07",
        },
    ];
    for test in tests {
        let intder = Intder::load_file(test.infile);
        let got = intder.convert_disps().unwrap();
        let want = load_geoms(test.wantfile);
        for i in 0..got.len() {
            assert_abs_diff_eq!(got[i], want[i], epsilon = 4e-8);
        }
    }
}

fn load_tensor(filename: &str) -> Tensor3 {
    let f = std::fs::File::open(filename).expect("failed to open tensor file");
    let lines = BufReader::new(f).lines();
    let mut hold = Vec::new();
    let mut buf = Vec::new();
    for line in lines.map_while(Result::ok) {
        let mut fields = line.split_whitespace().peekable();
        if fields.peek().is_none() {
            // in between chunks
            hold.push(buf);
            buf = Vec::new();
        } else {
            let row = fields
                .map(|s| s.parse::<f64>().unwrap())
                .collect::<Vec<_>>();
            buf.push(row);
        }
    }
    hold.push(buf);
    let a = hold.len();
    let b = hold[0].len();
    let c = hold[0][0].len();
    let mut ret = Tensor3::zeros((a, b, c));
    for i in 0..a {
        for j in 0..b {
            for k in 0..c {
                ret[(i, j, k)] = hold[i][j][k];
            }
        }
    }
    ret
}

#[test]
fn test_lintr_fc3() {
    let intder = Intder::load_file("testfiles/h2o.freq.in");
    let b_sym = intder.sym_b_matrix(&intder.geom);
    let got = intder.lintr_fc3(&b_sym);
    let want = load_tensor("testfiles/h2o.lintr.f3");
    assert_abs_diff_eq!(got, want, epsilon = 1e-6);
}

#[test]
fn test_convert_fcs() {
    struct Test {
        infile: &'static str,
        fcs: (String, String, String),
        sizes: (usize, usize, usize),
        eps: (f64, f64, f64),
    }
    let tests = vec![
        Test {
            infile: "testfiles/h2o.freq.in",
            fcs: (
                "testfiles/h2o.15".to_string(),
                "testfiles/h2o.30".to_string(),
                "testfiles/h2o.40".to_string(),
            ),
            sizes: (9, 165, 495),
            eps: (1e-7, 4e-7, 2e-6),
        },
        Test {
            infile: "testfiles/thoco.freq.in",
            fcs: (
                "testfiles/thoco.15".to_string(),
                "testfiles/thoco.30".to_string(),
                "testfiles/thoco.40".to_string(),
            ),
            sizes: (12, 364, 1365),
            eps: (2e-7, 6e-7, 3e-6),
        },
        Test {
            infile: "testfiles/c3h2.freq.in",
            fcs: (
                "testfiles/c3h2.15".to_string(),
                "testfiles/c3h2.30".to_string(),
                "testfiles/c3h2.40".to_string(),
            ),
            sizes: (15, 680, 3060),
            eps: (2e-7, 4e-7, 2e-6),
        },
        Test {
            infile: "testfiles/c2h2.freq.in",
            fcs: (
                "testfiles/c2h2.15".to_string(),
                "testfiles/c2h2.30".to_string(),
                "testfiles/c2h2.40".to_string(),
            ),
            sizes: (12, 364, 1365),
            eps: (3e-7, 7e-7, 4e-6),
        },
        Test {
            infile: "testfiles/h2co.freq.in",
            fcs: (
                "testfiles/h2co.15".to_string(),
                "testfiles/h2co.30".to_string(),
                "testfiles/h2co.40".to_string(),
            ),
            sizes: (12, 364, 1365),
            eps: (2e-7, 7e-7, 3e-6),
        },
        Test {
            infile: "testfiles/halnh.freq.in",
            fcs: (
                "testfiles/halnh.15".to_string(),
                "testfiles/halnh.30".to_string(),
                "testfiles/halnh.40".to_string(),
            ),
            sizes: (12, 364, 1365),
            eps: (2e-7, 7e-7, 3e-6),
        },
    ];
    for test in tests {
        let intder = Intder::load_file(test.infile);
        let (fc2, fc3, fc4) = intder.convert_fcs();

        let want_fc2 = DMat::from_row_slice(
            test.sizes.0,
            test.sizes.0,
            &load_vec(&test.fcs.0),
        );
        assert_abs_diff_eq!(fc2, want_fc2, epsilon = test.eps.0);

        let want_fc3 =
            DMat::from_row_slice(test.sizes.1, 1, &load_vec(&test.fcs.1));
        assert_abs_diff_eq!(
            DMat::from_row_slice(test.sizes.1, 1, &fc3),
            want_fc3,
            epsilon = test.eps.1
        );

        let want_fc4 =
            DMat::from_row_slice(test.sizes.2, 1, &load_vec(&test.fcs.2));
        let got_fc4 = DMat::from_row_slice(test.sizes.2, 1, &fc4);
        if abs_diff_ne!(got_fc4, want_fc4, epsilon = test.eps.2) {
            println!("max diff={:.2e}", (&got_fc4 - &want_fc4).abs().max());
            for i in 1..=test.sizes.0 {
                for j in 1..=i {
                    for k in 1..=j {
                        for l in 1..=k {
                            println!(
                                "{i:5}{j:5}{k:5}{l:5}{:20.12}{:20.12}{:20.12}",
                                got_fc4[fc4_index(i, j, k, l)],
                                want_fc4[fc4_index(i, j, k, l)],
                                got_fc4[fc4_index(i, j, k, l)]
                                    - want_fc4[fc4_index(i, j, k, l)]
                            );
                        }
                    }
                }
            }
            panic!("fc4 mismatch");
        }
    }
}
