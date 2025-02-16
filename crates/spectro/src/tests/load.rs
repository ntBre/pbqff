use std::str::FromStr;

use crate::dummy::DummyVal;
use crate::Curvil::*;
use crate::*;

use approx::assert_abs_diff_eq;
use na::matrix;
use nalgebra as na;
use symm::{Plane, PointGroup};

#[test]
fn load() {
    let got = Spectro::load("testfiles/spectro.in");
    let want = Spectro {
        header: vec![
            1, 1, 5, 2, 0, 0, 2, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
        ],
        geom: Molecule::from_str(
            "
C      0.00000000 -0.89784195  0.00000000
C      0.67631609  0.36939872  0.00000000
C     -0.67631609  0.36939872  0.00000000
H      1.55759065  0.94685790  0.00000000
H     -1.55759065  0.94685790  0.00000000
",
        )
        .unwrap(),
        weights: vec![
            (1, 12.0),
            (2, 12.0),
            (3, 12.0),
            (4, 1.007825),
            (5, 1.007825),
        ],
        curvils: vec![
            Bond(2, 3),
            Bond(1, 2),
            Bond(1, 3),
            Bond(2, 4),
            Bond(3, 5),
            Bend(2, 4, 1),
            Bend(3, 5, 1),
            Tors(4, 2, 1, 3),
            Tors(5, 3, 1, 2),
        ],
        degmodes: vec![],
        dummies: vec![],
        rotor: Rotor::AsymmTop,
        n3n: 15,
        i3n3n: 680,
        i4n3n: 3060,
        nvib: 9,
        i2vib: 45,
        i3vib: 165,
        i4vib: 495,
        natom: 5,
        axes: matrix![
        1.0, 0.0, 0.0;
        0.0, 1.0, 0.0;
        0.0, 0.0, 1.0;
        ],
        rotcon: vec![
            1.1424655406001178,
            1.0623779089716092,
            0.5504835966157277,
        ],
        primat: vec![14.755482612605023, 15.867828460536952, 30.62331107314198],
        iatom: 0,
        axis_order: 0,
        axis: Axis::Z,
        verbose: false,
        dump_fcs: false,
    };
    assert_eq!(got.header, want.header);
    assert_eq!(got.geom, want.geom);
    assert_eq!(got.weights, want.weights);
    assert_eq!(got.curvils, want.curvils);
    assert_eq!(got.degmodes, want.degmodes);
    assert_eq!(got.dummies, want.dummies);
    assert_eq!(got.rotor, want.rotor);
    assert_eq!(got.n3n, want.n3n);
    assert_eq!(got.i3n3n, want.i3n3n);
    assert_eq!(got.i4n3n, want.i4n3n);
    assert_eq!(got.nvib, want.nvib);
    assert_eq!(got.i2vib, want.i2vib);
    assert_eq!(got.i3vib, want.i3vib);
    assert_eq!(got.i4vib, want.i4vib);
    assert_eq!(got.natom, want.natom);
    assert_eq!(got.axes, want.axes);
    assert_eq!(got.primat, want.primat);
    assert_eq!(got.rotcon, want.rotcon);
    assert_eq!(got, want);
}

#[test]
fn load_dummy() {
    let got = Spectro::load("testfiles/dummy.in");
    let want = Spectro {
        header: vec![
            1, 1, 3, 0, 2, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 0,
        ],
        geom: Molecule::from_str(
            "
He    0.00000000  0.00000000  -1.76798207
H     0.00000000  0.00000000  -0.53078548
H     0.00000000  0.00000000   0.53078548
He    0.00000000  0.00000000   1.76798207
",
        )
        .unwrap(),
        weights: vec![
            (1, 4.00260325413),
            (2, 1.007825),
            (3, 1.007825),
            (4, 4.00260325413),
            (5, 0.0),
            (6, 0.0),
            (7, 0.0),
            (8, 0.0),
        ],
        curvils: vec![
            Bond(1, 2),
            Bond(2, 3),
            Bond(3, 4),
            Bend(2, 1, 3),
            Bend(3, 2, 4),
        ],
        degmodes: vec![vec![3, 2, 0], vec![1, 2, 3], vec![4, 6], vec![5, 7]],
        dummies: vec![
            Dummy {
                x: DummyVal::Value(1.1111111111),
                y: DummyVal::Atom(3),
                z: DummyVal::Atom(1),
            },
            Dummy {
                x: DummyVal::Atom(3),
                y: DummyVal::Value(1.1111111111),
                z: DummyVal::Atom(1),
            },
            Dummy {
                x: DummyVal::Value(1.1111111111),
                y: DummyVal::Atom(3),
                z: DummyVal::Atom(2),
            },
            Dummy {
                y: DummyVal::Value(1.1111111111),
                x: DummyVal::Atom(3),
                z: DummyVal::Atom(2),
            },
        ],
        rotor: Rotor::Linear,
        n3n: 12,
        i3n3n: 364,
        i4n3n: 1365,
        nvib: 7,
        i2vib: 28,
        i3vib: 84,
        i4vib: 210,
        natom: 4,
        axes: matrix![
        1.0, 0.0, 0.0;
        0.0, 1.0, 0.0;
        0.0, 0.0, -1.0;
        ],
        // TODO were these old values right or are the new ones?
        // primat: vec![0.0, 19.422657990598942, 19.422657990598942],
        // rotcon: vec![
        //     std::f64::INFINITY,
        //     0.8679363261189563,
        //     0.8679363261189563,
        // ],
        primat: vec![25.590234721060718, 25.590234721060718, 0.0],
        rotcon: vec![0.6587524734953539, 0.6587524734953539, 0.0],
        iatom: 0,
        axis_order: 0,
        axis: Axis::Z,
        verbose: false,
        dump_fcs: false,
    };
    assert_eq!(got.header, want.header);
    assert_eq!(got.geom, want.geom);
    assert_eq!(got.weights, want.weights);
    assert_eq!(got.curvils, want.curvils);
    assert_eq!(got.degmodes, want.degmodes);
    assert_eq!(got.dummies, want.dummies);
    assert_eq!(got.rotor, want.rotor);
    assert_eq!(got.n3n, want.n3n);
    assert_eq!(got.i3n3n, want.i3n3n);
    assert_eq!(got.i4n3n, want.i4n3n);
    assert_eq!(got.nvib, want.nvib);
    assert_eq!(got.i2vib, want.i2vib);
    assert_eq!(got.i3vib, want.i3vib);
    assert_eq!(got.i4vib, want.i4vib);
    assert_eq!(got.natom, want.natom);
    assert_eq!(got.axes, want.axes);
    assert_eq!(got.primat, want.primat);
    assert_eq!(got.rotcon, want.rotcon);
    assert_eq!(got, want);
}

struct Test {
    label: String,
    geom: Molecule,
    want_geom: Molecule,
    want_axes: Mat3,
}

impl Test {
    fn new(label: &str, geom: &str, want_geom: &str, want_axes: Mat3) -> Self {
        Self {
            label: String::from(label),
            geom: Molecule::from_str(geom).unwrap(),
            want_geom: Molecule::from_str(want_geom).unwrap(),
            want_axes,
        }
    }
}

/// test the geometry and axes produced by process_geom. output from
/// spectro2.out with input from spectro2.in
#[test]
fn geom_axes() {
    let tests = [
        // asymmetric tops
        Test::new(
            "h2o_sic",
            "
         H      0.0000000000      1.4313273344      0.9860352735
         O      0.0000000000      0.0000000000     -0.1242266321
         H      0.0000000000     -1.4313273344      0.9860352735
        ",
            "
        H     -0.7574256      0.5217723      0.0000000
        O      0.0000000     -0.0657528      0.0000000
        H      0.7574256      0.5217723      0.0000000
        ",
            matrix![
            0.00000000,0.00000000,1.00000000;
             1.00000000,0.00000000,0.00000000;
            0.00000000,1.00000000,0.00000000;
                            ],
        ),
        Test::new(
            "c3hcn",
            "
 H      0.0000000000      3.0019203461      3.8440336302
 C      0.0000000000      1.2175350291      2.8648781120
 C      0.0000000000     -1.4653360811      2.9704535522
 C      0.0000000000     -0.0243525793      0.6850185070
 C      0.0000000000     -0.0005362006     -1.9780266119
 N      0.0000000000      0.0178435988     -4.1716979516
",
            "
H     -2.03454900  1.58841460  0.00000000
C     -1.51635242  0.64418623  0.00000000
C     -1.57214543 -0.77553057  0.00000000
C     -0.36278597 -0.01293119  0.00000000
C      1.04643575 -0.00025357  0.00000000
N      2.20727579  0.00953399  0.00000000
",
            matrix![
            0.00000000,0.00000000,1.00000000;
            0.00005289,-1.00000000,0.00000000;
            1.00000000,0.00005289,0.00000000;
                ],
        ),
        Test::new(
            "c3hf",
            "
H   2.898322568200000   0.000000000000000  -3.153478939100000
C  -1.507649850300000   0.000000000000000  -2.021622563100000
C   1.224147945300000   0.000000000000000  -1.996566771000000
C  -0.024129098700000   0.000000000000000   0.153971073800000
F   0.040793403500000   0.000000000000000   2.610122623400000
        ",
            "
H      1.66901488  1.53379921  0.00000000
C      1.07022358 -0.79778171  0.00000000
C      1.05686473  0.64782210  0.00000000
C     -0.08110492 -0.01281612 -0.00000000
F     -1.38084634  0.02144956 -0.00000000
        ",
            matrix![
            0.00006911, 1.00000000, 0.00000000;
            0.00000000, 0.00000000, -1.00000000;
            1.00000000, -0.00006911, 0.00000000;
            ],
        ),
        // symmetric tops, index 3
        Test::new(
            "ph3",
            "
 H      2.2425420887      1.3150842078     -0.0000000000
 P      0.0000186413     -0.1314457979     -0.0000000000
 H     -1.1208854075      1.3156385211      1.9418758180
 H     -1.1208854075      1.3156385211     -1.9418758180
",
            "
H     1.1865658      0.0000000      0.6975731
P     0.0000000      0.0000000     -0.0680930
H    -0.5932829      1.0275961      0.6975731
H    -0.5932829     -1.0275961      0.6975731
",
            matrix![
            0.99999999,0.00000000, 0.00016482;
            -0.00016482,0.00000000, 0.99999999;
            0.00000000,1.00000000,0.00000000;
                            ],
        ),
        Test::new(
            "bipy",
            "
 C     -0.0000000000      0.0000000000      1.9792009664
 C      1.9835092315      0.0000000000      0.0000000000
 C     -0.9917546159     -1.7177693832      0.0000000000
 C     -0.9917546159      1.7177693832      0.0000000000
 C     -0.0000000000      0.0000000000     -1.9792009664
 H     -0.0000000000      0.0000000000      4.0059148040
 H     -0.0000000000      0.0000000000     -4.0059148040
",
            "
 C      0.0000000000      0.0000000000     -1.0473477000
 C      1.0496276000      0.0000000000      0.0000000000
 C     -0.5248138000     -0.9090042000      0.0000000000
 C     -0.5248138000      0.9090042000      0.0000000000
 C      0.0000000000      0.0000000000      1.0473477000
 H      0.0000000000      0.0000000000     -2.1198382000
 H      0.0000000000      0.0000000000      2.1198382000
",
            matrix![
            1.00000000,0.00000000,0.00000000;
            0.00000000,1.00000000,0.00000000;
            0.00000000,0.00000000,-1.00000000;
                                        ],
        ),
    ];
    for test in &tests[..] {
        let mut s = Spectro {
            geom: test.geom.clone(),
            ..Spectro::default()
        };
        process_geom(&mut s);
        check_eigen!(&s.axes, &test.want_axes, 1e-8, &test.label, &test.label);
        assert_abs_diff_eq!(s.geom, test.want_geom, epsilon = 1e-6);
    }
}

#[test]
fn test_make_iatl() {
    let mut s = Spectro {
        geom: Molecule::from_str(
            "
C      1.20584904 -0.69619727  0.00000000
C      0.00000000 -1.39239454  0.00000000
C     -1.20584904 -0.69619727  0.00000000
C     -1.20584904  0.69619727  0.00000000
C      0.00000000  1.39239454  0.00000000
C      1.20584904  0.69619727  0.00000000
H      0.00000000 -2.47385363  0.00000000
H     -2.14242009 -1.23692681  0.00000000
H     -2.14242009  1.23692681  0.00000000
H      0.00000000  2.47385363  0.00000000
H      2.14242009  1.23692681  0.00000000
H      2.14242009 -1.23692681  0.00000000
",
        )
        .unwrap(),
        ..Spectro::default()
    };

    s.make_iatl(PointGroup::D2h {
        axes: [Axis::Y, Axis::X, Axis::Z],
        planes: [
            Plane(Axis::X, Axis::Z),
            Plane(Axis::Y, Axis::Z),
            Plane(Axis::X, Axis::Y),
        ],
    });
}
