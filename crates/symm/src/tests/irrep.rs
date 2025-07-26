use crate::{
    Axis::*,
    Irrep::{self, *},
    Molecule,
    plane::Plane,
    point_group::PointGroup,
};
use std::str::FromStr;

#[test]
fn c3h2() {
    let mol_orig = Molecule::from_str(
        "
    C        0.000000   -0.888844    0.000000
    C       -0.662697    0.368254    0.000000
    C        0.662697    0.368254    0.000000
    H       -1.595193    0.906925    0.000000
    H        1.595193    0.906925    0.000000
",
    )
    .unwrap();
    let pg = mol_orig.point_group();
    let tests = vec![
        (
            vec![
                0.00, 0.03, 0.00, 0.22, -0.11, 0.00, -0.22, -0.11, 0.00, -0.57,
                0.33, 0.00, 0.57, 0.33, 0.00,
            ],
            A1,
        ),
        (
            vec![
                -0.01, 0.00, 0.00, 0.17, -0.10, 0.00, 0.17, 0.10, 0.00, -0.59,
                0.34, 0.00, -0.59, -0.34, 0.00,
            ],
            B2,
        ),
        (
            vec![
                0.00, 0.00, 0.00, 0.00, 0.00, 0.40, 0.00, 0.00, -0.40, 0.00,
                0.00, -0.58, 0.00, 0.00, 0.58,
            ],
            A2,
        ),
        (
            vec![
                0.00, 0.00, 0.16, 0.00, 0.00, -0.27, 0.00, 0.00, -0.27, 0.00,
                0.00, 0.64, 0.00, 0.00, 0.64,
            ],
            B1,
        ),
    ];
    for test in tests {
        let mol = mol_orig.clone() + test.0;
        assert_eq!(mol.irrep(&pg), test.1);
    }
}

#[test]
fn c2h4() {
    let mol = Molecule::from_str(
        "
C        0.0000000000       -0.0023898386        1.2600838751
H        0.0000000000        1.7483464088        2.3303799608
H        0.0000000000       -1.7425505916        2.3220592143
C        0.0000000000       -0.0014113004       -1.2600853510
H        0.0000000000        1.7444525133       -2.3255411215
H        0.0000000000       -1.7464471915       -2.3268965777
",
    )
    .unwrap();
    assert_eq!(
        mol.irrep(&PointGroup::C2v {
            axis: Z,
            planes: [Plane(X, Z), Plane(Y, Z)]
        }),
        B2
    );
}

#[test]
fn c2h4_10() {
    let mol = Molecule::from_str(
        "
C      0.00139327 -1.25400055 -0.00000000
H     -0.00036672 -2.38212447  1.70464007
H     -0.00036672 -2.38212448 -1.70464006
C     -0.00139327  1.25400055  0.00000000
H      0.00036672  2.38212448  1.70464006
H      0.00036672  2.38212448 -1.70464007
",
    )
    .unwrap();
    assert_eq!(
        mol.irrep(&PointGroup::C2v {
            axis: Y,
            planes: [Plane(X, Y), Plane(Y, Z)]
        }),
        B1
    );
}

#[test]
fn c2h4_again() {
    let mol = Molecule::from_str(
        "C     -0.00139327 -0.00000000 -1.25400055
H      0.00036672  1.70464007 -2.38212447
H      0.00036672 -1.70464006 -2.38212448
C      0.00139327  0.00000000  1.25400055
H     -0.00036672  1.70464006  2.38212448
H     -0.00036672 -1.70464007  2.38212448
",
    )
    .unwrap();
    assert_eq!(
        mol.irrep(&PointGroup::C2v {
            axis: Z,
            planes: [Plane(X, Z), Plane(Y, Z)]
        }),
        B1
    );
}

#[test]
fn c2h4_again_again() {
    let mol = Molecule::from_str(
        "
C     -0.00146023  0.00000002  1.36220367
H      0.00039049 -1.76535618  2.54697919
H      0.00039049  1.76535607  2.54697922
C      0.00146023 -0.00000001 -1.36220368
H     -0.00039049 -1.76535607 -2.54697922
H     -0.00039049  1.76535618 -2.54697919",
    )
    .unwrap();
    assert_eq!(
        mol.irrep_approx(
            &PointGroup::C2v {
                axis: Z,
                planes: [Plane(X, Z), Plane(Y, Z)]
            },
            1e-6
        )
        .unwrap(),
        B1
    );
}

#[test]
fn c2h4_d2h() {
    struct Test {
        mol: &'static str,
        pg: PointGroup,
        eps: f64,
        want: Irrep,
    }
    let tests = [
        Test {
            mol: "
C      0.66679330  0.00000000 -0.42664810
H      1.23098540 -0.92361100  0.39872610
H      1.23098540  0.92361100  0.39872610
C     -0.66679330  0.00000000  0.42665160
H     -1.23098540 -0.92361100 -0.39873220
H     -1.23098540  0.92361100 -0.39873220
",
            pg: PointGroup::D2h {
                axes: [Z, X, Y],
                planes: [Plane(X, Y), Plane(Y, Z), Plane(X, Z)],
            },
            eps: 1e-5,
            want: B3g,
        },
        Test {
            mol: "
C     -0.70962050  0.00000000 -0.26815500
H     -1.36750740 -0.95203650  0.46264510
H     -1.36750740  0.95203650  0.46265790
C      0.70962050  0.00000000 -0.26815500
H      1.36750740 -0.95203650  0.46265790
H      1.36750740  0.95203650  0.46264510
",
            pg: PointGroup::D2h {
                axes: [X, Y, Z],
                planes: [Plane(Y, Z), Plane(X, Z), Plane(X, Y)],
            },
            eps: 1e-4,
            want: B3u,
        },
        Test {
            mol: "
C     -0.75776430  0.00000000 -0.26815500
H     -1.39084240  0.94668350  0.46270250
H     -1.39084240 -0.94668350  0.46260040
C      0.75776430  0.00000000 -0.26815510
H      1.39084240  0.94668350  0.46260050
H      1.39084240 -0.94668350  0.46270260
",
            pg: PointGroup::D2h {
                axes: [X, Y, Z],
                planes: [Plane(Y, Z), Plane(X, Z), Plane(X, Y)],
            },
            eps: 1e-3,
            want: B3u,
        },
    ];
    for test in tests {
        let mol = Molecule::from_str(test.mol).unwrap();
        assert_eq!(mol.irrep_approx(&test.pg, test.eps,).unwrap(), test.want);
    }
}
