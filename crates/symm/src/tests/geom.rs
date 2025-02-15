//! tests for geometrical operations like the center of mass and the moment of
//! inertia

use std::str::FromStr;

use crate::*;
use approx::assert_abs_diff_eq;

#[test]
fn com() {
    let mol = Molecule::from_str(
        "
    			H 0.0000000000 1.4313901416 0.9860410955
			O 0.0000000000 0.0000000000 -0.1242384417
			H 0.0000000000 -1.4313901416 0.9860410955
",
    )
    .unwrap();
    let got = mol.com();
    let want = Vec3::from_row_slice(&[
        0.0000000,
        0.0000000,
        9.711_590_454_604_224e-6 / 0.52917706,
    ]);
    assert_abs_diff_eq!(got, want, epsilon = 1e-8);
}

#[test]
fn inertia_tensor() {
    let mut mol = Molecule::from_str(
        "
    			H 0.0000000000 1.4313901416 0.9860410955
			O 0.0000000000 0.0000000000 -0.1242384417
			H 0.0000000000 -1.4313901416 0.9860410955
",
    )
    .unwrap();
    let com = mol.com();
    mol.translate(-com);
    let got = mol.moi();
    let want = na::matrix![
        1.7743928167251328, 0.0, 0.0;
        0.0, 0.617_925_936_198_827_2, 0.0;
        0.0, 0.0, 1.1564668805263056;
    ] / 0.52917706
        / 0.52917706;
    assert_abs_diff_eq!(got, want, epsilon = 1e-7);
}

#[test]
fn normalize() {
    let mut got = molecule![
        C 0.0000000000 0.0000000000 -0.5592657284
        N 0.0000000000 0.0000000000 0.5966002840
        H 0.0000000000 0.0000000000 -1.6261489121
    ];
    let want = molecule![
    C      0.00000000  0.00000000  0.55942032
    N      0.00000000  0.00000000 -0.59644569
    H      0.00000000  0.00000000  1.62630350
        ];
    got.normalize();
    assert_eq!(got, want);
}

#[test]
fn normalize_diatomic() {
    let mut got = molecule![
        H      0.00000000  0.00000000  0.00000000
        O      0.00000000  0.00000000  1.30998562
    ];
    let want = molecule![
        H     -1.23233718  0.00000000  0.00000000
        O      0.07764845  0.00000000  0.00000000
    ];
    got.normalize();
    assert_eq!(got, want);

    got.normalize();
    assert_eq!(got, want);
}
