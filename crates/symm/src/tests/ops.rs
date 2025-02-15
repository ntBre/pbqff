//! tests for symmetry operations like rotations and reflections

use crate::*;

#[test]
fn rotate() {
    use Axis::*;
    let tests = vec![
        // X around all axes
        (
            vec![Atom::new(1, 1.0, 0.0, 0.0)],
            vec![Atom::new(1, 1.0, 0.0, 0.0)],
            180.0,
            X,
        ),
        (
            vec![Atom::new(1, 1.0, 0.0, 0.0)],
            vec![Atom::new(1, -1.0, 0.0, 0.0)],
            180.0,
            Y,
        ),
        (
            vec![Atom::new(1, 1.0, 0.0, 0.0)],
            vec![Atom::new(1, -1.0, 0.0, 0.0)],
            180.0,
            Z,
        ),
        // Y around all axes
        (
            vec![Atom::new(1, 0.0, 1.0, 0.0)],
            vec![Atom::new(1, 0.0, -1.0, 0.0)],
            180.0,
            X,
        ),
        (
            vec![Atom::new(1, 0.0, 1.0, 0.0)],
            vec![Atom::new(1, 0.0, 1.0, 0.0)],
            180.0,
            Y,
        ),
        (
            vec![Atom::new(1, 0.0, 1.0, 0.0)],
            vec![Atom::new(1, 0.0, -1.0, 0.0)],
            180.0,
            Z,
        ),
        // Z around all axes
        (
            vec![Atom::new(1, 0.0, 0.0, 1.0)],
            vec![Atom::new(1, 0.0, 0.0, -1.0)],
            180.0,
            X,
        ),
        (
            vec![Atom::new(1, 0.0, 0.0, 1.0)],
            vec![Atom::new(1, 0.0, 0.0, -1.0)],
            180.0,
            Y,
        ),
        (
            vec![Atom::new(1, 0.0, 0.0, 1.0)],
            vec![Atom::new(1, 0.0, 0.0, 1.0)],
            180.0,
            Z,
        ),
    ];
    for test in tests {
        let h = Molecule { atoms: test.0 };
        let want = Molecule { atoms: test.1 };
        let got = h.rotate(test.2, &test.3);
        assert_eq!(got, want);
    }
}

#[test]
fn reflect() {
    use Axis::*;
    let tests = vec![
        // X through all the planes
        (
            vec![Atom::new(1, 1.0, 0.0, 0.0)],
            vec![Atom::new(1, -1.0, 0.0, 0.0)],
            Plane(Y, Z),
        ),
        (
            vec![Atom::new(1, 1.0, 0.0, 0.0)],
            vec![Atom::new(1, 1.0, 0.0, 0.0)],
            Plane(X, Z),
        ),
        (
            vec![Atom::new(1, 1.0, 0.0, 0.0)],
            vec![Atom::new(1, 1.0, 0.0, 0.0)],
            Plane(X, Y),
        ),
        // Y through all the planes
        (
            vec![Atom::new(1, 0.0, 1.0, 0.0)],
            vec![Atom::new(1, 0.0, 1.0, 0.0)],
            Plane(Y, Z),
        ),
        (
            vec![Atom::new(1, 0.0, 1.0, 0.0)],
            vec![Atom::new(1, 0.0, -1.0, 0.0)],
            Plane(X, Z),
        ),
        (
            vec![Atom::new(1, 0.0, 1.0, 0.0)],
            vec![Atom::new(1, 0.0, 1.0, 0.0)],
            Plane(X, Y),
        ),
        // Z through all the planes
        (
            vec![Atom::new(1, 0.0, 0.0, 1.0)],
            vec![Atom::new(1, 0.0, 0.0, 1.0)],
            Plane(Y, Z),
        ),
        (
            vec![Atom::new(1, 0.0, 0.0, 1.0)],
            vec![Atom::new(1, 0.0, 0.0, 1.0)],
            Plane(X, Z),
        ),
        (
            vec![Atom::new(1, 0.0, 0.0, 1.0)],
            vec![Atom::new(1, 0.0, 0.0, -1.0)],
            Plane(X, Y),
        ),
    ];
    for test in tests {
        let h = Molecule { atoms: test.0 };
        let want = Molecule { atoms: test.1 };
        let got = h.reflect(&test.2);
        assert_eq!(got, want);
    }
}
