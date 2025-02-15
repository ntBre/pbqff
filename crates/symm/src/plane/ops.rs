use super::Plane;
use crate::Axis;
use std::fmt::Display;
use std::ops::BitXor;

impl BitXor<Axis> for Plane {
    type Output = Axis;

    fn bitxor(self, rhs: Axis) -> Self::Output {
        let Plane(ax, bx) = self;
        use Axis::*;
        match (ax, bx) {
            (X, Y) | (Y, X) => match rhs {
                X => Y,
                Y => X,
                Z => panic!("Z not in XY"),
            },
            (X, Z) | (Z, X) => match rhs {
                X => Z,
                Y => panic!("Y not in XZ"),
                Z => X,
            },
            (Y, Z) | (Z, Y) => match rhs {
                X => panic!("X not in YZ"),
                Y => Z,
                Z => Y,
            },
            _ => panic!("impossible Axis combination for Plane"),
        }
    }
}

impl Display for Plane {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Plane({}, {})", self.0, self.1)
    }
}
