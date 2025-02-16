use serde::{Deserialize, Serialize};

use crate::Axis;

mod ops;

// restrict these to combinations of cartesian axes for now. a more general
// plane is described by (a, b, c) in the equation ax + by + cz = 0
#[derive(Debug, PartialEq, Eq, Copy, Clone, Serialize, Deserialize)]
pub struct Plane(pub Axis, pub Axis);

impl Plane {
    /// return a normalized version of Plane
    pub(crate) fn new(ax: Axis, bx: Axis) -> Self {
        use Axis::*;
        match (ax, bx) {
            (X, Y) | (Y, X) => Plane(X, Y),
            (X, Z) | (Z, X) => Plane(X, Z),
            (Y, Z) | (Z, Y) => Plane(Y, Z),
            _ => panic!("impossible Axis combination for Plane"),
        }
    }

    /// return the axis perpendicular to `self`
    pub fn perp(&self) -> Axis {
        let Plane(ax, bx) = self;
        use Axis::*;
        match (ax, bx) {
            (X, Y) | (Y, X) => Z,
            (X, Z) | (Z, X) => Y,
            (Y, Z) | (Z, Y) => X,
            _ => panic!("impossible Axis combination for Plane"),
        }
    }
}
