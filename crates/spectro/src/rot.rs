use std::fmt::Display;

use crate::state::State;

/// a rotational constant. eventually this will probably be an enum depending on
/// the Rotor type
#[derive(Clone, Debug, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Rot {
    pub state: State,
    pub a: f64,
    pub b: f64,
    pub c: f64,
}

impl Rot {
    pub fn new(state: State, a: f64, b: f64, c: f64) -> Self {
        Self { state, a, b, c }
    }
}

#[cfg(test)]
impl approx::AbsDiffEq for Rot {
    type Epsilon = f64;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        self.state == other.state
            && self.a.abs_diff_eq(&other.a, epsilon)
            && self.b.abs_diff_eq(&other.b, epsilon)
            && self.c.abs_diff_eq(&other.c, epsilon)
    }
}

impl Display for Rot {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}{:12.7}{:12.7}{:12.7}",
            self.state, self.a, self.b, self.c
        )
    }
}
