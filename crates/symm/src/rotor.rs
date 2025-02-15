use std::fmt::Display;

use serde::{Deserialize, Serialize};

#[macro_export]
macro_rules! unset_rotor {
    () => {
        panic!("rotor type not set")
    };
}

// TODO probably need to get rid of default here. it will complicate the spectro
// usage a little bit but then I don't have to have the None variant and keep
// checking for it
#[derive(Clone, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum Rotor {
    Diatomic,
    Linear,
    SphericalTop,
    OblateSymmTop,
    ProlateSymmTop,
    AsymmTop,
    #[default]
    None,
}

impl Rotor {
    /// panics if `self` is not set
    pub fn is_prolate(&self) -> bool {
        assert!(*self != Rotor::None);
        *self == Rotor::ProlateSymmTop
    }

    /// Report whether or not `self` is either an `OblateSymmTop` or a
    /// `ProlateSymmTop`. panics if `self` is not set
    pub fn is_sym_top(&self) -> bool {
        use Rotor::*;
        match &self {
            Diatomic => false,
            OblateSymmTop | ProlateSymmTop | Linear => true,
            SphericalTop | AsymmTop => false,
            None => unset_rotor!(),
        }
    }

    /// Returns `true` if the rotor is [`Linear`].
    ///
    /// [`Linear`]: Rotor::Linear
    #[must_use]
    pub fn is_linear(&self) -> bool {
        matches!(self, Self::Linear | Self::Diatomic)
    }

    /// Returns `true` if the rotor is [`OblateSymmTop`].
    ///
    /// [`OblateSymmTop`]: Rotor::OblateSymmTop
    #[must_use]
    pub fn is_oblate(&self) -> bool {
        matches!(self, Self::OblateSymmTop)
    }

    /// Returns `true` if the rotor is [`ProlateSymmTop`].
    ///
    /// [`ProlateSymmTop`]: Rotor::ProlateSymmTop
    #[must_use]
    pub fn is_prolate_symm_top(&self) -> bool {
        matches!(self, Self::ProlateSymmTop)
    }

    /// Returns `true` if the rotor is [`Diatomic`].
    ///
    /// [`Diatomic`]: Rotor::Diatomic
    #[must_use]
    pub fn is_diatomic(&self) -> bool {
        matches!(self, Self::Diatomic)
    }

    /// Returns `true` if the rotor is [`SphericalTop`].
    ///
    /// [`SphericalTop`]: Rotor::SphericalTop
    #[must_use]
    pub fn is_spherical_top(&self) -> bool {
        matches!(self, Self::SphericalTop)
    }

    /// Returns `true` if the rotor is [`AsymmTop`].
    ///
    /// [`AsymmTop`]: Rotor::AsymmTop
    #[must_use]
    pub fn is_asymm_top(&self) -> bool {
        matches!(self, Self::AsymmTop)
    }
}

impl Display for Rotor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Rotor::Diatomic => "diatomic",
                Rotor::Linear => "linear",
                Rotor::SphericalTop => "a spherical top",
                Rotor::OblateSymmTop => "an oblate symmetric top",
                Rotor::ProlateSymmTop => "a prolate symmetric top",
                Rotor::AsymmTop => "an asymmetric top",
                Rotor::None => unset_rotor!(),
            }
        )
    }
}
