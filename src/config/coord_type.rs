use serde::{Deserialize, Serialize};
use std::fmt::Display;

#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Clone, Copy)]
pub enum CoordType {
    #[serde(alias = "cart")]
    Cart,
    #[serde(alias = "sic")]
    Sic,
    #[serde(alias = "normal")]
    Normal,
}

impl CoordType {
    /// Returns `true` if the coord type is [`Cart`].
    ///
    /// [`Cart`]: CoordType::Cart
    #[must_use]
    pub fn is_cart(&self) -> bool {
        matches!(self, Self::Cart)
    }

    /// Returns `true` if the coord type is [`Sic`].
    ///
    /// [`Sic`]: CoordType::Sic
    #[must_use]
    pub fn is_sic(&self) -> bool {
        matches!(self, Self::Sic)
    }

    /// Returns `true` if the coord type is [`Normal`].
    ///
    /// [`Normal`]: CoordType::Normal
    #[must_use]
    pub fn is_normal(&self) -> bool {
        matches!(self, Self::Normal)
    }
}

impl Display for CoordType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                CoordType::Cart => "cart",
                CoordType::Sic => "sic",
                CoordType::Normal => "normal",
            }
        )
    }
}
