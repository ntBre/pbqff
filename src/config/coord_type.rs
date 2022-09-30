use serde::Deserialize;
use std::fmt::Display;

#[derive(Deserialize, Debug, PartialEq, Eq, Clone, Copy)]
pub enum CoordType {
    #[serde(alias = "cart")]
    Cart,
    #[serde(alias = "sic")]
    Sic,
}

impl CoordType {
    /// Returns `true` if the coord type is [`cart`].
    ///
    /// [`cart`]: CoordType::cart
    #[must_use]
    pub fn is_cart(&self) -> bool {
        matches!(self, Self::Cart)
    }

    /// Returns `true` if the coord type is [`sic`].
    ///
    /// [`sic`]: CoordType::sic
    #[must_use]
    pub fn is_sic(&self) -> bool {
        matches!(self, Self::Sic)
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
            }
        )
    }
}
