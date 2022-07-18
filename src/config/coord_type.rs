use serde::Deserialize;
use std::fmt::Display;

#[allow(non_camel_case_types)]
#[derive(Deserialize, Debug, PartialEq)]
pub enum CoordType {
    cart,
    sic,
}

impl CoordType {
    /// Returns `true` if the coord type is [`cart`].
    ///
    /// [`cart`]: CoordType::cart
    #[must_use]
    pub fn is_cart(&self) -> bool {
        matches!(self, Self::cart)
    }

    /// Returns `true` if the coord type is [`sic`].
    ///
    /// [`sic`]: CoordType::sic
    #[must_use]
    pub fn is_sic(&self) -> bool {
        matches!(self, Self::sic)
    }
}

impl Display for CoordType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                CoordType::cart => "cart",
                CoordType::sic => "sic",
            }
        )
    }
}
