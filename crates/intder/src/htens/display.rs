use std::fmt::Display;

use super::Htens;

impl Display for Htens {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "h111={:8}", self.h111)?;
        write!(f, "h112={:8}", self.h112)?;
        write!(f, "h113={:8}", self.h113)?;
        write!(f, "h123={:8}", self.h123)?;
        write!(f, "h221={:8}", self.h221)?;
        write!(f, "h222={:8}", self.h222)?;
        write!(f, "h223={:8}", self.h223)?;
        write!(f, "h331={:8}", self.h331)?;
        write!(f, "h332={:8}", self.h332)?;
        write!(f, "h333={:8}", self.h333)?;
        write!(f, "h411={:8}", self.h411)?;
        write!(f, "h421={:8}", self.h421)?;
        write!(f, "h422={:8}", self.h422)?;
        write!(f, "h431={:8}", self.h431)?;
        write!(f, "h432={:8}", self.h432)?;
        write!(f, "h433={:8}", self.h433)?;
        write!(f, "h441={:8}", self.h441)?;
        write!(f, "h442={:8}", self.h442)?;
        write!(f, "h443={:8}", self.h443)?;
        write!(f, "h444={:8}", self.h444)?;

        Ok(())
    }
}
