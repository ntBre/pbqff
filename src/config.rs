use serde::Deserialize;

mod coord_type;
pub use coord_type::*;

#[derive(Deserialize, Debug, PartialEq)]
struct RawConfig {
    /// the geometry to start with
    geometry: String,

    /// whether or not to optimize the structure first
    optimize: bool,

    /// charge on the molecule
    charge: isize,

    /// distance in Å to displace the atoms
    step_size: f64,

    /// whether to use SICs or Cartesian coordinates
    coord_type: CoordType,

    /// the template to use for the quantum chemistry program
    template: String,
}

impl RawConfig {
    fn load(filename: &str) -> Self {
        let contents = std::fs::read_to_string(filename)
            .expect("failed to load config file");
        toml::from_str(&contents).expect("failed to deserialize config file")
    }
}

pub struct Config {
    pub geometry: psqs::geom::Geom,
    pub optimize: bool,
    pub charge: isize,
    pub step_size: f64,
    pub coord_type: CoordType,
    pub template: String,
}

impl Config {
    pub fn load(filename: &str) -> Self {
        let rc = RawConfig::load(filename);
        Self {
            geometry: rc.geometry.parse().unwrap(),
            optimize: rc.optimize,
            charge: rc.charge,
            step_size: rc.step_size,
            coord_type: rc.coord_type,
            template: rc.template,
        }
    }
}

impl std::fmt::Display for Config {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "
Configuration Options:
geometry = {{
{}
}}
optimize = {}
charge = {}
step_size = {}
coord_type = {}
template = {}
",
            self.geometry.to_string().trim(),
            self.optimize,
            self.charge,
            self.step_size,
            self.coord_type,
            self.template,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn config() {
        let got = RawConfig::load("testfiles/test.toml");
        let want = RawConfig {
            geometry: "C
C 1 CC
C 1 CC 2 CCC
H 2 CH 1 HCC 3 180.0
H 3 CH 1 HCC 2 180.0

CC =                  1.42101898
CCC =                55.60133141
CH =                  1.07692776
HCC =               147.81488230
"
            .to_string(),
            optimize: true,
            charge: 0,
            step_size: 0.005,
            coord_type: CoordType::sic,
	    template: String::from(
		"A0 scfcrt=1.D-21 aux(precision=14) PM6 external=testfiles/params.dat",
		),
        };
        assert_eq!(got, want);
    }
}
