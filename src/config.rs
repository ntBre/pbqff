use serde::Deserialize;

#[derive(Deserialize, Debug, PartialEq)]
struct RawConfig {
    /// the geometry to start with
    geometry: String,

    /// whether or not to optimize the structure first
    optimize: bool,

    /// the path to the wrapper SPECTRO program to be called, eventually this
    /// will be deprecated when I have a Rust version of gspectro
    gspectro_cmd: String,

    /// the path to the real SPECTRO executable to be passed to gspectro,
    /// deprecated when I rewrite SPECTRO
    spectro_cmd: String,

    /// charge on the molecule
    charge: isize,

    /// distance in Å to displace the atoms
    step_size: f64,
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
    pub gspectro_cmd: String,
    pub spectro_cmd: String,
    pub charge: isize,
    pub step_size: f64,
}

impl Config {
    pub fn load(filename: &str) -> Self {
        let rc = RawConfig::load(filename);
        Self {
            geometry: rc.geometry.parse().unwrap(),
            optimize: rc.optimize,
            gspectro_cmd: rc.gspectro_cmd,
            spectro_cmd: rc.spectro_cmd,
            charge: rc.charge,
            step_size: rc.step_size,
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
gspectro_cmd = {}
spectro_cmd = {}
charge = {}
step_size = {}
",
            self.geometry.to_string().trim(),
            self.optimize,
            self.gspectro_cmd,
            self.spectro_cmd,
            self.charge,
            self.step_size
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
            gspectro_cmd:
                "/home/brent/Projects/chemutils/spectro/spectro/spectro"
                    .to_string(),
            spectro_cmd: "/home/brent/Projects/pbqff/bin/spectro".to_string(),
            charge: 0,
            step_size: 0.005,
        };
        assert_eq!(got, want);
    }
}
