use serde::Deserialize;

// TODO geometry should probably be processed more in the future. call this
// RawConfig for deserializing and then do some conversion after to a real
// Config
#[derive(Deserialize, Debug, PartialEq)]
pub struct Config {
    pub geometry: String,
    pub optimize: bool,
    pub spectro: String,
}

impl Config {
    pub fn load(filename: &str) -> Self {
        let contents = std::fs::read_to_string(filename)
            .expect("failed to load config file");
        toml::from_str(&contents).expect("failed to deserialize config file")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn config() {
        let got = Config::load("testfiles/test.toml");
        let want = Config {
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
            spectro: "/home/brent/Downloads/spec3jm/backup/spectro.x"
                .to_string(),
        };
        assert_eq!(got, want);
    }
}
