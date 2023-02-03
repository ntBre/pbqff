use serde::{Deserialize, Serialize};

mod coord_type;
pub use coord_type::*;

#[derive(Deserialize, Debug, PartialEq)]
#[serde(deny_unknown_fields)]
struct RawConfig {
    /// the geometry to start with
    geometry: String,

    /// whether or not to optimize the structure first
    optimize: bool,

    /// charge on the molecule
    charge: isize,

    /// distance in Å to displace the atoms
    step_size: f64,

    /// how long to sleep between intervals polling running jobs
    sleep_int: usize,

    /// limit for the number of jobs to run at once
    job_limit: usize,

    /// the number of jobs to include a single submit script
    chunk_size: usize,

    /// whether to use SICs or Cartesian coordinates
    coord_type: CoordType,

    /// the template to use for the quantum chemistry program
    template: String,

    /// the quantum chemistry program to use. options supported currently are
    /// mopac and molpro
    program: Program,

    /// the type of queuing system to use. options supported currently are pbs
    /// and slurm
    queue: Queue,

    /// whether to use finite differences or a least-squares fitting to obtain
    /// the force constants. currently this only affects normal coordinates.
    /// defaults to false
    findiff: Option<bool>,

    /// interval for dumping checkpoint files. 0 means no checkpoints
    check_int: usize,
}

#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Clone, Copy)]
pub enum Program {
    #[serde(alias = "mopac")]
    Mopac,
    #[serde(alias = "molpro")]
    Molpro,
}

#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Clone, Copy)]
pub enum Queue {
    #[serde(alias = "pbs")]
    Pbs,
    #[serde(alias = "slurm")]
    Slurm,
    #[serde(alias = "local")]
    Local,
}

#[derive(Clone, Serialize, Deserialize, PartialEq, Debug)]
#[serde(from = "RawConfig")]
pub struct Config {
    pub geometry: psqs::geom::Geom,
    pub optimize: bool,
    pub charge: isize,
    pub step_size: f64,
    pub coord_type: CoordType,
    pub template: String,
    pub program: Program,
    pub queue: Queue,
    pub sleep_int: usize,
    pub job_limit: usize,
    pub chunk_size: usize,
    pub findiff: bool,
    pub check_int: usize,
}

impl From<RawConfig> for Config {
    fn from(rc: RawConfig) -> Self {
        Self {
            geometry: rc.geometry.parse().unwrap(),
            optimize: rc.optimize,
            charge: rc.charge,
            step_size: rc.step_size,
            coord_type: rc.coord_type,
            template: rc.template,
            program: rc.program,
            sleep_int: rc.sleep_int,
            job_limit: rc.job_limit,
            chunk_size: rc.chunk_size,
            queue: rc.queue,
            findiff: rc.findiff.unwrap_or(false),
            check_int: rc.check_int,
        }
    }
}

impl Config {
    pub fn load(filename: &str) -> Self {
        let contents = std::fs::read_to_string(filename)
            .expect("failed to load config file");
        toml::from_str(&contents).unwrap_or_else(|e| {
            panic!("failed to deserialize config file '{filename}' with {e}")
        })
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
        let got = Config::load("testfiles/test.toml");
        let want = Config {
            geometry: psqs::geom::Geom::Zmat(
                "C
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
            ),
            optimize: true,
            charge: 0,
            step_size: 0.005,
            coord_type: CoordType::Sic,
            template: String::from(
                "scfcrt=1.D-21 aux(precision=14 comp xp xs xw) PM6 THREADS=1 \
		 external=testfiles/params.dat",
            ),
            program: Program::Mopac,
            sleep_int: 2,
            job_limit: 2048,
            chunk_size: 1,
            queue: Queue::Slurm,
            findiff: false,
            check_int: 100,
        };
        assert_eq!(got, want);
    }
}
