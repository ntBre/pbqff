//! Configuration settings for running a pbqff

use std::fmt::Display;

use serde::{Deserialize, Serialize};

mod coord_type;
pub use coord_type::*;

#[cfg(test)]
mod tests;

#[derive(Deserialize, Debug, PartialEq)]
#[serde(deny_unknown_fields)]
struct RawConfig {
    geometry: String,
    optimize: bool,
    charge: isize,
    step_size: f64,
    sleep_int: usize,
    job_limit: usize,
    chunk_size: usize,
    coord_type: CoordType,
    template: String,
    hybrid_template: Option<String>,
    queue_template: Option<String>,
    program: Program,
    queue: Queue,
    findiff: Option<bool>,
    check_int: usize,
}

#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Clone, Copy)]
pub enum Program {
    #[serde(alias = "mopac")]
    Mopac,
    #[serde(alias = "molpro")]
    Molpro,
}

impl Display for Program {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Program::Mopac => write!(f, "mopac"),
            Program::Molpro => write!(f, "molpro"),
        }
    }
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

impl Display for Queue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Queue::Pbs => "pbs",
                Queue::Slurm => "slurm",
                Queue::Local => "local",
            }
        )
    }
}

#[derive(Clone, Serialize, Deserialize, PartialEq, Debug)]
#[serde(from = "RawConfig")]
pub struct Config {
    /// the geometry to start with. Parsed from a string using
    /// [psqs::geom::Geom]'s implementation of [std::str::FromStr]
    pub geometry: psqs::geom::Geom,

    /// whether or not to optimize the structure before running the QFF
    pub optimize: bool,

    /// charge on the molecule
    pub charge: isize,

    /// distance in Å to displace the atoms
    pub step_size: f64,

    /// whether to use SICs or Cartesian coordinates
    pub coord_type: CoordType,

    /// the template to use for the quantum chemistry program
    pub template: String,

    /// optional quantum chemistry program template to use for the cubic and
    /// quartic force constants. only used for finite-difference normal
    /// coordinate QFFs. If not supplied, this defaults to the same value as
    /// `template`
    pub hybrid_template: String,

    /// the optional template to use for the queuing system. If this is not
    /// provided, the queue's implementation of
    /// [psqs::queue::Queue::default_submit_script] will be used
    pub queue_template: Option<String>,

    /// the quantum chemistry program to use. options supported currently are
    /// mopac and molpro, as deserialized from [Program]
    pub program: Program,

    /// the type of queuing system to use. options supported currently are pbs
    /// and slurm, as deserialized from [Queue]
    pub queue: Queue,

    /// how long to sleep between intervals polling running jobs
    pub sleep_int: usize,

    /// limit for the number of jobs to run at once
    pub job_limit: usize,

    /// the number of jobs to include a single submit script
    pub chunk_size: usize,

    /// whether to use finite differences or a least-squares fitting to obtain
    /// the force constants. currently this only affects normal coordinates.
    /// defaults to false
    pub findiff: bool,

    /// interval for dumping checkpoint files. 0 means no checkpoints
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
            hybrid_template: rc
                .hybrid_template
                .unwrap_or_else(|| rc.template.clone()),
            template: rc.template,
            program: rc.program,
            sleep_int: rc.sleep_int,
            job_limit: rc.job_limit,
            chunk_size: rc.chunk_size,
            queue: rc.queue,
            findiff: rc.findiff.unwrap_or(false),
            check_int: rc.check_int,
            queue_template: rc.queue_template,
        }
    }
}

impl Config {
    /// load a [Config] from the TOML file specified by `filename`. panics on
    /// failure to read the file and on failure to deserialize it
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
        let Config {
            geometry,
            optimize,
            charge,
            step_size,
            coord_type,
            template,
            hybrid_template,
            queue_template,
            program,
            queue,
            sleep_int,
            job_limit,
            chunk_size,
            findiff,
            check_int,
        } = self;
        write!(
            f,
            "
Configuration Options:
geometry = {{
{geometry}
}}
optimize = {optimize}
charge = {charge}
step_size = {step_size}
coord_type = {coord_type}
template = {template}
hybrid_template = {hybrid_template}
queue_template = {}
program = {program}
queue = {queue}
sleep_int = {sleep_int}
job_limit = {job_limit}
chunk_size = {chunk_size}
findiff = {findiff}
check_int = {check_int}
",
            queue_template.as_ref().unwrap_or(&String::new()),
        )
    }
}
