//! Configuration settings for running a pbqff

use std::{
    fmt::{Debug, Display},
    fs::read_to_string,
    path::Path,
};

use serde::{Deserialize, Serialize};

mod coord_type;
pub use coord_type::*;

#[cfg(test)]
mod tests;

/// Templates can either be literal strings in the config file, or the name of a
/// file to be loaded
#[derive(Clone, Serialize, Deserialize, PartialEq, Debug)]
#[serde(untagged)]
enum TemplateSrc {
    Literal(String),
    File { file: String },
}

impl From<TemplateSrc> for String {
    fn from(value: TemplateSrc) -> Self {
        match value {
            TemplateSrc::Literal(s) => s,
            TemplateSrc::File { file } => {
                read_to_string(&file).unwrap_or_else(|_| {
                    panic!("failed to locate template file {file}")
                })
            }
        }
    }
}

#[derive(Deserialize, Debug, PartialEq)]
#[serde(deny_unknown_fields)]
struct RawConfig {
    /// The initial geometry to use for the computation. Both XYZ and Z-matrix
    /// geometries are accepted.
    geometry: String,

    /// Whether or not to perform a geometry optimization on the input
    optimize: bool,

    /// The molecular charge. This value can be spliced into the template using
    /// the {{.charge}} directive.
    charge: isize,

    /// The size of the displacement to take in the QFF coordinates.
    step_size: f64,

    /// The interval in seconds to wait between loops checking if any jobs have
    /// finished.
    sleep_int: usize,

    /// The maximum number of jobs to submit at once, as determined by the
    /// number of individual input files. This distinction is important when
    /// chunk_size is greater than 1 because the maximum number of jobs
    /// submitted to the queue will be job_limit / chunk_size .
    job_limit: usize,

    /// The number of individual calculations to bundle into a single queue
    /// submission.
    chunk_size: usize,

    /// The type of coordinate to use in the QFF. Currently-supported values are
    /// "sic", "cart", and "normal". Note that SIC QFFs require an additional
    /// input file called intder.in to define the internal coordinates.
    coord_type: CoordType,

    /// The template input file for the quantum chemistry program. Supported
    /// formatting directives depend on the program in question. Molpro supports
    /// {{.geom}} for the geometry and {{.charge}} for the molecular charge,
    /// while Mopac expects a static template.
    template: TemplateSrc,

    /// An optional template file for the cubic and quartic portion of the QFF.
    /// If this is provided, the regular template is used only for the harmonic
    /// portion of the QFF, and this is used for the rest of the points.
    hybrid_template: Option<TemplateSrc>,

    /// The template input file for the queuing system. Supported formatting
    /// directives include {{.basename}} for the base name of the submit script,
    /// which is useful for naming the job in the queue, and {{.filename}} for
    /// the name of the quantum chemistry program input file. Not all program
    /// and queue combinations expand these, however.
    queue_template: Option<TemplateSrc>,

    /// The quantum chemistry program to use in running the QFF.
    /// Currently-supported values are "cfour", "dftb+", "molpro", "mopac".
    program: Program,

    /// The queuing system to use in running the QFF. Currently-supported values
    /// are "local", which uses bash to run computations directly, "pbs", and
    /// "slurm".
    queue: Queue,

    /// Whether to use finite differences or least-squares fitting for the
    /// potential energy surface. Currently normal coordinates are the only
    /// coord_type to use this option, so it has a default value of false,
    /// meaning use the fitted version of normal coordinates. Setting this
    /// option to true forces the use of finite differences for the normal
    /// coordinate QFF.
    findiff: Option<bool>,

    /// The interval at which to write checkpoint files. Every coordinate type
    /// will write an initial checkpoint (res.chk), but this interval determines
    /// whether or not checkpoints are written while the single-point energies
    /// are being run. A value of 0 will disable checkpoints entirely. This
    /// interval refers to the number of polling iterations that have occurred,
    /// not the number of jobs that have completed. The iteration count is shown
    /// in the log file when the number of jobs remaining is printed.
    check_int: usize,

    /// An optional vector of atomic masses to use for normal coordinate
    /// generation and Spectro.
    weights: Option<Vec<f64>>,

    /// An optional number of atoms to hold constant in the QFF displacements.
    /// These must come at the end of the geometry. Experimental
    dummy_atoms: Option<usize>,
}

#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Clone, Copy)]
pub enum Program {
    #[serde(alias = "dftb+")]
    DFTBPlus,
    #[serde(alias = "mopac")]
    Mopac,
    #[serde(alias = "molpro")]
    Molpro,
    #[serde(alias = "cfour", alias = "CFOUR")]
    Cfour,
}

impl Display for Program {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Program::Mopac => write!(f, "mopac"),
            Program::Molpro => write!(f, "molpro"),
            Program::DFTBPlus => write!(f, "dftb+"),
            Program::Cfour => write!(f, "cfour"),
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

/// Construct a full `Config` using [Config::load] on a TOML file or use
/// [Config::new] and the Builder pattern
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

    /// atomic masses to use for the QFF. really only matters for normal
    /// coordinates. must have the same length as the number of atoms, but this
    /// is not currently validated
    pub weights: Option<Vec<f64>>,

    pub dummy_atoms: Option<usize>,
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
                .unwrap_or_else(|| rc.template.clone())
                .into(),
            template: rc.template.into(),
            program: rc.program,
            sleep_int: rc.sleep_int,
            job_limit: rc.job_limit,
            chunk_size: rc.chunk_size,
            queue: rc.queue,
            findiff: rc.findiff.unwrap_or(false),
            check_int: rc.check_int,
            queue_template: rc.queue_template.map(TemplateSrc::into),
            weights: rc.weights,
            dummy_atoms: rc.dummy_atoms,
        }
    }
}

macro_rules! int_builders {
    ($($name: ident$(,)*)*) => {
        $(pub fn $name(mut self, i: usize) -> Self {
            self.$name = i;
            self
        })*
    }
}

impl Config {
    /// Construct a [Config] with default values for `step_size` (0.005 Å),
    /// `hybrid_template` (`template`), `queue_template` (`None`), `sleep_int`
    /// (1 second), `job_limit` (128), `chunk_size` (1), `findiff` (`true`), and
    /// `check_int` (0; disabled)
    pub fn new(
        geometry: psqs::geom::Geom,
        optimize: bool,
        charge: isize,
        coord_type: CoordType,
        template: String,
        program: Program,
        queue: Queue,
    ) -> Self {
        Self {
            geometry,
            optimize,
            charge,
            step_size: 0.005,
            coord_type,
            template: template.clone(),
            hybrid_template: template,
            queue_template: None,
            program,
            queue,
            sleep_int: 1,
            job_limit: 128,
            chunk_size: 1,
            findiff: true,
            check_int: 0,
            weights: None,
            dummy_atoms: None,
        }
    }

    pub fn step_size(mut self, step_size: f64) -> Self {
        self.step_size = step_size;
        self
    }

    pub fn hybrid_template(mut self, t: String) -> Self {
        self.hybrid_template = t;
        self
    }

    pub fn queue_template(mut self, t: Option<String>) -> Self {
        self.queue_template = t;
        self
    }

    int_builders!(sleep_int, job_limit, chunk_size, check_int);

    pub fn findiff(mut self, b: bool) -> Self {
        self.findiff = b;
        self
    }

    /// load a [Config] from the TOML file specified by `filename`. panics on
    /// failure to read the file and on failure to deserialize it.
    ///
    /// TODO return a result with real error handling
    pub fn load<P>(filename: P) -> Self
    where
        P: AsRef<Path> + Debug,
    {
        let contents = std::fs::read_to_string(&filename)
            .expect("failed to load config file");
        let ret: Self = toml::from_str(&contents).unwrap_or_else(|e| {
            panic!("failed to deserialize config file '{filename:?}' with {e}")
        });

        ret.validate();

        ret
    }

    /// check that the settings in `self` make any sense. currently only check
    /// that job_limit is greater than chunk size
    fn validate(&self) {
        if self.job_limit < self.chunk_size {
            eprintln!(
                "In pbqff.toml: Your job_limit ({}) is TOO LOW. \
                 Must be greater than chunk_size ({}), exiting",
                self.job_limit, self.chunk_size
            );
            std::process::exit(1);
        }
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
            weights,
            dummy_atoms,
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
dummy_atoms = {dummy_atoms:?}
",
            queue_template.as_ref().unwrap_or(&String::new()),
        )?;
        if let Some(ws) = weights {
            write!(f, "weights = [ ")?;
            for w in ws {
                write!(f, "{w}, ")?;
            }
            writeln!(f, "]")?;
        }
        Ok(())
    }
}
