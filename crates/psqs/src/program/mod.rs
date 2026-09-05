use std::{
    error::Error,
    fmt::{Debug, Display},
    path::Path,
    str::FromStr,
    time::SystemTime,
};

use serde::{Deserialize, Serialize};
use symm::Atom;

use crate::geom::Geom;

pub mod cfour;
pub mod dftbplus;
pub mod molpro;
pub mod mopac;
pub mod orca;

#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct ProgramResult {
    pub energy: f64,
    pub cart_geom: Option<Vec<Atom>>,
    pub time: f64,
}

#[derive(Debug, PartialEq, Eq)]
pub enum ProgramError {
    FileNotFound(String),
    ErrorInOutput(String),
    EnergyNotFound(String),
    EnergyParseError(String),
    GeomNotFound(String),
    ReadFileError(String, std::io::ErrorKind),
}

impl ProgramError {
    /// Returns `true` if the program error is [`ErrorInOutput`].
    ///
    /// [`ErrorInOutput`]: ProgramError::ErrorInOutput
    #[must_use]
    pub fn is_error_in_output(&self) -> bool {
        matches!(self, Self::ErrorInOutput(..))
    }
}

impl Display for ProgramError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{self:?}")
    }
}

impl Error for ProgramError {}

#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub enum Procedure {
    Opt,
    Freq,
    SinglePt,
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Template {
    pub header: String,
}

impl Template {
    pub fn from(s: &str) -> Self {
        Self {
            header: s.to_string(),
        }
    }
}

impl From<String> for Template {
    fn from(header: String) -> Self {
        Self { header }
    }
}

impl FromStr for Template {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(Self {
            header: s.to_string(),
        })
    }
}

/// A program backend runnable on a [crate::queue::Queue].
///
/// Calculation-specific data lives in [`Job`]; implementations provide the
/// shared behavior for writing inputs and reading outputs.
pub trait Program: Sync + Debug {
    /// Render the shell command that runs the input rooted at `filename`.
    fn command(&self, filename: &str) -> String;

    /// Whether PBS must submit the script from its parent directory.
    ///
    /// Molpro 2022 requires this on the cluster where PBQFF was developed. The
    /// queue uses the same relative path when rendering [`Self::command`].
    fn submit_from_script_dir(&self) -> bool {
        false
    }

    /// Return the output associated with `job`.
    fn outfile(&self, job: &Job) -> String {
        job.filename.clone() + ".out"
    }

    /// Return the input file associated with `job`.
    fn infile(&self, job: &Job) -> String;

    /// the file extension for the input file
    fn extension(&self) -> &'static str;

    /// Write the input file for `job`.
    fn write_input(&self, job: &Job, proc: Procedure);

    /// read the output file `filename`
    fn read_output(
        &self,
        filename: &str,
    ) -> Result<ProgramResult, ProgramError>;

    /// Return all filenames associated with `job` for deletion when it
    /// finishes.
    fn associated_files(&self, job: &Job) -> Vec<String>;

    /// Build the jobs described by `moles` in memory, but don't write any of
    /// their files yet
    #[allow(clippy::too_many_arguments)]
    fn build_jobs(
        &self,
        moles: Vec<Geom>,
        dir: &Path,
        start_index: usize,
        coeff: f64,
        job_num: usize,
        charge: isize,
        tmpl: Template,
    ) -> Vec<Job> {
        let mut count: usize = start_index;
        let mut job_num = job_num;
        let mut jobs = Vec::new();
        for mol in moles {
            let filename = format!("job.{job_num:08}");
            let filename = dir.join(filename).to_str().unwrap().to_string();
            job_num += 1;
            let mut job = Job::new(filename, tmpl.clone(), charge, mol, count);
            job.coeff = coeff;
            jobs.push(job);
            count += 1;
        }
        jobs
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Job {
    /// Filename without the program's input extension.
    pub filename: String,

    /// Template used to write the program input.
    pub template: Template,

    /// Molecular charge.
    pub charge: isize,

    /// Input geometry.
    pub geom: Geom,

    pub pbs_file: String,
    pub job_id: String,

    /// the index in the output array to store the result
    pub index: usize,

    /// the coefficient to multiply by when storing the result
    pub coeff: f64,

    /// the last modified time of the program's output file
    pub(crate) modtime: SystemTime,
}

impl Job {
    pub fn new(
        filename: String,
        template: Template,
        charge: isize,
        geom: Geom,
        index: usize,
    ) -> Self {
        Self {
            filename,
            template,
            charge,
            geom,
            pbs_file: String::new(),
            job_id: String::new(),
            index,
            coeff: 1.0,
            modtime: SystemTime::UNIX_EPOCH,
        }
    }

    /// Return the current modification time of `outfile`, or `self.modtime` if
    /// its metadata cannot be read.
    pub fn modtime(&self, outfile: impl AsRef<Path>) -> SystemTime {
        if let Ok(meta) = std::fs::metadata(outfile) {
            meta.modified().unwrap()
        } else {
            self.modtime
        }
    }
}

/// parses the `nth` field of `line` into a float and returns
/// [ProgramError::EnergyParseError] containing `outname` if it fails. a string
/// containing `outname` is allocated in the Err case
#[inline]
fn parse_energy(
    line: &str,
    nth: usize,
    outname: &str,
) -> Result<Option<f64>, ProgramError> {
    line.split_whitespace()
        .nth(nth)
        .map(str::parse::<f64>)
        .transpose()
        .map_err(|_| ProgramError::EnergyParseError(outname.to_owned()))
}

#[cfg(test)]
mod tests {
    use super::Program;
    use super::mopac::Mopac;

    #[test]
    fn program_is_object_safe() {
        let _: &dyn Program = &Mopac;
    }
}
