use std::{
    error::Error, fmt::Display, path::Path, str::FromStr, time::SystemTime,
};

use serde::{Deserialize, Serialize};
use symm::Atom;

use crate::geom::Geom;

pub mod cfour;
pub mod dftbplus;
pub mod molpro;
pub mod mopac;

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

/// A trait for describing programs runnable on a [crate::queue::Queue]
pub trait Program {
    /// returns the file associated with the program's input. it should not
    /// include an extension
    fn filename(&self) -> String;

    /// return the output of `self.filename()` with ".out" appended
    fn outfile(&self) -> String {
        self.filename() + ".out"
    }

    /// return the input file associated with `self`
    fn infile(&self) -> String;

    /// set `filename`
    fn set_filename(&mut self, filename: &str);

    /// the template for writing input files
    fn template(&self) -> &Template;

    /// the file extension for the input file
    fn extension(&self) -> String;

    /// molecular charge
    fn charge(&self) -> isize;

    /// write the input file to the name returned by `filename`
    fn write_input(&mut self, proc: Procedure);

    /// read the output file `filename`
    fn read_output(filename: &str) -> Result<ProgramResult, ProgramError>;

    /// Return all the filenames associated with the Program for deletion when
    /// it finishes
    fn associated_files(&self) -> Vec<String>;

    fn new(
        filename: String,
        template: Template,
        charge: isize,
        geom: Geom,
    ) -> Self;

    /// Build the jobs described by `moles` in memory, but don't write any of
    /// their files yet
    fn build_jobs(
        moles: Vec<Geom>,
        dir: impl AsRef<Path>,
        start_index: usize,
        coeff: f64,
        job_num: usize,
        charge: isize,
        tmpl: Template,
    ) -> Vec<Job<Self>>
    where
        Self: std::marker::Sized,
    {
        let mut count: usize = start_index;
        let mut job_num = job_num;
        let mut jobs = Vec::new();
        for mol in moles {
            let filename = format!("job.{job_num:08}");
            let filename =
                dir.as_ref().join(filename).to_str().unwrap().to_string();
            job_num += 1;
            let mut job = Job::new(
                Self::new(filename, tmpl.clone(), charge, mol.clone()),
                count,
            );
            job.coeff = coeff;
            jobs.push(job);
            count += 1;
        }
        jobs
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Job<P: Program> {
    pub program: P,
    pub pbs_file: String,
    pub job_id: String,

    /// the index in the output array to store the result
    pub index: usize,

    /// the coefficient to multiply by when storing the result
    pub coeff: f64,

    /// the last modified time of `program`'s output file
    pub(crate) modtime: SystemTime,
}

impl<P: Program> Job<P> {
    pub fn new(program: P, index: usize) -> Self {
        Self {
            program,
            pbs_file: String::new(),
            job_id: String::new(),
            index,
            coeff: 1.0,
            modtime: SystemTime::UNIX_EPOCH,
        }
    }

    /// return the current modtime of `self.program`'s output file, or
    /// `self.modtime` if there is an error accessing the metadata
    pub fn modtime(&self) -> SystemTime {
        let p = self.program.outfile();
        if let Ok(meta) = std::fs::metadata(p) {
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
