use crate::geom::geom_string;
use crate::program::{Program, ProgramError};
use regex::Regex;
use symm::Atom;

use super::{Job, Procedure, ProgramResult};
use std::fs::{File, read_to_string};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::sync::OnceLock;

/// kcal/mol per hartree
pub const KCALHT: f64 = 627.5091809;

pub use self::params::*;
pub mod params;

#[cfg(test)]
mod tests;

/// Shared MOPAC input/output backend.
#[derive(Debug, Clone, Copy, Default)]
pub struct Mopac;

impl Program for Mopac {
    fn command(&self, filename: &str) -> String {
        format!("$MOPAC_CMD {filename}.mop")
    }

    fn extension(&self) -> &'static str {
        "mop"
    }

    fn write_input(&self, job: &Job, proc: Procedure) {
        use std::fmt::Write;
        // header should look like
        //   scfcrt=1.D-21 aux(precision=14) PM6
        // so that the charge, and optionally XYZ, A0, and 1SCF can be added
        let mut header = job.template.clone().header;
        write!(header, " charge={}", job.charge).unwrap();
        match proc {
            Procedure::Opt => {
                // optimization is the default, so just don't add 1SCF
            }
            Procedure::Freq => todo!(),
            Procedure::SinglePt => {
                header.push_str(" 1SCF");
            }
        }
        if job.geom.is_xyz() {
            header.push_str(" XYZ");
        }
        let geom = geom_string(&job.geom);
        let filename = format!("{}.mop", job.filename);
        let mut file = match File::create(&filename) {
            Ok(f) => f,
            Err(e) => panic!("failed to create {filename} with {e}"),
        };
        write!(
            file,
            "{header}
Comment line 1
Comment line 2
{geom}
",
        )
        .expect("failed to write input file");
    }

    /// Reads a MOPAC output file. If normal termination occurs, also try
    /// reading the `.aux` file to extract the energy from there. This function
    /// panics if an error is found in the output file. If a non-fatal error
    /// occurs (file not found, not written to yet, etc) None is returned.
    fn read_output(
        &self,
        filename: &str,
    ) -> Result<ProgramResult, ProgramError> {
        let res = Self::read_aux(filename);
        if res.is_ok() {
            return res;
        }
        let outfile = format!("{}.out", &filename);
        let contents = match read_to_string(&outfile) {
            Ok(s) => s,
            Err(_) => {
                return Err(ProgramError::FileNotFound(outfile));
            }
        };

        let [panic, error] = READ_OUT_CELL.get_or_init(|| {
            [
                Regex::new("(?i)panic").unwrap(),
                Regex::new("(?i)error").unwrap(),
            ]
        });

        if error.is_match(&contents) {
            return Err(ProgramError::ErrorInOutput(filename.to_owned()));
        } else if panic.is_match(&contents) {
            panic!("panic requested in read_output");
        }
        res
    }

    fn associated_files(&self, job: &Job) -> Vec<String> {
        let fname = &job.filename;
        vec![
            format!("{fname}.mop"),
            format!("{fname}.out"),
            format!("{fname}.arc"),
            format!("{fname}.aux"),
        ]
    }

    fn infile(&self, job: &Job) -> String {
        job.filename.clone() + ".mop"
    }
}

static READ_OUT_CELL: OnceLock<[Regex; 2]> = OnceLock::new();
static READ_AUX_CELL: OnceLock<[Regex; 6]> = OnceLock::new();

impl Mopac {
    /// write the `params` to `filename`
    pub fn write_params(params: &Params, path: impl AsRef<Path>) {
        let path = path.as_ref();
        let mut file = match File::create(path) {
            Ok(f) => f,
            Err(e) => {
                eprintln!("failed to create {path:?} with {e}");
                std::process::exit(1);
            }
        };
        write!(file, "{params}").expect("failed to write params file");
    }

    /// return the heat of formation from a MOPAC aux file in Hartrees.
    /// `filename` should not include the .aux extension
    pub fn read_aux(filename: &str) -> Result<ProgramResult, ProgramError> {
        let auxfile = format!("{}.aux", &filename);
        let Ok(f) = File::open(&auxfile) else {
            return Err(ProgramError::FileNotFound(auxfile));
        };
        let mut energy = None;

        let [heat_re, atom_re, elt_re, core_re, charge_re, time_re] =
            READ_AUX_CELL.get_or_init(|| {
                [
                    Regex::new("^ HEAT_OF_FORMATION").unwrap(),
                    Regex::new("^ ATOM_X_OPT").unwrap(),
                    Regex::new("^ ATOM_EL").unwrap(),
                    Regex::new("^ ATOM_CORE").unwrap(),
                    Regex::new("^ ATOM_CHARGES").unwrap(),
                    Regex::new("^ CPU_TIME:SEC=").unwrap(),
                ]
            });
        #[derive(PartialEq)]
        enum State {
            Geom,
            Labels,
            Done,
            None,
        }
        /// don't look for these after they've been found
        struct Guard {
            heat: bool,
            atom: bool,
            element: bool,
            time: bool,
        }
        let mut state = State::None;
        let mut guard = Guard {
            heat: false,
            atom: false,
            element: false,
            time: false,
        };
        // atomic labels
        let mut labels = Vec::new();
        // coordinates
        let mut coords = Vec::new();
        let mut time = 0.0;
        for line in BufReader::new(f).lines().map_while(Result::ok) {
            if !guard.element && elt_re.is_match(&line) {
                state = State::Labels;
                guard.element = true;
            } else if state == State::Labels && core_re.is_match(&line) {
                state = State::None;
            } else if state == State::Labels {
                labels
                    .extend(line.split_ascii_whitespace().map(str::to_string));
            // line like HEAT_OF_FORMATION:KCAL/MOL=+0.97127947459164715838D+02
            } else if !guard.heat && heat_re.is_match(&line) {
                let fields: Vec<&str> = line.trim().split('=').collect();
                match fields[1].replace('D', "E").parse::<f64>() {
                    Ok(f) => {
                        energy = Some(f / KCALHT);
                    }
                    Err(_) => {
                        return Err(ProgramError::EnergyParseError(auxfile));
                    }
                }
                guard.heat = true;
            } else if !guard.time && time_re.is_match(&line) {
                time = line
                    .split('=')
                    .nth(1)
                    .unwrap()
                    .replace('D', "E")
                    .parse()
                    .unwrap();
                guard.time = true;
            } else if !guard.atom && atom_re.is_match(&line) {
                state = State::Geom;
                guard.atom = true;
            } else if state == State::Geom && charge_re.is_match(&line) {
                state = State::Done;
                break;
            } else if state == State::Geom {
                coords.extend(
                    line.split_ascii_whitespace()
                        .map(|s| s.parse::<f64>().unwrap()),
                );
            }
        }
        if state != State::Done {
            return Err(ProgramError::GeomNotFound(auxfile));
        }
        assert_eq!(coords.len() / 3, labels.len());
        let ret = coords
            .chunks_exact(3)
            .zip(labels)
            .map(|(coord, l)| {
                Atom::new_from_label(&l, coord[0], coord[1], coord[2])
            })
            .collect();
        if let Some(energy) = energy {
            Ok(ProgramResult {
                energy,
                cart_geom: Some(ret),
                time,
            })
        } else {
            Err(ProgramError::EnergyNotFound(auxfile))
        }
    }
}
