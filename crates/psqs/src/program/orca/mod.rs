use std::{
    fs::{File, read_to_string},
    io::Write,
    sync::LazyLock,
};

use regex::Regex;
use symm::Atom;

use crate::geom::{Geom, geom_string};

use super::{Job, Procedure, Program, ProgramError, ProgramResult};

#[cfg(test)]
mod tests;

/// ORCA input and output support.
#[derive(Debug, Clone, Copy, Default)]
pub struct Orca;

impl Program for Orca {
    fn command(&self, filename: &str) -> String {
        format!("$ORCA_CMD {filename}.inp > {filename}.out 2>&1")
    }

    fn infile(&self, job: &Job) -> String {
        job.filename.clone() + ".inp"
    }

    fn extension(&self) -> &'static str {
        "inp"
    }

    /// Write an ORCA input file.
    ///
    /// The template must contain `{{.procedure}}`, which expands to `Opt`,
    /// `Freq`, or an empty string for a single-point calculation. It may also
    /// contain `{{.geom}}` and `{{.charge}}`. For example:
    ///
    /// ```text
    /// ! B3LYP def2-SVP VeryTightSCF {{.procedure}}
    /// * xyz {{.charge}} 1
    /// {{.geom}}
    /// *
    /// ```
    fn write_input(&self, job: &Job, proc: Procedure) {
        let procedure = match proc {
            Procedure::Opt => "Opt",
            Procedure::Freq => "Freq",
            Procedure::SinglePt => "",
        };
        let mut body = job.template.header.clone();
        assert!(
            body.contains("{{.procedure}}"),
            "ORCA templates must contain {{{{.procedure}}}} so PBQFF can \
             select optimization and single-point calculations"
        );
        let geom = match &job.geom {
            Geom::Xyz(_) => geom_string(&job.geom),
            Geom::Zmat(_) => {
                panic!("PBQFF does not support ORCA Z-matrix input")
            }
        };
        body = body
            .replace("{{.procedure}}", procedure)
            .replace("{{.charge}}", &job.charge.to_string())
            .replace("{{.geom}}", &geom);

        let filename = self.infile(job);
        let mut file = File::create(&filename)
            .unwrap_or_else(|e| panic!("failed to create {filename} with {e}"));
        write!(file, "{body}").expect("failed to write ORCA input file");
    }

    /// Read the last `FINAL SINGLE POINT ENERGY` from a completed output.
    /// An explicitly printed `PBQFF = <energy>` value takes precedence, which
    /// allows ORCA Compound workflows to select a derived energy.
    fn read_output(
        &self,
        filename: &str,
    ) -> Result<ProgramResult, ProgramError> {
        let outfile = format!("{filename}.out");
        let contents = match read_to_string(&outfile) {
            Ok(contents) => contents,
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => {
                return Err(ProgramError::FileNotFound(outfile));
            }
            Err(e) => {
                return Err(ProgramError::ReadFileError(outfile, e.kind()));
            }
        };

        parse_output(&contents, &outfile)
    }

    fn associated_files(&self, job: &Job) -> Vec<String> {
        let filename = &job.filename;
        [
            ".inp",
            ".out",
            ".bibtex",
            ".engrad",
            ".gbw",
            ".hess",
            ".opt",
            ".property.txt",
            ".vibspectrum",
            ".xtbrestart",
            ".xyz",
            "_trj.xyz",
        ]
        .map(|extension| format!("{filename}{extension}"))
        .to_vec()
    }
}

fn parse_output(
    contents: &str,
    outfile: &str,
) -> Result<ProgramResult, ProgramError> {
    static ERROR_REGEXES: LazyLock<[Regex; 2]> = LazyLock::new(|| {
        [
        Regex::new("(?i)panic").unwrap(),
        Regex::new(
            r"(?i)input error|error \(orca_main\)|orca finished by error termination|optimization did not converge",
        )
        .unwrap(),
    ]
    });

    let [panic_re, error_re] = &*ERROR_REGEXES;
    if panic_re.is_match(contents) {
        panic!("panic requested in read_output");
    }
    if error_re.is_match(contents) {
        return Err(ProgramError::ErrorInOutput(outfile.to_owned()));
    }

    let mut energy = None;
    let mut pbqff_energy = None;
    let mut time = None;
    let mut terminated = false;
    let mut geom_state = GeomState::None;
    let mut atoms = Vec::new();

    for line in contents.lines() {
        let trimmed = line.trim();
        if trimmed == "CARTESIAN COORDINATES (ANGSTROEM)" {
            atoms.clear();
            geom_state = GeomState::Header;
            continue;
        }

        match geom_state {
            GeomState::Header => {
                if trimmed.starts_with('-') {
                    geom_state = GeomState::Atoms;
                }
                continue;
            }
            GeomState::Atoms if trimmed.is_empty() => {
                geom_state = GeomState::None;
                continue;
            }
            GeomState::Atoms => {
                let mut fields = line.split_ascii_whitespace();
                let (Some(label), Some(x), Some(y), Some(z)) = (
                    fields.next(),
                    fields.next(),
                    fields.next(),
                    fields.next(),
                ) else {
                    atoms.clear();
                    geom_state = GeomState::None;
                    continue;
                };
                let (Ok(x), Ok(y), Ok(z)) = (x.parse(), y.parse(), z.parse())
                else {
                    atoms.clear();
                    geom_state = GeomState::None;
                    continue;
                };
                atoms.push(Atom::new_from_label(label, x, y, z));
                continue;
            }
            GeomState::None => {}
        }

        if line.trim_start().starts_with("PBQFF =") {
            pbqff_energy = parse_last_float(line, outfile)?;
        } else if line.trim_start().starts_with("FINAL SINGLE POINT ENERGY") {
            energy = parse_last_float(line, outfile)?;
        } else if line.trim_start().starts_with("TOTAL RUN TIME:") {
            time = parse_run_time(line);
        } else if line.contains("****ORCA TERMINATED NORMALLY****") {
            terminated = true;
        }
    }

    let energy = pbqff_energy.or(energy);
    let (Some(energy), Some(time)) = (energy, time) else {
        return Err(ProgramError::EnergyNotFound(outfile.to_owned()));
    };
    if !terminated {
        return Err(ProgramError::EnergyNotFound(outfile.to_owned()));
    }

    Ok(ProgramResult {
        energy,
        cart_geom: (!atoms.is_empty()).then_some(atoms),
        time,
    })
}

fn parse_last_float(
    line: &str,
    outfile: &str,
) -> Result<Option<f64>, ProgramError> {
    line.split_ascii_whitespace()
        .next_back()
        .map(str::parse)
        .transpose()
        .map_err(|_| ProgramError::EnergyParseError(outfile.to_owned()))
}

#[derive(Clone, Copy)]
enum GeomState {
    None,
    Header,
    Atoms,
}

fn parse_run_time(line: &str) -> Option<f64> {
    let fields: Vec<_> = line.split_ascii_whitespace().skip(3).collect();
    let mut total = 0.0;
    for pair in fields.chunks_exact(2) {
        let value: f64 = pair[0].parse().ok()?;
        total += match pair[1] {
            "days" => value * 86_400.0,
            "hours" => value * 3_600.0,
            "minutes" => value * 60.0,
            "seconds" => value,
            "msec" => value / 1_000.0,
            _ => return None,
        };
    }
    Some(total)
}
