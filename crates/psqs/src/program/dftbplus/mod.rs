use std::{
    fs::{read_to_string, File},
    path::Path,
    sync::OnceLock,
};

use log::{trace, warn};
use regex::Regex;
use serde::{Deserialize, Serialize};
use symm::Atom;

use crate::{geom::Geom, program::Procedure};

use super::{parse_energy, Program, ProgramError, ProgramResult, Template};

#[cfg(test)]
mod tests;

static INPUT_CELL: OnceLock<[Regex; 3]> = OnceLock::new();
static CELL: OnceLock<[Regex; 5]> = OnceLock::new();

#[derive(Clone, Deserialize, Serialize)]
pub struct DFTBPlus {
    /// in this case, `filename` is actually a directory name because every
    /// DFTB+ input file has to have the same name
    filename: String,
    template: Template,
    charge: isize,
    geom: Geom,
}

impl Program for DFTBPlus {
    fn filename(&self) -> String {
        self.filename.clone()
    }

    fn infile(&self) -> String {
        todo!()
    }

    fn set_filename(&mut self, filename: &str) {
        self.filename = filename.into();
    }

    fn template(&self) -> &Template {
        &self.template
    }

    /// every file has to have the same name, so I don't actually need to match
    /// up extensions
    fn extension(&self) -> String {
        String::new()
    }

    fn charge(&self) -> isize {
        self.charge
    }

    /// Example [Template]:
    /// ```text
    /// Geometry = xyzFormat {
    /// {{.geom}}
    /// }
    ///
    /// Hamiltonian = DFTB {
    ///   Scc = Yes
    ///   SlaterKosterFiles = Type2FileNames {
    ///     Prefix = "/opt/dftb+/slako/mio/mio-1-1/"
    ///     Separator = "-"
    ///     Suffix = ".skf"
    ///   }
    ///   MaxAngularMomentum {
    ///     O = "p"
    ///     H = "s"
    ///   }
    ///   Charge = {{.charge}}
    /// }
    ///
    /// Options {}
    ///
    /// Analysis {
    ///   CalculateForces = Yes
    /// }
    ///
    /// ParserOptions {
    ///   ParserVersion = 12
    /// }
    /// ```
    fn write_input(&mut self, proc: Procedure) {
        use std::io::Write;
        let mut body = self.template().clone().header;
        // skip optgrad but accept optg at the end of a line
        let [opt, charge, geom_re] = INPUT_CELL.get_or_init(|| {
            [
                Regex::new(r"(?i)Driver = GeometryOptimization").unwrap(),
                Regex::new(r"\{\{.charge\}\}").unwrap(),
                Regex::new(r"\{\{.geom\}\}").unwrap(),
            ]
        });
        let mut found_opt = false;
        if opt.is_match(&body) {
            found_opt = true;
        }
        {
            use std::fmt::Write;
            match proc {
                Procedure::Opt => {
                    if !found_opt {
                        writeln!(
                            body,
                            r#"Driver = GeometryOptimization {{
  Optimizer = Rational {{}}
  MovedAtoms = 1:-1
  MaxSteps = 100
  OutputPrefix = "geom.out"
  Convergence {{
        Energy = 1e-8
        GradElem = 1e-7
        GradNorm = 1e-7
        DispElem = 1e-7
        DispNorm = 1e-7
}}
}}"#,
                        )
                        .unwrap();
                    }
                }
                Procedure::Freq => todo!(),
                Procedure::SinglePt => {
                    if found_opt {
                        let mut braces: Vec<char> = Vec::new();
                        let mut in_geom = false;
                        let mut new_body = Vec::new();
                        for line in body.lines() {
                            if opt.is_match(line) {
                                in_geom = true;
                                braces.push('{');
                            } else if in_geom {
                                for c in line.chars() {
                                    if c == '{' {
                                        braces.push(c);
                                    } else if c == '}' {
                                        braces.pop();
                                    }
                                }
                            } else {
                                new_body.push(line);
                            }
                            if in_geom && braces.is_empty() {
                                // trailing newline to make tests easier
                                new_body.push("");
                                in_geom = false;
                            }
                        }
                        body = new_body.join("\n");
                    }
                }
            }
        }
        let geom = match &self.geom {
            Geom::Zmat(_) => {
                panic!("don't know how to handle a Z-matrix in dftb+");
            }
            geom @ Geom::Xyz(atoms) => format!("{}\n\n{geom}\n", atoms.len()),
        };
        body = geom_re.replace(&body, geom).to_string();
        body = charge
            .replace(&body, &format!("{}", self.charge))
            .to_string();

        let dir = Path::new(&self.filename);
        std::fs::create_dir_all(dir).unwrap_or_else(|e| {
            panic!("failed to create {} with {e}", self.filename)
        });
        let mut file =
            File::create(dir.join("dftb_in.hsd")).unwrap_or_else(|e| {
                panic!(
                    "failed to create dftb input in {} with {e}",
                    self.filename
                )
            });
        write!(file, "{body}").expect("failed to write input file");
    }

    fn read_output(filename: &str) -> Result<ProgramResult, ProgramError> {
        let path = Path::new(filename);

        let outfile = path.join("out");
        let outname = outfile.to_string_lossy().to_string();
        let contents = match read_to_string(&outfile) {
            Ok(s) => s,
            Err(_) => {
                return Err(ProgramError::FileNotFound(outname));
            }
        };

        let [panic_re, error_re, time_re, energy_re, geom_warn] = CELL
            .get_or_init(|| {
                trace!("initializing dftb+ output regexes");
                [
                    Regex::new("(?i)panic").unwrap(),
                    Regex::new(r"\bERROR\b").unwrap(),
                    Regex::new(r"^Total\s+=\s+").unwrap(),
                    Regex::new(r"^Total Energy: ").unwrap(),
                    Regex::new(r"Geometry did NOT converge!").unwrap(),
                ]
            });

        if panic_re.is_match(&contents) {
            panic!("panic requested in read_output");
        } else if error_re.is_match(&contents) {
            return Err(ProgramError::ErrorInOutput(outname));
        }

        if geom_warn.is_match(&contents) {
            warn!("geometry did not converge, results may be unreliable");
        }

        // main output
        let mut energy = None;
        let mut time = None;
        for line in contents.lines() {
            if time_re.is_match(line) {
                time = Some(
                    line.split_ascii_whitespace()
                        .nth(4)
                        .unwrap()
                        .parse()
                        .unwrap_or_else(|e| panic!("{e:#?}")),
                );
            } else if energy_re.is_match(line) {
                energy = parse_energy(line, 2, &outname)?;
            }
        }

        // read xyz. TODO we only need to do this if it's an optimization
        let geomfile = path.join("geom.out.xyz");
        let cart_geom = if let Ok(s) = std::fs::read_to_string(geomfile) {
            // always a proper XYZ file, so skip n atoms and comment lines
            let mut atoms = Vec::new();
            for line in s.lines().skip(2) {
                let mut sp = line.split_ascii_whitespace();
                // lines look like this with the charge(?) at the end, so take
                // the first four fields:
                // O   0.00000000  -0.71603315   0.00000000    6.59260702
                atoms.push(Atom::new_from_label(
                    sp.next().unwrap(),
                    sp.next().unwrap().parse().unwrap(),
                    sp.next().unwrap().parse().unwrap(),
                    sp.next().unwrap().parse().unwrap(),
                ));
            }
            Some(atoms)
        } else {
            None
        };

        let Some(energy) = energy else {
            return Err(ProgramError::EnergyNotFound(outname));
        };

        let Some(time) = time else {
            // the time is the last thing printed, so don't trust the energy if
            // we don't find the time. we could have read an earlier energy in a
            // geometry optimization, for example
            return Err(ProgramError::EnergyNotFound(outname));
        };

        Ok(ProgramResult {
            energy,
            cart_geom,
            time,
        })
    }

    fn associated_files(&self) -> Vec<String> {
        vec![
            "charges.bin".to_owned(),
            "detailed.out".to_owned(),
            "geom.out.gen".to_owned(),
            "geom.out.xyz".to_owned(),
            "band.out".to_owned(),
            "dftb_pin.hsd".to_owned(),
            "dftb_in.hsd".to_owned(),
        ]
    }

    fn new(
        filename: String,
        template: Template,
        charge: isize,
        geom: Geom,
    ) -> Self {
        Self {
            filename,
            template,
            charge,
            geom,
        }
    }
}
