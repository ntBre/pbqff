use std::path::Path;
use std::time::Duration;
use std::{collections::HashSet, process::Command};

use crate::program::Program;
use crate::queue::Queue;

/// Pbs is a type for holding the information for submitting a pbs job.
/// `filename` is the name of the Pbs submission script
#[derive(Debug)]
pub struct Pbs {
    pub chunk_size: usize,
    pub job_limit: usize,
    pub sleep_int: usize,
    pub dir: &'static str,
    pub no_del: bool,
    pub template: String,
}

impl Pbs {
    pub fn new(
        chunk_size: usize,
        job_limit: usize,
        sleep_int: usize,
        dir: &'static str,
        no_del: bool,
        template: String,
    ) -> Self {
        Self {
            chunk_size,
            job_limit,
            sleep_int,
            dir,
            no_del,
            template,
        }
    }
}

/// Submit a PBS command with bounded retries.
fn submit_inner(
    cmd: &mut Command,
    sleep_int: usize,
) -> std::io::Result<String> {
    let mut retries = 5;
    loop {
        match cmd.output() {
            Ok(s) => {
                if !s.status.success() {
                    if retries > 0 {
                        eprintln!(
                            "qsub failed with output: {s:#?}, \
				   retrying {retries} more times"
                        );
                        retries -= 1;
                        std::thread::sleep(Duration::from_secs(
                            sleep_int as u64,
                        ));
                        continue;
                    }
                    panic!("qsub failed with output: {s:#?}");
                }
                let raw =
                    std::str::from_utf8(&s.stdout).unwrap().trim().to_string();
                return Ok(raw
                    .split_whitespace()
                    .last()
                    .unwrap_or("no jobid")
                    .to_string());
            }
            Err(e) => return Err(e),
        }
    }
}

impl Queue for Pbs {
    fn script_ext(&self) -> &'static str {
        "pbs"
    }

    fn template(&self) -> &str {
        &self.template
    }

    /// Molpro 2022 requires submission from the script's directory. Other
    /// programs use the normal `qsub path/to/script` form.
    fn submit(&self, program: &dyn Program, filename: &str) -> String {
        let mut cmd = Command::new(self.submit_command());
        if program.submit_from_script_dir() {
            let path = Path::new(filename);
            let dir = path.parent().unwrap();
            let base = path.file_name().unwrap();
            cmd.arg(base).current_dir(dir);
        } else {
            cmd.arg(filename);
        }
        submit_inner(&mut cmd, self.sleep_int).unwrap()
    }

    fn program_filename(
        &self,
        program: &dyn Program,
        filename: &str,
    ) -> String {
        if program.submit_from_script_dir() {
            Path::new(filename)
                .file_name()
                .unwrap()
                .to_string_lossy()
                .into_owned()
        } else {
            filename.to_owned()
        }
    }

    fn submit_command(&self) -> &str {
        "qsub"
    }

    fn chunk_size(&self) -> usize {
        self.chunk_size
    }

    fn job_limit(&self) -> usize {
        self.job_limit
    }

    fn sleep_int(&self) -> usize {
        self.sleep_int
    }

    fn dir(&self) -> &str {
        self.dir
    }

    /// run `qstat -u $USER`. form of the output is:
    ///
    /// maple:
    ///                                                     Req'd  Req'd   Elap
    /// Job ID  Username Queue    Jobname    SessID NDS TSK Memory Time  S Time
    /// ------- -------- -------- ---------- ------ --- --- ------ ----- - -----
    /// 819446  user     queue    C6HNpts      5085   1   1    8gb 26784 R 00:00
    fn stat_cmd(&self) -> String {
        let user = std::env::var("USER").expect("couldn't find $USER env var");
        let status = match Command::new("qstat").args(["-u", &user]).output() {
            Ok(status) => status,
            Err(e) => panic!("failed to run `qstat -u {user}` with {e}"),
        };
        String::from_utf8(status.stdout).expect("failed to parse qstat output")
    }

    fn status(&self) -> HashSet<String> {
        let mut ret = HashSet::new();
        let lines = self.stat_cmd();
        // skip to end of header
        let lines = lines.lines().skip_while(|l| !l.contains("-----------"));
        for line in lines {
            let fields: Vec<_> = line.split_whitespace().collect();
            assert!(fields.len() == 11);
            ret.insert(fields[0].to_string());
        }
        ret
    }

    fn no_del(&self) -> bool {
        self.no_del
    }
}

#[cfg(test)]
mod tests {
    use insta::assert_snapshot;

    use crate::{
        program::{
            cfour::Cfour, dftbplus::DFTBPlus, molpro::Molpro, mopac::Mopac,
            orca::Orca,
        },
        queue::templates,
    };

    use super::*;

    fn pbs(template: &str) -> Pbs {
        Pbs {
            chunk_size: 1,
            job_limit: 1,
            sleep_int: 1,
            dir: "/tmp",
            no_del: false,
            template: template.to_owned(),
        }
    }

    macro_rules! make_tests {
        ($($name:ident, $queue:expr => $program:expr$(,)*)*) => {
            $(
            #[test]
            fn $name() {
                let tmp = tempfile::NamedTempFile::new().unwrap();
                Queue::write_submit_script(
                    $queue,
                    &$program,
                    &["pts/opt0", "pts/opt1", "pts/opt2", "pts/opt3"]
                    .map(str::to_owned),
                    tmp.path().to_str().unwrap(),
                );
                let got = std::fs::read_to_string(tmp).unwrap();
                let got: Vec<&str> = got
                    .lines()
                    .filter(|l| {
                        !(l.starts_with("#PBS -N")
                            || l.starts_with("#PBS -o"))
                    })
                    .collect();
                let got = got.join("\n");
                assert_snapshot!(got);
            }
            )*
        }
    }

    make_tests! {
        mopac_pbs, &pbs(templates::PBS_MOPAC) => Mopac,
        molpro_pbs, &pbs(templates::PBS_MOLPRO) => Molpro,
        cfour_pbs, &pbs(templates::PBS_CFOUR) => Cfour,
        dftb_pbs, &pbs(templates::PBS_DFTBPLUS) => DFTBPlus,
        orca_pbs, &pbs(templates::PBS_ORCA) => Orca,
    }
}
