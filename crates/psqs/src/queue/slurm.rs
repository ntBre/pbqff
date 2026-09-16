use std::collections::HashSet;

use crate::queue::Queue;

/// Slurm is a type for holding the information for submitting a slurm job.
/// `filename` is the name of the Slurm submission script
#[derive(Debug)]
pub struct Slurm {
    chunk_size: usize,
    job_limit: usize,
    sleep_int: usize,
    dir: &'static str,
    no_del: bool,
    pub(crate) template: String,
}

impl Slurm {
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

impl Queue for Slurm {
    fn script_ext(&self) -> &'static str {
        "slurm"
    }

    fn template(&self) -> &str {
        &self.template
    }

    fn submit_command(&self) -> &str {
        "sbatch"
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

    /// Run `squeue -u $USER --format "%.18i %.2t" --noheader` and return the
    /// output.
    ///
    /// The form of the output is:
    ///
    /// ```text
    /// 30627992  R
    /// ```
    ///
    /// where the first field is the JobID, and the second field is its state.
    ///
    /// See <https://man.archlinux.org/man/squeue.1.en> for other format
    /// options.
    fn stat_cmd(&self) -> String {
        let user = std::env::var("USER").expect("couldn't find $USER env var");
        let status = match std::process::Command::new("squeue")
            .args(["--user", &user])
            .args(["--format", "%.18i %.2t"])
            .arg("--noheader")
            .output()
        {
            Ok(status) => status,
            Err(e) => panic!("failed to run squeue with {e}"),
        };
        String::from_utf8(status.stdout)
            .expect("failed to convert squeue output to String")
    }

    fn status(&self) -> HashSet<String> {
        let mut ret = HashSet::new();
        let lines = self.stat_cmd();
        for line in lines.lines() {
            let fields: Vec<_> = line.split_whitespace().collect();
            let [job_id, state] = fields.as_slice() else {
                panic!("unexpected line in squeue output: `{line}`");
            };
            // exclude completing jobs to combat stuck completing bug
            if *state != "CG" {
                ret.insert(job_id.to_string());
            }
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

    fn slurm(template: &str) -> Slurm {
        Slurm {
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
                    &["opt0", "opt1", "opt2", "opt3"].map(str::to_owned),
                    tmp.path().to_str().unwrap(),
                );
                let got = std::fs::read_to_string(tmp).unwrap();
                let got: Vec<&str> = got
                    .lines()
                    .filter(|l| {
                        !(l.starts_with("#SBATCH --job-name")
                            || l.starts_with("#SBATCH -o"))
                    })
                    .collect();
                let got = got.join("\n");
                assert_snapshot!(got);
            }
            )*
        }
    }

    make_tests! {
        mopac_slurm, &slurm(templates::SLURM_MOPAC) => Mopac,
        molpro_slurm, &slurm(templates::SLURM_MOLPRO) => Molpro,
        cfour_slurm, &slurm(templates::SLURM_CFOUR) => Cfour,
        dftb_slurm, &slurm(templates::SLURM_DFTBPLUS) => DFTBPlus,
        orca_slurm, &slurm(templates::SLURM_ORCA) => Orca,
    }
}
