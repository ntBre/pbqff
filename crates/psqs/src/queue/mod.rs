use std::{
    cell::LazyCell,
    collections::{HashMap, HashSet},
    fmt::Debug,
    path::Path,
    process::Command,
    str,
    time::Duration,
};

use crate::{
    NO_RESUB,
    geom::Geom,
    program::{Procedure, Program, ProgramError},
};
use crate::{
    program::{Job, ProgramResult},
    time,
};

pub mod local;
pub mod pbs;
pub mod slurm;
pub mod templates;
use drain::*;
mod drain;

pub use drain::Check;

#[derive(PartialEq, Eq, Debug)]
pub struct Resubmit {
    pub inp_file: String,
    pub pbs_file: String,
    pub job_id: String,
}

/// A queue backend that writes, submits, and monitors calculation jobs.
///
/// Program-specific input, output, and launch behavior lives in [`Program`],
/// so adding a program does not require another implementation for every
/// queue.
pub trait Queue: Sync + Debug {
    /// The extension to append to submission scripts for this queue.
    fn script_ext(&self) -> &'static str;

    fn dir(&self) -> &str;

    fn submit_command(&self) -> &str;

    fn chunk_size(&self) -> usize;

    fn job_limit(&self) -> usize;

    fn sleep_int(&self) -> usize;

    /// The command to check the status of jobs in the queue.
    fn stat_cmd(&self) -> String;

    /// Return the jobs found in the queue based on [`Self::stat_cmd`].
    fn status(&self) -> HashSet<String>;

    /// Return `true` if all output files should be preserved.
    fn no_del(&self) -> bool;

    /// The submission-script template selected during configuration.
    fn template(&self) -> &str;

    /// Submit `filename` to the queue and return the job ID.
    fn submit(&self, _program: &dyn Program, filename: &str) -> String {
        loop {
            match Command::new(self.submit_command()).arg(filename).output() {
                Ok(s) => {
                    if s.status.success() {
                        let raw = str::from_utf8(&s.stdout)
                            .unwrap()
                            .trim()
                            .to_string();
                        return raw
                            .split_whitespace()
                            .last()
                            .unwrap_or("")
                            .to_string();
                    }
                    log::warn!(
                        "failed to submit {filename} with `{}`",
                        String::from_utf8_lossy(&s.stderr)
                    );
                    if *NO_RESUB {
                        std::process::exit(1);
                    }
                    std::thread::sleep(Duration::from_secs(1));
                }
                Err(e) => panic!("{e:?}"),
            };
        }
    }

    fn write_submit_script(
        &self,
        program: &dyn Program,
        infiles: &[String],
        filename: &str,
    ) {
        use std::fmt::Write;
        let path = Path::new(filename);
        let basename = path.file_name().unwrap();
        let mut body = self
            .template()
            .to_owned()
            .replace("{{.basename}}", basename.to_str().unwrap())
            .replace("{{.filename}}", filename);
        for f in infiles {
            let f = self.program_filename(program, f);
            writeln!(body, "{}", program.command(&f)).unwrap();
        }
        if std::fs::write(filename, body).is_err() {
            panic!("write_submit_script: failed to create {filename}");
        };
    }

    /// Adjust the input filename used inside a submission script.
    fn program_filename(
        &self,
        _program: &dyn Program,
        filename: &str,
    ) -> String {
        filename.to_owned()
    }

    /// take a name of a Program input file with the extension attached, replace
    /// the extension (ext) with _redo.ext and write _redo.SCRIPT_EXT, then
    /// submit the redo script
    fn resubmit(&self, program: &dyn Program, path: &Path) -> Resubmit {
        let dir = path.parent().unwrap().to_str().unwrap();
        let base = path.file_stem().unwrap().to_str().unwrap();
        {
            let ext = path.extension().unwrap().to_str().unwrap();
            let inp_file = format!("{dir}/{base}_redo.{ext}");
            if let Err(e) = std::fs::copy(path, &inp_file) {
                panic!("failed to copy {path:?} to {inp_file} with `{e}`")
            }
        }
        // nothing but the copy needs the name with extension
        let inp_name = format!("{dir}/{base}_redo");
        let pbs_file = format!("{}/{}_redo.{}", dir, base, self.script_ext());
        self.write_submit_script(
            program,
            std::slice::from_ref(&inp_name),
            &pbs_file,
        );
        let job_id = self.submit(program, &pbs_file);
        Resubmit {
            inp_file: inp_name,
            pbs_file,
            job_id,
        }
    }

    /// Build a chunk of jobs by writing the Program input file and the
    /// corresponding submission script and then submitting the script. returns
    /// the total durations spent writing input files, writing the submit
    /// script, and submitting the script
    fn build_chunk(
        &self,
        program: &dyn Program,
        dir: &str,
        jobs: &mut [Job],
        chunk_num: usize,
        proc: Procedure,
    ) -> (HashMap<String, usize>, Duration, Duration, Duration) {
        self.build_chunk_inner(program, dir, "main", chunk_num, jobs, proc)
    }

    fn build_chunk_inner(
        &self,
        program: &dyn Program,
        dir: &str,
        base: &str,
        chunk_num: usize,
        jobs: &mut [Job],
        proc: Procedure,
    ) -> (HashMap<String, usize>, Duration, Duration, Duration) {
        let mut input = Duration::default();
        let mut script = Duration::default();
        let mut submit = Duration::default();
        let queue_file =
            format!("{}/{base}{}.{}", dir, chunk_num, self.script_ext());
        let jl = jobs.len();
        let mut slurm_jobs = HashMap::new();
        let filenames: Vec<_> = jobs
            .iter_mut()
            .map(|job| {
                time!(e, {
                    program.write_input(job, proc);
                });
                input += e;
                job.pbs_file = queue_file.to_string();
                job.filename.clone()
            })
            .collect();
        slurm_jobs.insert(queue_file.clone(), jl);
        time!(e, {
            self.write_submit_script(program, &filenames, &queue_file);
        });
        script += e;
        // run jobs
        let job_id;
        time!(e, {
            job_id = self.submit(program, &queue_file);
        });
        submit += e;
        for job in jobs {
            job.job_id = job_id.clone();
        }
        (slurm_jobs, input, script, submit)
    }

    fn drain_err_case(
        &self,
        program: &dyn Program,
        e: ProgramError,
        qstat: &mut HashSet<String>,
        slurm_jobs: &mut HashMap<String, usize>,
        job: &mut Job,
    ) {
        let no_resub = LazyCell::new(|| std::env::var("SEMP_RESUB").is_ok());
        // just overwrite the existing job with the resubmitted
        // version
        if !qstat.contains(&job.job_id) {
            let time = job.modtime(program.outfile(job));
            if time > job.modtime {
                // file has been updated since we last looked at it, so need to
                // look again
                job.modtime = time;
                return;
            }
            eprintln!(
                "resubmitting {} (id={}) for {:?}",
                job.filename, job.job_id, e
            );
            if *no_resub {
                eprintln!(
                    "resubmission disabled by SEMP_RESUB environment variable, exiting"
                );
                std::process::exit(1);
            }
            let resub = format!("{}.{}", job.filename, program.extension());
            let Resubmit {
                inp_file,
                pbs_file,
                job_id,
            } = self.resubmit(program, Path::new(&resub));
            job.filename = inp_file;
            job.pbs_file = pbs_file.clone();
            slurm_jobs.insert(pbs_file, 1);
            qstat.insert(job_id.clone());
            job.job_id = job_id;
        }
    }

    /// optimize is a copy of drain for optimizing jobs
    fn optimize(
        &self,
        program: &dyn Program,
        dir: &str,
        jobs: Vec<Job>,
        dst: &mut [Geom],
    ) -> Result<f64, Vec<usize>> {
        Opt.drain(program, dir, self, jobs, dst, Check::None)
    }

    /// resume draining from the checkpoint file in `checkpoint`
    fn resume(
        &self,
        program: &dyn Program,
        dir: &str,
        checkpoint: &str,
        dst: &mut [f64],
        check: Check,
    ) -> Result<f64, Vec<usize>> {
        let jobs = Single::load_checkpoint(checkpoint, dst);
        eprintln!(
            "resuming from checkpoint in '{checkpoint}' with {} jobs remaining",
            jobs.len()
        );
        self.drain(program, dir, jobs, dst, check)
    }

    /// run the single-point energy calculations in `jobs`, storing the results
    /// in `dst`. if `check_int` > 0, write checkpoint files at that interval
    fn drain(
        &self,
        program: &dyn Program,
        dir: &str,
        jobs: Vec<Job>,
        dst: &mut [f64],
        check: Check,
    ) -> Result<f64, Vec<usize>> {
        Single.drain(program, dir, self, jobs, dst, check)
    }

    fn energize(
        &self,
        program: &dyn Program,
        dir: &str,
        jobs: Vec<Job>,
        dst: &mut [ProgramResult],
    ) -> Result<f64, Vec<usize>> {
        Both.drain(program, dir, self, jobs, dst, Check::None)
    }
}

#[cfg(test)]
mod tests {
    use super::{Queue, local::Local};

    #[test]
    fn queue_is_object_safe() {
        let _: &dyn Queue = &Local::default();
    }
}
