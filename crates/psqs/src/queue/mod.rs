use std::{
    cell::LazyCell,
    collections::{HashMap, HashSet},
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
use drain::*;
use serde::{Deserialize, Serialize};
mod drain;

pub use drain::Check;

#[derive(PartialEq, Eq, Debug)]
pub struct Resubmit {
    pub inp_file: String,
    pub pbs_file: String,
    pub job_id: String,
}

pub trait Submit<P>: SubQueue<P>
where
    P: Program + Clone + Serialize + for<'a> Deserialize<'a>,
{
    /// submit `filename` to the queue and return the jobid
    fn submit(&self, filename: &str) -> String {
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
}

/// a trait for all of the program-independent parts of a [Queue]
pub trait SubQueue<P>
where
    P: Program + Clone + Serialize + for<'a> Deserialize<'a>,
{
    /// the extension to append to submit scripts for this type of Queue
    const SCRIPT_EXT: &'static str;

    fn dir(&self) -> &str;

    fn submit_command(&self) -> &str;

    fn chunk_size(&self) -> usize;

    fn job_limit(&self) -> usize;

    fn sleep_int(&self) -> usize;

    /// the command to check the status of jobs in the queue
    fn stat_cmd(&self) -> String;

    /// return a HashSet of jobs found in the queue based on the output of
    /// `stat_cmd`
    fn status(&self) -> HashSet<String>;

    /// return `true` if all output files should be preserved
    fn no_del(&self) -> bool;
}

pub trait Queue<P>: SubQueue<P> + Submit<P>
where
    P: Program + Clone + Send + Sync + Serialize + for<'a> Deserialize<'a>,
{
    fn default_submit_script(&self) -> String;

    fn template(&self) -> &Option<String>;

    fn program_cmd(&self, filename: &str) -> String;

    fn write_submit_script(
        &self,
        infiles: impl IntoIterator<Item = String>,
        filename: &str,
    ) {
        use std::fmt::Write;
        let path = Path::new(filename);
        let basename = path.file_name().unwrap();
        let mut body = self
            .template()
            .clone()
            .unwrap_or_else(|| <Self as Queue<P>>::default_submit_script(self))
            .replace("{{.basename}}", basename.to_str().unwrap())
            .replace("{{.filename}}", filename);
        for f in infiles {
            writeln!(body, "{}", self.program_cmd(&f)).unwrap();
        }
        if std::fs::write(filename, body).is_err() {
            panic!("write_submit_script: failed to create {filename}");
        };
    }

    /// take a name of a Program input file with the extension attached, replace
    /// the extension (ext) with _redo.ext and write _redo.SCRIPT_EXT, then
    /// submit the redo script
    fn resubmit(&self, path: impl AsRef<Path>) -> Resubmit {
        let path = path.as_ref();
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
        let pbs_file = format!("{}/{}_redo.{}", dir, base, Self::SCRIPT_EXT);
        self.write_submit_script([inp_name.clone()], &pbs_file);
        let job_id = self.submit(&pbs_file);
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
        dir: &str,
        jobs: &mut [Job<P>],
        chunk_num: usize,
        proc: Procedure,
    ) -> (HashMap<String, usize>, Duration, Duration, Duration) {
        self.build_chunk_inner(dir, "main", chunk_num, jobs, proc)
    }

    fn build_chunk_inner(
        &self,
        dir: &str,
        base: &str,
        chunk_num: usize,
        jobs: &mut [Job<P>],
        proc: Procedure,
    ) -> (HashMap<String, usize>, Duration, Duration, Duration) {
        let mut input = Duration::default();
        let mut script = Duration::default();
        let mut submit = Duration::default();
        let queue_file =
            format!("{}/{base}{}.{}", dir, chunk_num, Self::SCRIPT_EXT);
        let jl = jobs.len();
        let mut slurm_jobs = HashMap::new();
        let filenames = jobs.iter_mut().map(|job| {
            time!(e, {
                job.program.write_input(proc);
            });
            input += e;
            job.pbs_file = queue_file.to_string();
            job.program.filename()
        });
        slurm_jobs.insert(queue_file.clone(), jl);
        time!(e, {
            self.write_submit_script(filenames, &queue_file);
        });
        script += e;
        // run jobs
        let job_id;
        time!(e, {
            job_id = self.submit(&queue_file);
        });
        submit += e;
        for job in jobs {
            job.job_id = job_id.clone();
        }
        (slurm_jobs, input, script, submit)
    }

    fn drain_err_case(
        &self,
        e: ProgramError,
        qstat: &mut HashSet<String>,
        slurm_jobs: &mut HashMap<String, usize>,
        job: &mut Job<P>,
    ) {
        let no_resub = LazyCell::new(|| std::env::var("SEMP_RESUB").is_ok());
        // just overwrite the existing job with the resubmitted
        // version
        if !qstat.contains(&job.job_id) {
            let time = job.modtime();
            if time > job.modtime {
                // file has been updated since we last looked at it, so need to
                // look again
                job.modtime = time;
                return;
            }
            eprintln!(
                "resubmitting {} (id={}) for {:?}",
                job.program.filename(),
                job.job_id,
                e
            );
            if *no_resub {
                eprintln!(
                    "resubmission disabled by SEMP_RESUB environment variable, exiting"
                );
                std::process::exit(1);
            }
            let resub = format!(
                "{}.{}",
                job.program.filename(),
                job.program.extension()
            );
            let Resubmit {
                inp_file,
                pbs_file,
                job_id,
            } = self.resubmit(&resub);
            job.program.set_filename(&inp_file);
            job.pbs_file = pbs_file.clone();
            slurm_jobs.insert(pbs_file, 1);
            qstat.insert(job_id.clone());
            job.job_id = job_id;
        }
    }

    /// optimize is a copy of drain for optimizing jobs
    fn optimize(
        &self,
        dir: &str,
        jobs: Vec<Job<P>>,
        dst: &mut [Geom],
    ) -> Result<f64, Vec<usize>>
    where
        Self: Sync,
    {
        Opt.drain(dir, self, jobs, dst, Check::None)
    }

    /// resume draining from the checkpoint file in `checkpoint`
    fn resume(
        &self,
        dir: &str,
        checkpoint: &str,
        dst: &mut [f64],
        check: Check,
    ) -> Result<f64, Vec<usize>>
    where
        Self: Sync,
    {
        let jobs = Single::load_checkpoint(checkpoint, dst);
        eprintln!(
            "resuming from checkpoint in '{checkpoint}' with {} jobs remaining",
            jobs.len()
        );
        self.drain(dir, jobs, dst, check)
    }

    /// run the single-point energy calculations in `jobs`, storing the results
    /// in `dst`. if `check_int` > 0, write checkpoint files at that interval
    fn drain(
        &self,
        dir: &str,
        jobs: Vec<Job<P>>,
        dst: &mut [f64],
        check: Check,
    ) -> Result<f64, Vec<usize>>
    where
        Self: Sync,
    {
        Single.drain(dir, self, jobs, dst, check)
    }

    fn energize(
        &self,
        dir: &str,
        jobs: Vec<Job<P>>,
        dst: &mut [ProgramResult],
    ) -> Result<f64, Vec<usize>>
    where
        Self: Sync,
    {
        Both.drain(dir, self, jobs, dst, Check::None)
    }
}
