use std::{collections::HashMap, path::Path, time::Duration};

use serde::{Deserialize, Serialize};

use crate::{
    program::{Job, Procedure, Program},
    queue::Queue,
};

pub(crate) struct Resub<
    'a,
    P: Program + Clone + Send + Sync + Serialize + for<'d> Deserialize<'d>,
    Q: Queue<P> + ?Sized,
> {
    jobs: Vec<Job<P>>,
    queue: &'a Q,
    dir: &'a str,
    counter: usize,
    proc: Procedure,
}

pub(crate) struct ResubOutput<P: Program + Clone + Send + Sync> {
    pub(crate) jobs: Vec<Job<P>>,
    pub(crate) slurm_jobs: HashMap<String, usize>,
    pub(crate) job_id: String,
    pub(crate) writing_input: Duration,
    pub(crate) writing_script: Duration,
    pub(crate) submitting: Duration,
}

impl<P: Program + Clone + Send + Sync> ResubOutput<P> {
    fn new(
        jobs: Vec<Job<P>>,
        slurm_jobs: HashMap<String, usize>,
        job_id: String,
        writing_input: Duration,
        writing_script: Duration,
        submitting: Duration,
    ) -> Self {
        Self {
            jobs,
            slurm_jobs,
            job_id,
            writing_input,
            writing_script,
            submitting,
        }
    }
}

impl<
    'a,
    P: Program + Clone + Send + Sync + Serialize + for<'d> Deserialize<'d>,
    Q: Queue<P> + ?Sized,
> Resub<'a, P, Q>
{
    pub(crate) fn new(queue: &'a Q, dir: &'a str, proc: Procedure) -> Self {
        Self {
            jobs: Vec::new(),
            queue,
            dir,
            counter: 0,
            proc,
        }
    }

    pub(crate) fn push(&mut self, job: Job<P>) {
        self.jobs.push(job)
    }

    pub(crate) fn resubmit(&mut self) -> Vec<ResubOutput<P>> {
        // this is inlined from Queue::resubmit minus actually submitting the
        // job. copy all of the original jobs to job_redo.ext
        for job in &mut self.jobs {
            let filename = format!(
                "{}.{}",
                job.program.filename(),
                job.program.extension()
            );
            let path = Path::new(&filename);
            let dir = path.parent().unwrap().to_str().unwrap();
            let base = path.file_stem().unwrap().to_str().unwrap();
            // nothing but the copy needs the name with extension
            let inp_name = format!("{dir}/{base}_redo");
            job.program.set_filename(&inp_name);
        }
        let mut jobs = std::mem::take(&mut self.jobs);
        jobs.chunks_mut(self.queue.chunk_size())
            .map(|jobs| {
                let (sj, wi, ws, ss) = self.queue.build_chunk_inner(
                    self.dir,
                    "redo",
                    self.counter,
                    jobs,
                    self.proc,
                );
                self.counter += 1;
                let job_id = jobs[0].job_id.clone();
                ResubOutput::new(jobs.to_vec(), sj, job_id, wi, ws, ss)
            })
            .collect()
    }
}
