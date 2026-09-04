use std::{collections::HashMap, path::Path, time::Duration};

use crate::{
    program::{Job, Procedure, Program},
    queue::Queue,
};

pub(crate) struct Resub<'a, Q: Queue + ?Sized> {
    jobs: Vec<Job>,
    program: &'a dyn Program,
    queue: &'a Q,
    dir: &'a str,
    counter: usize,
    proc: Procedure,
}

pub(crate) struct ResubOutput {
    pub(crate) jobs: Vec<Job>,
    pub(crate) slurm_jobs: HashMap<String, usize>,
    pub(crate) job_id: String,
    pub(crate) writing_input: Duration,
    pub(crate) writing_script: Duration,
    pub(crate) submitting: Duration,
}

impl ResubOutput {
    fn new(
        jobs: Vec<Job>,
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

impl<'a, Q: Queue + ?Sized> Resub<'a, Q> {
    pub(crate) fn new(
        program: &'a dyn Program,
        queue: &'a Q,
        dir: &'a str,
        proc: Procedure,
    ) -> Self {
        Self {
            jobs: Vec::new(),
            program,
            queue,
            dir,
            counter: 0,
            proc,
        }
    }

    pub(crate) fn push(&mut self, job: Job) {
        self.jobs.push(job)
    }

    pub(crate) fn resubmit(&mut self) -> Vec<ResubOutput> {
        // this is inlined from Queue::resubmit minus actually submitting the
        // job. copy all of the original jobs to job_redo.ext
        for job in &mut self.jobs {
            let filename =
                format!("{}.{}", job.filename, self.program.extension());
            let path = Path::new(&filename);
            let dir = path.parent().unwrap().to_str().unwrap();
            let base = path.file_stem().unwrap().to_str().unwrap();
            // nothing but the copy needs the name with extension
            let inp_name = format!("{dir}/{base}_redo");
            job.filename = inp_name;
        }
        let mut jobs = std::mem::take(&mut self.jobs);
        jobs.chunks_mut(self.queue.chunk_size())
            .map(|jobs| {
                let (sj, wi, ws, ss) = self.queue.build_chunk_inner(
                    self.program,
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
