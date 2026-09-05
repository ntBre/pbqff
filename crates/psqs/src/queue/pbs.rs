use std::path::Path;
use std::time::Duration;
use std::{collections::HashSet, process::Command};

use crate::program::Program;
use crate::program::dftbplus::DFTBPlus;
use crate::program::molpro::Molpro;
use crate::program::mopac::Mopac;
use crate::queue::Queue;

use super::{SubQueue, Submit};

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

impl Submit<Mopac> for Pbs {
    /// submit `filename` to the queue and return the jobid
    fn submit(&self, filename: &str) -> String {
        let mut cmd =
            Command::new(<Self as SubQueue<Mopac>>::submit_command(self));
        let cmd = cmd.arg("-f").arg(filename);
        submit_inner(cmd, self.sleep_int).unwrap()
    }
}

// Molpro 2022 submit script requires submission from the current directory, so
// we have to override the default impl
impl Submit<Molpro> for Pbs {
    fn submit(&self, filename: &str) -> String {
        let path = Path::new(filename);
        let dir = path.parent().unwrap();
        let base = path.file_name().unwrap();
        let mut cmd =
            Command::new(<Self as SubQueue<Molpro>>::submit_command(self));
        let cmd = cmd.arg(base).current_dir(dir);
        submit_inner(cmd, self.sleep_int).unwrap()
    }
}

/// helper function to consolidate error handling between the two submit
/// implementations
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

impl Queue<Molpro> for Pbs {
    fn template(&self) -> &str {
        &self.template
    }

    /// This must be consistent with the `Submit<Molpro>` implementation, which
    /// currently changes to the parent directory of the PBS script before
    /// submitting. This also assumes, then, that the PBS script is in the same
    /// directory as the input files, but I think that's a safe assumption.
    fn program_cmd(&self, filename: &str) -> String {
        let basename = Path::new(&filename).file_name().unwrap();
        format!("$MOLPRO_CMD {basename:?}.inp")
    }
}

impl Queue<Mopac> for Pbs {
    fn template(&self) -> &str {
        &self.template
    }

    fn program_cmd(&self, filename: &str) -> String {
        format!("$MOPAC_CMD {filename}.mop")
    }
}

impl Queue<DFTBPlus> for Pbs {
    fn template(&self) -> &str {
        &self.template
    }

    fn program_cmd(&self, filename: &str) -> String {
        format!("(cd {filename} && $DFTB_CMD > out)")
    }
}

impl Submit<DFTBPlus> for Pbs {
    fn submit(&self, filename: &str) -> String {
        let mut cmd =
            Command::new(<Self as SubQueue<DFTBPlus>>::submit_command(self));
        let cmd = cmd.arg("-f").arg(filename);
        submit_inner(cmd, self.sleep_int).unwrap()
    }
}

impl<P> SubQueue<P> for Pbs
where
    P: Program,
{
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

    const SCRIPT_EXT: &'static str = "pbs";

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
        let lines = <Pbs as SubQueue<P>>::stat_cmd(self);
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

    use crate::{program::cfour::Cfour, queue::templates};

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
        ($($name:ident, $queue:expr => $p:ty$(,)*)*) => {
            $(
            #[test]
            fn $name() {
                let tmp = tempfile::NamedTempFile::new().unwrap();
                <Pbs as Queue<$p>>::write_submit_script(
                    $queue,
                    ["pts/opt0.inp", "pts/opt1.inp", "pts/opt2.inp", "pts/opt3.inp"]
                    .map(|s| s.into()),
                    tmp.path().to_str().unwrap(),
                );
                let got = std::fs::read_to_string(tmp).unwrap();
                let got: Vec<&str> = got.lines().filter(|l|
                    !(l.starts_with("#PBS -N")
                        || l.starts_with("#PBS -o"))).collect();
                let got = got.join("\n");
                assert_snapshot!(got);
            }
            )*
        }
    }

    make_tests! {
        mopac_pbs, &pbs(templates::PBS_MOPAC) =>  Mopac,
        molpro_pbs, &pbs(templates::PBS_MOLPRO) =>  Molpro,
        cfour_pbs, &pbs(templates::PBS_CFOUR) => Cfour,
        dftb_pbs, &pbs(templates::PBS_DFTBPLUS) => DFTBPlus,
    }
}
