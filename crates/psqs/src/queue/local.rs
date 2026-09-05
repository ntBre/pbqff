use std::collections::HashSet;

use crate::queue::Queue;

/// Minimal implementation for running programs locally, primarily for tests.
#[derive(Debug)]
pub struct Local {
    pub dir: String,
    pub chunk_size: usize,
    pub template: String,
}

impl Default for Local {
    fn default() -> Self {
        Self {
            dir: ".".to_string(),
            chunk_size: 128,
            template: String::new(),
        }
    }
}

impl Local {
    pub fn new(
        chunk_size: usize,
        _job_limit: usize,
        _sleep_int: usize,
        dir: &'static str,
        _no_del: bool,
        template: String,
    ) -> Self {
        Self {
            dir: dir.to_string(),
            chunk_size,
            template,
        }
    }
}

impl Queue for Local {
    fn script_ext(&self) -> &'static str {
        "slurm"
    }

    fn template(&self) -> &str {
        &self.template
    }

    fn submit_command(&self) -> &str {
        "bash"
    }

    fn chunk_size(&self) -> usize {
        self.chunk_size
    }

    fn job_limit(&self) -> usize {
        1600
    }

    fn sleep_int(&self) -> usize {
        1
    }

    fn dir(&self) -> &str {
        &self.dir
    }

    fn stat_cmd(&self) -> String {
        todo!()
    }

    fn status(&self) -> HashSet<String> {
        for dir in ["opt", "pts", "freqs"] {
            let Ok(d) = std::fs::read_dir(dir) else {
                log::error!("{dir} not found for status");
                continue;
            };
            for f in d {
                eprintln!("contents of {:?}", f.as_ref().unwrap());
                eprintln!(
                    "{}",
                    std::fs::read_to_string(f.unwrap().path()).unwrap()
                );
                eprintln!("================");
            }
        }
        panic!("no status available for Local queue");
    }

    fn no_del(&self) -> bool {
        false
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

    fn local(template: &str) -> Local {
        Local {
            dir: String::new(),
            chunk_size: 0,
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
                let got: Vec<&str> = got.lines().filter(|l|
                    !l.contains("/tmp")).collect();
                let got = got.join("\n");
                assert_snapshot!(got);
            }
            )*
        }
    }

    make_tests! {
        mopac_local, &local(templates::LOCAL_MOPAC) => Mopac,
        molpro_local, &local(templates::LOCAL_MOLPRO) => Molpro,
        cfour_local, &local(templates::LOCAL_CFOUR) => Cfour,
        dftb_local, &local(templates::LOCAL_DFTBPLUS) => DFTBPlus,
        orca_local, &local(templates::LOCAL_ORCA) => Orca,
    }
}
