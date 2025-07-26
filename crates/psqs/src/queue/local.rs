use serde::{Deserialize, Serialize};
use std::collections::HashSet;

use crate::program::dftbplus::DFTBPlus;
use crate::program::molpro::Molpro;
use crate::program::{Program, mopac::Mopac};
use crate::queue::Queue;

use super::{SubQueue, Submit};

/// Minimal implementation for testing MOPAC locally
#[derive(Debug)]
pub struct Local {
    pub dir: String,
    pub chunk_size: usize,
    pub template: Option<String>,
}

impl Default for Local {
    fn default() -> Self {
        Self {
            dir: ".".to_string(),
            chunk_size: 128,
            template: None,
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
        template: Option<String>,
    ) -> Self {
        Self {
            dir: dir.to_string(),
            chunk_size,
            template,
        }
    }
}

impl Submit<Molpro> for Local {}

impl Queue<Molpro> for Local {
    fn template(&self) -> &Option<String> {
        &self.template
    }

    fn program_cmd(&self, filename: &str) -> String {
        format!("$MOLPRO_CMD {filename}.inp")
    }

    fn default_submit_script(&self) -> String {
        String::new()
    }
}

impl Submit<Mopac> for Local {}

impl Queue<Mopac> for Local {
    fn template(&self) -> &Option<String> {
        &self.template
    }

    fn program_cmd(&self, filename: &str) -> String {
        format!("$MOPAC_CMD {filename}.mop")
    }

    fn default_submit_script(&self) -> String {
        "export MOPAC_CMD=/opt/mopac/mopac
export LD_LIBRARY_PATH=/opt/mopac/\n"
            .into()
    }
}

impl Submit<DFTBPlus> for Local {}

impl Queue<DFTBPlus> for Local {
    fn template(&self) -> &Option<String> {
        &self.template
    }

    fn program_cmd(&self, filename: &str) -> String {
        format!("(cd {filename} && $DFTB_CMD > out)")
    }

    fn default_submit_script(&self) -> String {
        "DFTB_CMD=/opt/dftb+/dftb+\n".into()
    }
}

impl<P: Program + Clone + Serialize + for<'a> Deserialize<'a>> SubQueue<P>
    for Local
{
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

    const SCRIPT_EXT: &'static str = "slurm";

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

    use crate::program::cfour::Cfour;

    use super::*;

    fn local() -> Local {
        Local {
            dir: String::new(),
            chunk_size: 0,
            template: None,
        }
    }

    macro_rules! make_tests {
        ($($name:ident, $queue:expr => $p:ty$(,)*)*) => {
            $(
            #[test]
            fn $name() {
                let tmp = tempfile::NamedTempFile::new().unwrap();
                <Local as Queue<$p>>::write_submit_script(
                    $queue,
                    ["opt0.inp", "opt1.inp", "opt2.inp", "opt3.inp"].map(|s| s.into()),
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
        mopac_local, &local() =>  Mopac,
        molpro_local, &local() =>  Molpro,
        cfour_local, &local() => Cfour,
        dftb_local, &local() => DFTBPlus,
    }
}
