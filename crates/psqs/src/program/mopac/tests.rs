use std::collections::HashSet;
use std::fs;
use std::ops::{Deref, DerefMut};

use insta::{assert_debug_snapshot, assert_snapshot};
use tempfile::TempDir;

use crate::string;

use crate::queue::{Queue, SubQueue, Submit};

use super::*;

struct TestMopac {
    tempdir: TempDir,
    mopac: Mopac,
}

impl TestMopac {
    #[must_use]
    fn with_params(mut self, params: Option<Params>) -> Self {
        self.mopac.params = params;
        self
    }

    /// Set [`Mopac::param_dir`] to `Some(self.tempdir)`.
    #[must_use]
    fn with_param_dir(mut self) -> Self {
        self.param_dir =
            Some(self.tempdir.path().to_string_lossy().to_string());
        self
    }
}

impl Deref for TestMopac {
    type Target = Mopac;

    fn deref(&self) -> &Self::Target {
        &self.mopac
    }
}

impl DerefMut for TestMopac {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.mopac
    }
}

fn test_mopac() -> TestMopac {
    let names = vec![
        "USS", "ZS", "BETAS", "GSS", "USS", "UPP", "ZS", "ZP", "BETAS",
        "BETAP", "GSS", "GPP", "GSP", "GP2", "HSP",
    ];
    let atoms = vec![
        "H", "H", "H", "H", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C",
        "C",
    ];
    #[rustfmt::skip]
    let values = vec![
        -11.246958000000, 1.268641000000, -8.352984000000,
        14.448686000000, -51.089653000000, -39.937920000000,
        2.047558000000, 1.702841000000, -15.385236000000,
        -7.471929000000, 13.335519000000, 10.778326000000,
        11.528134000000, 9.486212000000, 0.717322000000,
    ];
    let tempdir = TempDir::new().unwrap();
    let mopac = Mopac::new_full(
        tempdir.path().join("test").to_string_lossy().to_string(),
        Some(Params::from(
            names.iter().map(|s| s.to_string()).collect(),
            atoms.iter().map(|s| s.to_string()).collect(),
            values,
        )),
        Geom::Xyz(Vec::new()),
        0,
        Template::from("scfcrt=1.D-21 aux(precision=14) PM6 A0"),
    );

    TestMopac { tempdir, mopac }
}

#[test]
fn test_write_input() {
    let mut tm = test_mopac().with_params(None);
    tm.write_input(Procedure::SinglePt);
    let got = fs::read_to_string(tm.infile()).expect("file not found");
    let want = "scfcrt=1.D-21 aux(precision=14) PM6 A0 charge=0 1SCF XYZ
Comment line 1
Comment line 2

"
    .to_string();
    assert_eq!(got, want);
}

#[test]
fn test_write_input_with_params() {
    let mut tm = test_mopac().with_param_dir();
    tm.write_input(Procedure::SinglePt);

    insta::with_settings!({filters => vec![(
        tm.param_file.as_ref().unwrap().as_str(), "[PARAM_FILE]",
    )]}, {
        assert_snapshot!(read_to_string(tm.infile()).expect("file not found"), @r"
        scfcrt=1.D-21 aux(precision=14) PM6 A0 charge=0 1SCF external=[PARAM_FILE] XYZ
        Comment line 1
        Comment line 2
        ");
    });
}

#[test]
fn test_write_params() {
    let tm = test_mopac();
    let param_file = tm.tempdir.path().join("params.dat");
    Mopac::write_params(tm.params.as_ref().unwrap(), &param_file);

    assert_snapshot!(read_to_string(&param_file).expect("file not found"), @r"
    USS H -11.246958000000
    ZS H 1.268641000000
    BETAS H -8.352984000000
    GSS H 14.448686000000
    USS C -51.089653000000
    UPP C -39.937920000000
    ZS C 2.047558000000
    ZP C 1.702841000000
    BETAS C -15.385236000000
    BETAP C -7.471929000000
    GSS C 13.335519000000
    GPP C 10.778326000000
    GSP C 11.528134000000
    GP2 C 9.486212000000
    HSP C 0.717322000000
    ");
}

#[test]
fn test_read_output() {
    assert_debug_snapshot!(Mopac::read_output("testfiles/job"), @r"
    Ok(
        ProgramResult {
            energy: 0.15478330901845888,
            cart_geom: Some(
                [
                    Atom {
                        atomic_number: 6,
                        x: 0.0,
                        y: 0.001986048947850556,
                        z: -0.8876484030621448,
                        weight: None,
                    },
                    Atom {
                        atomic_number: 6,
                        x: 0.0,
                        y: 0.6656848670392688,
                        z: 0.36485436700254253,
                        weight: None,
                    },
                    Atom {
                        atomic_number: 6,
                        x: 0.0,
                        y: -0.6648176873614392,
                        z: 0.3712115291247976,
                        weight: None,
                    },
                    Atom {
                        atomic_number: 1,
                        x: 0.0,
                        y: 1.6004641296483408,
                        z: 0.9066925741481017,
                        weight: None,
                    },
                    Atom {
                        atomic_number: 1,
                        x: 0.0,
                        y: -1.6033173582740208,
                        z: 0.9065798886662334,
                        weight: None,
                    },
                ],
            ),
            time: 0.015625,
        },
    )
    ");

    // opt success
    assert_debug_snapshot!(Mopac::read_output("testfiles/opt"), @r"
    Ok(
        ProgramResult {
            energy: 0.20175470737510037,
            cart_geom: Some(
                [
                    Atom {
                        atomic_number: 6,
                        x: 0.0,
                        y: 0.0,
                        z: 0.0,
                        weight: None,
                    },
                    Atom {
                        atomic_number: 6,
                        x: 1.4361996438838212,
                        y: 0.0,
                        z: 0.0,
                        weight: None,
                    },
                    Atom {
                        atomic_number: 6,
                        x: 0.7993316223304503,
                        y: 1.1932050849014117,
                        z: 0.0,
                        weight: None,
                    },
                    Atom {
                        atomic_number: 1,
                        x: 2.360710453618393,
                        y: -0.5060383602977097,
                        z: 2.6804e-14,
                        weight: None,
                    },
                    Atom {
                        atomic_number: 1,
                        x: 0.8934572415091369,
                        y: 2.2429362062954086,
                        z: -2.6804e-14,
                        weight: None,
                    },
                ],
            ),
            time: 0.03125,
        },
    )
    ");

    // failure (no termination message) in output - now catches noaux error
    // instead
    let f = String::from("testfiles/nojob");
    let got = Mopac::read_output(&f);
    assert_eq!(got.err().unwrap(), ProgramError::FileNotFound(f + ".aux"));

    // failure in aux
    let f = String::from("testfiles/noaux");
    let got = Mopac::read_output(&f);
    assert_eq!(got.err().unwrap(), ProgramError::FileNotFound(f + ".aux"));

    // this is passing but for some reason on maple it's throwing an error
    let got = Mopac::read_output("testfiles/bad");
    assert!(got.is_ok());
    assert!(got.unwrap().cart_geom.is_some());
}

#[test]
fn read_multi_el() {
    let got = Mopac::read_output("testfiles/mopac/multi_atom_el")
        .map(|pr| pr.cart_geom.unwrap().len());
    assert_eq!(got, Ok(46));
}

/// minimal queue for testing general submission
struct TestQueue;

impl Submit<Mopac> for TestQueue {}

impl Queue<Mopac> for TestQueue {
    fn template(&self) -> &Option<String> {
        static S: Option<String> = Some(String::new());
        &S
    }

    fn program_cmd(&self, filename: &str) -> String {
        format!("echo {filename}")
    }

    fn default_submit_script(&self) -> String {
        todo!()
    }
}

impl SubQueue<Mopac> for TestQueue {
    fn submit_command(&self) -> &str {
        "bash"
    }

    fn chunk_size(&self) -> usize {
        128
    }

    fn job_limit(&self) -> usize {
        1600
    }

    fn sleep_int(&self) -> usize {
        1
    }

    const SCRIPT_EXT: &'static str = "pbs";

    fn dir(&self) -> &str {
        "inp"
    }

    fn stat_cmd(&self) -> String {
        todo!()
    }

    fn status(&self) -> HashSet<String> {
        todo!()
    }

    fn no_del(&self) -> bool {
        false
    }
}

#[test]
fn test_submit() {
    let tmp = TempDir::new().unwrap();
    let pbs = tmp.path().join("main.pbs");
    let path = pbs.to_string_lossy();
    let tq = TestQueue;
    tq.write_submit_script(
        string!["input1.mop", "input2.mop", "input3.mop"],
        &path,
    );
    let got = tq.submit(&path);
    let want = "input3.mop";
    assert_eq!(got, want);
}

#[test]
fn test_resubmit() -> anyhow::Result<()> {
    let tmp = TempDir::new()?;
    let mop = tmp.path().join("job.mop");
    let redo_mop = tmp.path().join("job_redo.mop");
    let redo_pbs = tmp.path().join("job_redo.pbs");

    std::fs::copy("testfiles/job.mop", &mop)?;

    let got = TestQueue.resubmit(&mop);

    assert!(redo_mop.exists());
    assert!(redo_pbs.exists());
    assert_eq!(read_to_string(&mop)?, read_to_string(&redo_mop)?);

    insta::with_settings!({filters => vec![
        (tmp.path().to_str().unwrap(), "[TMP]"),
    ]}, {
        assert_debug_snapshot!(got, @r#"
        Resubmit {
            inp_file: "[TMP]/job_redo",
            pbs_file: "[TMP]/job_redo.pbs",
            job_id: "[TMP]/job_redo",
        }
        "#);
    });

    Ok(())
}
