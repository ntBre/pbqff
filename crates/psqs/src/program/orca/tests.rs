use std::{fs::read_to_string, str::FromStr};

use insta::{assert_debug_snapshot, assert_snapshot};
use tempfile::NamedTempFile;

use crate::{geom::Geom, program::Template};

use super::*;

fn test_job(template: &str) -> (NamedTempFile, Job) {
    let file = NamedTempFile::new().unwrap();
    let job = Job::new(
        file.path().to_string_lossy().to_string(),
        Template::from(template),
        -1,
        Geom::from_str(
            "3
water
O  0.00000000 -0.71603315  0.00000000
H  0.00000000 -0.14200298  0.77844804
H -0.00000000 -0.14200298 -0.77844804
",
        )
        .unwrap(),
        0,
    );
    (file, job)
}

#[test]
fn write_input() {
    let template = "! B3LYP def2-SVP VeryTightSCF {{.procedure}}

* xyz {{.charge}} 1
{{.geom}}*
";
    let (_file, job) = test_job(template);

    Orca.write_input(&job, Procedure::Opt);
    assert_snapshot!(read_to_string(Orca.infile(&job)).unwrap(), @r"
    ! B3LYP def2-SVP VeryTightSCF Opt

    * xyz -1 1
    O 0.000000000000 -0.716033150000 0.000000000000
    H 0.000000000000 -0.142002980000 0.778448040000
    H -0.000000000000 -0.142002980000 -0.778448040000
    *
    ");

    Orca.write_input(&job, Procedure::SinglePt);
    assert_snapshot!(read_to_string(Orca.infile(&job)).unwrap(), @r"
    ! B3LYP def2-SVP VeryTightSCF 

    * xyz -1 1
    O 0.000000000000 -0.716033150000 0.000000000000
    H 0.000000000000 -0.142002980000 0.778448040000
    H -0.000000000000 -0.142002980000 -0.778448040000
    *
    ");
}

#[test]
fn read_single_point_output() {
    let got = Orca.read_output("testfiles/orca/single").unwrap();
    assert_debug_snapshot!(got, @r"
    ProgramResult {
        energy: -385.369335800098,
        cart_geom: Some(
            [
                Atom {
                    atomic_number: 6,
                    x: -0.02274,
                    y: -0.007134,
                    z: -0.0,
                    weight: None,
                },
                Atom {
                    atomic_number: 1,
                    x: -0.0154,
                    y: -1.0799,
                    z: 0.0,
                    weight: None,
                },
            ],
        ),
        time: 48.811,
    }
    ");
}

#[test]
fn read_optimization_output_uses_final_energy_and_geometry() {
    let got = Orca.read_output("testfiles/orca/opt").unwrap();
    assert_debug_snapshot!(got, @r"
    ProgramResult {
        energy: -385.369299609594,
        cart_geom: Some(
            [
                Atom {
                    atomic_number: 6,
                    x: -0.022683,
                    y: -0.007064,
                    z: 0.0,
                    weight: None,
                },
                Atom {
                    atomic_number: 1,
                    x: -0.015291,
                    y: -1.079747,
                    z: 0.0,
                    weight: None,
                },
            ],
        ),
        time: 610.813,
    }
    ");
}

#[test]
fn rejects_failed_optimization_despite_normal_termination() {
    let got = Orca.read_output("testfiles/orca/failed_opt");
    assert!(matches!(got, Err(ProgramError::ErrorInOutput(_))));
}

#[test]
fn rejects_known_fatal_errors() {
    for output in [
        "INPUT ERROR\nUNRECOGNIZED KEYWORD",
        "Error (ORCA_MAIN): ... aborting the run",
        "ORCA finished by error termination",
    ] {
        assert!(matches!(
            parse_output(output, "error.out"),
            Err(ProgramError::ErrorInOutput(_))
        ));
    }
}

#[test]
fn does_not_accept_energy_before_termination() {
    let got = parse_output(
        "FINAL SINGLE POINT ENERGY      -385.369335800098\n",
        "unfinished.out",
    );
    assert!(matches!(got, Err(ProgramError::EnergyNotFound(_))));
}

#[test]
fn pbqff_energy_overrides_standard_energy() {
    let got = parse_output(
        "FINAL SINGLE POINT ENERGY (Solute-CPCM) -100.0
PBQFF = -101.25
****ORCA TERMINATED NORMALLY****
TOTAL RUN TIME: 0 days 0 hours 0 minutes 1 seconds 0 msec
",
        "compound.out",
    )
    .unwrap();
    assert_eq!(got.energy, -101.25);
}

#[test]
#[should_panic(expected = "ORCA templates must contain {{.procedure}}")]
fn requires_procedure_directive() {
    let (_file, job) = test_job("! B3LYP def2-SVP\n{{.geom}}");
    Orca.write_input(&job, Procedure::SinglePt);
}
