use super::*;

#[test]
fn config() {
    let got = Config::load("testfiles/test.toml");
    let want = Config {
        geometry: psqs::geom::Geom::Zmat(
            "C
C 1 CC
C 1 CC 2 CCC
H 2 CH 1 HCC 3 180.0
H 3 CH 1 HCC 2 180.0

CC =                  1.42101898
CCC =                55.60133141
CH =                  1.07692776
HCC =               147.81488230
"
            .to_string(),
        ),
        optimize: true,
        charge: 0,
        step_size: 0.005,
        coord_type: CoordType::Sic,
        template: String::from(
            "scfcrt=1.D-21 aux(precision=14 comp xp xs xw) PM6 THREADS=1",
        ),
        program: Program::Mopac,
        sleep_int: 2,
        job_limit: 2048,
        chunk_size: 1,
        queue: Queue::Slurm,
        findiff: false,
        check_int: 100,
        queue_template: None,
        hybrid_template: String::from(
            "scfcrt=1.D-21 aux(precision=14 comp xp xs xw) PM6 THREADS=1",
        ),
        weights: None,
        dummy_atoms: None,
        norm_resume_hff: false,
    };
    assert_eq!(got, want);
}

#[test]
fn path_template() {
    let got = Config::load("testfiles/path.toml");
    let want = Config {
        geometry: psqs::geom::Geom::Zmat(
            "C
C 1 CC
C 1 CC 2 CCC
H 2 CH 1 HCC 3 180.0
H 3 CH 1 HCC 2 180.0

CC =                  1.42101898
CCC =                55.60133141
CH =                  1.07692776
HCC =               147.81488230
"
            .to_string(),
        ),
        optimize: true,
        charge: 0,
        step_size: 0.005,
        coord_type: CoordType::Sic,
        template: String::from(
            "scfcrt=1.D-21 aux(precision=14 comp xp xs xw) PM6 THREADS=1",
        ),
        program: Program::Mopac,
        sleep_int: 2,
        job_limit: 2048,
        chunk_size: 1,
        queue: Queue::Slurm,
        findiff: false,
        check_int: 100,
        queue_template: None,
        hybrid_template: String::from(
            "scfcrt=1.D-21 aux(precision=14 comp xp xs xw) PM6 THREADS=1",
        ),
        weights: None,
        dummy_atoms: None,
        norm_resume_hff: false,
    };
    assert_eq!(got, want);
}

#[test]
fn normal() {
    let got = Config::load("testfiles/normal.toml");
    let want = Config {
        geometry: psqs::geom::Geom::Zmat(
            "O\nH 1 OH\nH 1 OH 2 HOH\n\nOH = 1.0\nHOH = 109.5\n".to_string(),
        ),
        optimize: true,
        charge: 0,
        step_size: 0.005,
        coord_type: CoordType::Normal,
        template: "scfcrt=1.D-21 aux(precision=14 comp xp xs xw) PM6 THREADS=1"
            .to_string(),
        hybrid_template:
            "scfcrt=1.D-21 aux(precision=14 comp xp xs xw) PM6 THREADS=1"
                .to_string(),
        queue_template: None,
        program: Program::Mopac,
        queue: Queue::Slurm,
        sleep_int: 2,
        job_limit: 2048,
        chunk_size: 1,
        findiff: false,
        check_int: 0,
        weights: None,
        dummy_atoms: None,
        norm_resume_hff: true,
    };
    assert_eq!(got, want);
}
