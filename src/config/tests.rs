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
    };
    assert_eq!(got, want);
}
