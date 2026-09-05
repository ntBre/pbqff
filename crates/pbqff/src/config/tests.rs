use insta::{assert_debug_snapshot, with_settings};
use test_case::test_case;

use super::*;

#[test_case("testfiles/test.toml" ; "basic sic")]
#[test_case("testfiles/cart.toml" ; "basic cart")]
#[test_case("testfiles/normal.toml" ; "basic norm")]
#[test_case("testfiles/path.toml" ; "path templates")]
fn load_config(path: &str) {
    with_settings!({ snapshot_suffix => path }, {
        assert_debug_snapshot!(Config::load(path));
    });
}

#[test]
fn resolves_bundled_queue_templates() {
    use psqs::queue::templates;

    let cases = [
        (Queue::Pbs, Program::Mopac, templates::PBS_MOPAC),
        (Queue::Pbs, Program::Molpro, templates::PBS_MOLPRO),
        (Queue::Pbs, Program::DFTBPlus, templates::PBS_DFTBPLUS),
        (Queue::Pbs, Program::Cfour, templates::PBS_CFOUR),
        (Queue::Pbs, Program::Orca, templates::PBS_ORCA),
        (Queue::Slurm, Program::Mopac, templates::SLURM_MOPAC),
        (Queue::Slurm, Program::Molpro, templates::SLURM_MOLPRO),
        (Queue::Slurm, Program::DFTBPlus, templates::SLURM_DFTBPLUS),
        (Queue::Slurm, Program::Cfour, templates::SLURM_CFOUR),
        (Queue::Slurm, Program::Orca, templates::SLURM_ORCA),
        (Queue::Local, Program::Mopac, templates::LOCAL_MOPAC),
        (Queue::Local, Program::Molpro, templates::LOCAL_MOLPRO),
        (Queue::Local, Program::DFTBPlus, templates::LOCAL_DFTBPLUS),
        (Queue::Local, Program::Cfour, templates::LOCAL_CFOUR),
        (Queue::Local, Program::Orca, templates::LOCAL_ORCA),
    ];

    for (queue, program, expected) in cases {
        let config = Config {
            queue,
            program,
            ..Config::default()
        };
        assert_eq!(config.resolved_queue_template(), expected);
    }
}

#[test]
fn explicit_queue_template_takes_precedence() {
    let config = Config {
        queue_template: Some("explicit template".to_owned()),
        ..Config::default()
    };

    assert_eq!(config.resolved_queue_template(), "explicit template");
}
