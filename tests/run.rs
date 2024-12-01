use std::{fs::read_to_string, path::Path};

use assert_cmd::Command;
use insta::{assert_snapshot, with_settings};
use tempfile::tempdir;
use test_case::test_case;

#[test_case("testfiles/test.toml", &["testfiles/intder.in"]; "sic")]
#[test_case("testfiles/cart.toml", &[]; "cart")]
#[test_case("testfiles/norm.toml", &[]; "norm")]
fn run(path: &str, other_files: &[&str]) -> std::io::Result<()> {
    let config_file = Path::new(path);
    let dir = tempdir()?;
    std::fs::copy(config_file, dir.path().join("pbqff.toml"))?;
    for file in other_files {
        let p = Path::new(file);
        let filename = p.file_name().unwrap();
        std::fs::copy(file, dir.path().join(filename))?;
    }
    let mut cmd = Command::cargo_bin("pbqff").unwrap();
    let assert = cmd.arg("pbqff.toml").current_dir(&dir).assert();
    let output = assert.get_output();

    assert!(
        output.status.success(),
        "stderr: {}\nlog: {}",
        String::from_utf8_lossy(&output.stderr),
        read_to_string(dir.path().join("pbqff.log"))?,
    );

    // filter out essentially all of the numerical results just to check that
    // all of the sections are there
    with_settings!({
        filters => vec![
            (r"(?m)^PID: \d+$", "PID: [PID]"),
            (r"(?m)^version: [a-z0-9]+$", "version: [version]"),
            (r"(?s)(normalized geometry:).*(point group)", "$1[Normalized Geometry]\n$2"),
            (r"(?ms)^Normal Coordinates:.*?^$", "[Normal Coordinates]"),
            (r"(?ms)^STATE NO.*?^$", "[Vibrational States]"),
            (r"(?s)(Transformed Geometry \(Å\):).*(Equilibrium)", "$1\n[Transformed Geometry]\n$2"),
            (r"(?ms)(Equilibrium Rotational Constants \(cm-1\):).*?^$", "$1\n[Rotational Constants]\n"),
            (r"(?s)(Geometry:).*(Vibrational)", "$1\n[Geometry]\n$2"),
            (r"(?s)(Vibrational Frequencies \(cm-1\):).*?(ZPT)", "$1\n[Freqs]\n$2"),
            (r"(?ms)(^Rotational Constants \(cm-1\):).*?(Coriolis Resonances)", "$1\n[Rotational Constants]\n$2"),
            (r"(?m)^(anpass sum of squared residuals):(.+)$", "$1: [residual]"),
        ],
        snapshot_suffix => path

    }, {
        assert_snapshot!(read_to_string(dir.path().join("pbqff.out")).unwrap());
    });

    Ok(())
}
