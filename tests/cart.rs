use std::{fs::read_to_string, path::Path};

use assert_cmd::Command;
use insta::{assert_snapshot, with_settings};
use tempfile::tempdir;

#[test]
fn cart() -> std::io::Result<()> {
    let config_file = Path::new("testfiles/cart.toml");
    let dir = tempdir()?;
    std::fs::copy(config_file, dir.path().join("cart.toml"))?;
    let mut cmd = Command::cargo_bin("pbqff").unwrap();
    let assert = cmd.arg("cart.toml").current_dir(&dir).assert();
    let output = assert.get_output();

    assert!(
        output.status.success(),
        "stderr: {}\nlog: {}",
        String::from_utf8_lossy(&output.stderr),
        read_to_string(dir.path().join("pbqff.log"))?,
    );

    // filter out essentially all of the numerical results just to check that
    // all of the sections are there
    with_settings!(
        {filters => vec![
            (r"(?m)^PID: \d+$", "PID: [PID]"),
            (r"(?s)(normalized geometry:).*(point group)", "$1[Normalized Geometry]\n$2"),
            (r"(?ms)^STATE NO.*?^$", "[Vibrational States]"),
            (r"(?s)(Transformed Geometry \(Å\):).*(Equilibrium)", "$1\n[Transformed Geometry]\n$2"),
            (r"(?ms)(Equilibrium Rotational Constants \(cm-1\):).*?^$", "$1\n[Rotational Constants]\n"),
            (r"(?s)(Geometry:).*(Vibrational)", "$1\n[Geometry]\n$2"),
            (r"(?s)(Vibrational Frequencies \(cm-1\):).*?(ZPT)", "$1\n[Freqs]\n$2"),
            (r"(?ms)(^Rotational Constants \(cm-1\):).*?(Coriolis Resonances)", "$1\n[Rotational Constants]\n$2"),
        ]}, {
            assert_snapshot!(read_to_string(dir.path().join("pbqff.out")).unwrap());
        }
    );

    Ok(())
}
