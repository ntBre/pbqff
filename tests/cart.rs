use std::{fs::read_to_string, path::Path};

use assert_cmd::Command;
use insta::assert_snapshot;
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

    assert_snapshot!(read_to_string(dir.path().join("pbqff.out"))?);

    Ok(())
}
