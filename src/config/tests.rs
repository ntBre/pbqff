use insta::assert_debug_snapshot;

use super::*;

#[test]
fn config() {
    assert_debug_snapshot!(Config::load("testfiles/test.toml"));
}

#[test]
fn path_template() {
    assert_debug_snapshot!(Config::load("testfiles/path.toml"));
}

#[test]
fn normal() {
    assert_debug_snapshot!(Config::load("testfiles/normal.toml"));
}
