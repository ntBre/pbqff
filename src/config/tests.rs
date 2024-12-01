use insta::assert_debug_snapshot;
use test_case::test_case;

use super::*;

#[test_case("testfiles/test.toml")]
#[test_case("testfiles/path.toml")]
#[test_case("testfiles/normal.toml")]
fn load_config(path: &str) {
    assert_debug_snapshot!(Config::load(path));
}
