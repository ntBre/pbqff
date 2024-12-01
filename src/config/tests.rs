use insta::{assert_debug_snapshot, with_settings};
use test_case::test_case;

use super::*;

#[test_case("testfiles/test.toml" ; "cart")]
#[test_case("testfiles/path.toml" ; "path templates")]
#[test_case("testfiles/normal.toml" ; "norm")]
fn load_config(path: &str) {
    with_settings!({ snapshot_suffix => path }, {
        assert_debug_snapshot!(Config::load(path));
    });
}
