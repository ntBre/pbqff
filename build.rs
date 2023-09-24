use std::{env, ffi::OsString, fs, path::Path};

#[cfg(not(feature = "vers"))]
fn make_id() -> &'static str {
    "deadbeef"
}

#[cfg(feature = "vers")]
fn make_id() -> String {
    let repo = git2::Repository::discover(".").unwrap();
    let head = repo.head().unwrap();
    let id = &head.peel_to_commit().unwrap().id().to_string()[..8];
    id.to_owned()
}

fn version(out_dir: OsString) {
    let dest_path = Path::new(&out_dir).join("version.rs");
    let id = make_id();
    fs::write(
        dest_path,
        format!(
            "pub fn version() -> &'static str {{
	    \"{id}\"
	}}
	"
        ),
    )
    .unwrap();
}

fn main() {
    println!("cargo:rerun-if-changed=.git/index");
    let out_dir = env::var_os("OUT_DIR").unwrap();
    version(out_dir);
}
