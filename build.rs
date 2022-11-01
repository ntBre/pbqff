use std::{env, fs, path::Path};

fn main() {
    println!("cargo:rerun-if-changed=.git/index");
    let repo = git2::Repository::discover(".").unwrap();
    let head = repo.head().unwrap();
    let id = &head.peel_to_commit().unwrap().id().to_string()[..8];
    // revparse("HEAD");
    let out_dir = env::var_os("OUT_DIR").unwrap();
    let dest_path = Path::new(&out_dir).join("version.rs");
    fs::write(
        &dest_path,
        format!(
            "pub fn version() -> &'static str {{
	    \"{}\"
	}}
	",
            id
        ),
    )
    .unwrap();
}
