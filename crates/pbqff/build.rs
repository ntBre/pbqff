use std::{env, error::Error, ffi::OsString, fs, path::Path, process::Command};

fn main() {
    println!("cargo:rerun-if-changed=.git/index");
    let out_dir = env::var_os("OUT_DIR").unwrap();
    version(&out_dir);
}

fn make_id() -> Result<String, Box<dyn Error>> {
    let cmd = Command::new("git").arg("rev-parse").arg("HEAD").output()?;
    Ok(String::from_utf8(cmd.stdout[..8].to_vec())?)
}

fn version(out_dir: &OsString) {
    let dest_path = Path::new(&out_dir).join("version.rs");
    let id = make_id().unwrap_or_else(|_| "deadbeef".to_string());
    fs::write(
        dest_path,
        format!("pub fn version() -> &'static str {{ \"{id}\" }}"),
    )
    .unwrap();
}
