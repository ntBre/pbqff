use std::{env, error::Error, ffi::OsString, fs, path::Path, process::Command};

fn main() {
    println!("cargo:rerun-if-changed=.git/index");
    let out_dir = env::var_os("OUT_DIR").unwrap();
    version(&out_dir);
    dispatch(&out_dir);
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

/// construct a big match expression handling all known combinations of
/// coordinate types (Normal, Cart, and Sic), programs (Mopac and Molpro), and
/// queues (Pbs, Slurm, and Local)
pub fn dispatch(out_dir: &OsString) {
    let coord_types = ["Normal", "Cart", "Sic"];
    let coord_builders = [
        "Normal::findiff(config.findiff)",
        "Cart",
        "Sic::new(Intder::load_file(\"intder.in\"))",
    ];
    let programs = ["Mopac", "Molpro", "DFTBPlus", "Cfour"];
    let queues = ["Pbs", "Slurm", "Local"];
    use std::fmt::Write;
    let mut s = String::from(
        "fn dispatch(config: &Config, args: Args) -> (Spectro, Output) {
    let m = (config.coord_type, config.program, config.queue, args.checkpoint);
match m {
",
    );
    for (i, coord) in coord_types.iter().enumerate() {
        for program in programs {
            for queue in queues {
                for (val, fun) in [(true, "resume"), (false, "run")] {
                    let mut resume = String::from("config,");
                    if val {
                        write!(
                            resume,
                            "<{coord} as CoordType<::std::io::Stdout, \
			     {queue}, {program}>>::Resume::load(\"res.chk\")"
                        )
                        .unwrap();
                    }
                    write!(
                        s,
                        "(
config::CoordType::{coord},
config::Program::{program},
config::Queue::{queue},
{val},
) => <{coord} as CoordType<_, _, {program}>>::{fun}(
{},
\".\",
&mut std::io::stdout(),
&{queue}::new(
    config.chunk_size,
    config.job_limit,
    config.sleep_int,
    \"pts\",
    args.no_del,
    config.queue_template.clone(),
),
{resume}
),",
                        coord_builders[i],
                    )
                    .unwrap();
                }
            }
        }
    }
    write!(s, "}}}}").unwrap();
    let dest_path = Path::new(&out_dir).join("dispatch.rs");
    fs::write(dest_path, s).unwrap();
}
