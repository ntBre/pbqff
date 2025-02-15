use std::{fs::File, os::unix::prelude::AsRawFd, path::Path};

use pbqff::{
    cleanup,
    config::{self, Config},
    coord_type::{normal::Normal, Cart, CoordType, Load, Sic},
    die, Intder,
};
use psqs::{
    program::{cfour::Cfour, dftbplus::DFTBPlus, molpro::Molpro, mopac::Mopac},
    queue::{local::Local, pbs::Pbs, slurm::Slurm},
};

include!(concat!(env!("OUT_DIR"), "/version.rs"));
include!(concat!(env!("OUT_DIR"), "/dispatch.rs"));

use clap::Parser;

/// push button quartic force fields
#[derive(Parser, Debug)]
#[command(author, about, long_about = None)]
struct Args {
    /// input file
    #[arg(value_parser, default_value_t = String::from("pbqff.toml"))]
    infile: String,

    /// Resume from the checkpoint files in the current directory (chk.json and
    /// res.chk). Defaults to false.
    #[arg(short, long, default_value_t = false)]
    checkpoint: bool,

    /// Overwrite existing output from a previous run. Defaults to false.
    #[arg(short, long, default_value_t = false)]
    overwrite: bool,

    /// Print the git commit hash and exit. Defaults to false.
    #[arg(short, long, default_value_t = false)]
    version: bool,

    /// Set the maximum number of threads to use. Defaults to 0, which means to
    /// use as many threads as there are CPUS.
    #[arg(short, long, default_value_t = 0)]
    threads: usize,

    /// Serialize the input file to JSON and exit. For use by qffbuddy.
    #[arg(short, default_value_t = false, hide = true)]
    json: bool,

    /// Don't delete any files when running the single-point energies. Defaults
    /// to false.
    #[arg(short, long, default_value_t = false)]
    no_del: bool,
}

use spectro::{Output, Spectro};

fn main() -> Result<(), std::io::Error> {
    env_logger::init();
    let args = Args::parse();
    if args.version {
        println!("version: {}", version());
        return Ok(());
    }
    if args.json {
        let config = Config::load(&args.infile);
        match serde_json::to_string(&config) {
            Ok(s) => println!("{}", s),
            Err(e) => {
                die!("failed to deserialize {} with {e}", args.infile);
            }
        };
        return Ok(());
    }
    let path = Path::new("pbqff.out");
    if path.exists() && !args.overwrite {
        die!("existing pbqff output. overwrite with -o/--overwrite");
    }
    let outfile = File::create(path).expect("failed to create outfile");
    let logfile = File::create("pbqff.log").expect("failed to create log file");
    let out_fd = outfile.as_raw_fd();
    let log_fd = logfile.as_raw_fd();
    // redirect stdout to outfile and stderr to logfile
    unsafe {
        libc::dup2(out_fd, 1);
        libc::dup2(log_fd, 2);
    }
    let config = Config::load(&args.infile);
    println!("PID: {}", std::process::id());
    println!("version: {}", version());
    psqs::max_threads(args.threads);
    cleanup(".");
    let _ = std::fs::create_dir("pts");

    let (spectro, output) = dispatch(&config, args);

    spectro.write_output(&mut std::io::stdout(), &output)?;

    let mut f = std::fs::File::create("spectro.json")?;
    use std::io::Write;
    writeln!(f, "{}", serde_json::to_string_pretty(&output)?)?;

    println!("normal termination of pbqff");

    Ok(())
}
