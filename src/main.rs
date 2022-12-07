use std::{fs::File, os::unix::prelude::AsRawFd, path::Path};

use psqs::{
    program::{molpro::Molpro, mopac::Mopac},
    queue::{pbs::Pbs, slurm::Slurm},
};
use rust_pbqff::{
    cleanup,
    config::{self, Config},
    coord_type::{normal::Normal, Cart, CoordType, SIC},
    Intder,
};

include!(concat!(env!("OUT_DIR"), "/version.rs"));

macro_rules! queue {
    ($q: ty, $config: ident) => {
        <$q>::new(
            $config.chunk_size,
            $config.job_limit,
            $config.sleep_int,
            "pts",
        )
    };
}

use clap::Parser;

/// push button quartic force fields
#[derive(Parser, Debug)]
#[command(author, about, long_about = None)]
struct Args {
    /// input file
    #[arg(value_parser, default_value_t = String::from("pbqff.toml"))]
    infile: String,

    /// whether or not to overwrite existing output
    #[arg(short, long, default_value_t = false)]
    overwrite: bool,

    /// print the git commit hash and exit
    #[arg(short, long, default_value_t = false)]
    version: bool,

    /// the maximum number of threads to use by rayon
    #[arg(short, long, default_value_t = 0)]
    threads: usize,
}

fn main() -> Result<(), std::io::Error> {
    let args = Args::parse();
    if args.version {
        println!("version: {}", version());
        return Ok(());
    }
    let path = Path::new("pbqff.out");
    if path.exists() && !args.overwrite {
        eprintln!("existing pbqff output. overwrite with -o/--overwrite");
        return Ok(());
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
    println!("PID: {}", std::process::id());
    println!("version: {}", version());
    psqs::max_threads(args.threads);
    cleanup();
    let _ = std::fs::create_dir("pts");
    let config = Config::load(&args.infile);
    let m = (config.coord_type, config.program, config.queue);
    let (spectro, output) = match m {
        (
            config::CoordType::Normal,
            config::Program::Molpro,
            config::Queue::Pbs,
        ) => <Normal as CoordType<_, _, Molpro>>::run(
            Normal::default(),
            &mut std::io::stdout(),
            &queue!(Pbs, config),
            &config,
        ),
        (
            config::CoordType::Cart,
            config::Program::Mopac,
            config::Queue::Pbs,
        ) => <Cart as CoordType<_, _, Mopac>>::run(
            Cart,
            &mut std::io::stdout(),
            &queue!(Pbs, config),
            &config,
        ),
        (
            config::CoordType::Cart,
            config::Program::Mopac,
            config::Queue::Slurm,
        ) => <Cart as CoordType<_, _, Mopac>>::run(
            Cart,
            &mut std::io::stdout(),
            &queue!(Slurm, config),
            &config,
        ),
        (
            config::CoordType::Cart,
            config::Program::Molpro,
            config::Queue::Pbs,
        ) => <Cart as CoordType<_, _, Molpro>>::run(
            Cart,
            &mut std::io::stdout(),
            &queue!(Pbs, config),
            &config,
        ),
        (
            config::CoordType::Cart,
            config::Program::Molpro,
            config::Queue::Slurm,
        ) => <Cart as CoordType<_, _, Molpro>>::run(
            Cart,
            &mut std::io::stdout(),
            &queue!(Slurm, config),
            &config,
        ),
        (
            config::CoordType::Sic,
            config::Program::Mopac,
            config::Queue::Pbs,
        ) => {
            let sic = SIC::new(Intder::load_file("intder.in"));
            <SIC as CoordType<_, _, Mopac>>::run(
                sic,
                &mut std::io::stdout(),
                &queue!(Pbs, config),
                &config,
            )
        }
        (
            config::CoordType::Sic,
            config::Program::Mopac,
            config::Queue::Slurm,
        ) => {
            let sic = SIC::new(Intder::load_file("intder.in"));
            <SIC as CoordType<_, _, Mopac>>::run(
                sic,
                &mut std::io::stdout(),
                &queue!(Slurm, config),
                &config,
            )
        }
        (
            config::CoordType::Sic,
            config::Program::Molpro,
            config::Queue::Pbs,
        ) => {
            let sic = SIC::new(Intder::load_file("intder.in"));
            <SIC as CoordType<_, _, Molpro>>::run(
                sic,
                &mut std::io::stdout(),
                &queue!(Pbs, config),
                &config,
            )
        }
        (
            config::CoordType::Sic,
            config::Program::Molpro,
            config::Queue::Slurm,
        ) => {
            let sic = SIC::new(Intder::load_file("intder.in"));
            <SIC as CoordType<_, _, Molpro>>::run(
                sic,
                &mut std::io::stdout(),
                &queue!(Slurm, config),
                &config,
            )
        }
        _ => panic!("unsupported combination {m:#?}"),
    };

    spectro.write_output(&mut std::io::stdout(), &output)?;

    let mut f = std::fs::File::create("spectro.json")?;
    use std::io::Write;
    writeln!(f, "{}", serde_json::to_string_pretty(&output)?)?;

    println!("normal termination of pbqff");

    Ok(())
}
