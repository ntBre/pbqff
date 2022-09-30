use std::{fs::File, io::Stdout, os::unix::prelude::AsRawFd};

use psqs::{
    program::{molpro::Molpro, mopac::Mopac},
    queue::slurm::Slurm,
};
use rust_pbqff::{
    config::{self, Config},
    coord_type::{Cart, CoordType, SIC},
    Intder,
};

/// clean up from a previous run, emitting warnings on failures
fn cleanup() {
    std::fs::remove_dir_all("opt")
        .unwrap_or_else(|e| eprintln!("failed to remove 'opt' with {e}"));
    std::fs::remove_dir_all("pts")
        .unwrap_or_else(|e| eprintln!("failed to remove 'pts' with {e}"));
    std::fs::remove_dir_all("freqs")
        .unwrap_or_else(|e| eprintln!("failed to remove 'freqs' with {e}"));
}

fn main() -> Result<(), std::io::Error> {
    let outfile = File::create("pbqff.out").expect("failed to create outfile");
    let logfile = File::create("pbqff.log").expect("failed to create log file");
    let out_fd = outfile.as_raw_fd();
    let log_fd = logfile.as_raw_fd();
    // redirect stdout to outfile and stderr to logfile
    unsafe {
        libc::dup2(out_fd, 1);
        libc::dup2(log_fd, 2);
    }
    cleanup();
    let _ = std::fs::create_dir("pts");
    let config = Config::load("pbqff.toml");
    let queue: Slurm = Slurm::new(
        config.chunk_size,
        config.job_limit,
        config.sleep_int,
        "pts",
    );
    let (spectro, output) = match (config.coord_type, config.program) {
        (config::CoordType::Cart, config::Program::Mopac) => {
            <Cart as CoordType<Stdout, Slurm, Mopac>>::run(
                &Cart,
                &mut std::io::stdout(),
                &queue,
                &config,
            )
        }
        (config::CoordType::Sic, config::Program::Mopac) => {
            let sic = SIC::new(Intder::load_file("intder.in"));
            <SIC as CoordType<Stdout, Slurm, Mopac>>::run(
                &sic,
                &mut std::io::stdout(),
                &queue,
                &config,
            )
        }
        (config::CoordType::Cart, config::Program::Molpro) => {
            <Cart as CoordType<Stdout, Slurm, Molpro>>::run(
                &Cart,
                &mut std::io::stdout(),
                &queue,
                &config,
            )
        }
        (config::CoordType::Sic, config::Program::Molpro) => {
            let sic = SIC::new(Intder::load_file("intder.in"));
            <SIC as CoordType<Stdout, Slurm, Molpro>>::run(
                &sic,
                &mut std::io::stdout(),
                &queue,
                &config,
            )
        }
    };

    spectro.write_output(&mut std::io::stdout(), output)?;
    println!("normal termination of pbqff");

    Ok(())
}
