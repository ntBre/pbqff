use std::{fs::File, os::unix::prelude::AsRawFd};

use psqs::queue::slurm::Slurm;
use rust_pbqff::{
    config::{self, Config},
    coord_type::{Cart, CoordType, SIC},
    Intder, Spectro,
};

fn cleanup() {
    let _ = std::fs::remove_dir("opt");
    let _ = std::fs::remove_dir("pts");
    let _ = std::fs::remove_dir("freqs");
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
    let spectro = Spectro::nocurvil();
    let queue = Slurm::new(32, 2048, 2, "pts");
    let output =
        match config.coord_type {
            config::CoordType::cart => {
                Cart.run(&mut std::io::stdout(), &queue, &config, &spectro)
            }
            config::CoordType::sic => SIC::new(Intder::load_file("intder.in"))
                .run(&mut std::io::stdout(), &queue, &config, &spectro),
        };

    spectro.write_output(&mut std::io::stdout(), output)?;
    println!("normal termination of pbqff");

    Ok(())
}
