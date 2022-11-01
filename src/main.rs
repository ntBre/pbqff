use std::{fs::File, os::unix::prelude::AsRawFd};

use psqs::{
    program::{molpro::Molpro, mopac::Mopac},
    queue::{pbs::Pbs, slurm::Slurm},
};
use rust_pbqff::{
    cleanup,
    config::{self, Config},
    coord_type::{Cart, CoordType, SIC},
    Intder,
};

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
    println!("PID: {}", std::process::id());
    cleanup();
    let _ = std::fs::create_dir("pts");
    let config = Config::load("pbqff.toml");
    let (spectro, output) =
        match (config.coord_type, config.program, config.queue) {
            (
                config::CoordType::Cart,
                config::Program::Mopac,
                config::Queue::Pbs,
            ) => <Cart as CoordType<_, _, Mopac>>::run(
                &Cart,
                &mut std::io::stdout(),
                &queue!(Pbs, config),
                &config,
            ),
            (
                config::CoordType::Cart,
                config::Program::Mopac,
                config::Queue::Slurm,
            ) => <Cart as CoordType<_, _, Mopac>>::run(
                &Cart,
                &mut std::io::stdout(),
                &queue!(Slurm, config),
                &config,
            ),
            (
                config::CoordType::Cart,
                config::Program::Molpro,
                config::Queue::Pbs,
            ) => <Cart as CoordType<_, _, Molpro>>::run(
                &Cart,
                &mut std::io::stdout(),
                &queue!(Pbs, config),
                &config,
            ),
            (
                config::CoordType::Cart,
                config::Program::Molpro,
                config::Queue::Slurm,
            ) => <Cart as CoordType<_, _, Molpro>>::run(
                &Cart,
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
                    &sic,
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
                    &sic,
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
                    &sic,
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
                    &sic,
                    &mut std::io::stdout(),
                    &queue!(Slurm, config),
                    &config,
                )
            }
        };

    spectro.write_output(&mut std::io::stdout(), &output)?;

    let mut f = std::fs::File::create("spectro.json")?;
    use std::io::Write;
    writeln!(f, "{}", serde_json::to_string_pretty(&output)?)?;

    println!("normal termination of pbqff");

    Ok(())
}
