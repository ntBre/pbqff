use psqs::queue::slurm::Slurm;
use rust_pbqff::{
    config::Config,
    coord_type::{CoordType, SIC},
    Intder, Spectro,
};

fn main() {
    let _ = std::fs::create_dir("pts");
    let config = Config::load("pbqff.toml");
    let coord = SIC::new(Intder::load_file("intder.in"));
    let spectro = Spectro::load("spectro.in");
    let queue = Slurm::new(32, 2048, 2, "pts");
    coord.run(&mut std::io::stdout(), &queue, &config, &spectro);
}
