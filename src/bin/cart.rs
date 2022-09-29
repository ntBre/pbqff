use std::time::Instant;

use psqs::queue::local::LocalQueue;
use rust_pbqff::{
    config::Config,
    coord_type::{Cart, CoordType},
    Spectro,
};

fn main() -> Result<(), std::io::Error> {
    // inlined cleanup from main
    let _ = std::fs::remove_dir("opt");
    let _ = std::fs::remove_dir("pts");
    let _ = std::fs::remove_dir("freqs");
    // end inline
    let _ = std::fs::create_dir("pts");
    let config = Config::load("testfiles/cart.toml");
    let spectro = Spectro::nocurvil();
    let queue = LocalQueue {
        chunk_size: 128,
        dir: "pts".to_string(),
    };

    let now = Instant::now();
    let output = Cart.run(&mut std::io::stdout(), &queue, &config, &spectro);
    println!(
        "Finished Cart run after {} seconds.",
        now.elapsed().as_secs()
    );

    spectro.write_output(&mut std::io::stdout(), output)?;
    println!("normal termination of pbqff");

    Ok(())
}
