use std::{io::Stdout, time::Instant};

use psqs::{program::mopac::Mopac, queue::local::LocalQueue};
use rust_pbqff::{
    config::Config,
    coord_type::{Cart, CoordType},
};

fn main() -> Result<(), std::io::Error> {
    // inlined cleanup from main
    let _ = std::fs::remove_dir("opt");
    let _ = std::fs::remove_dir("pts");
    let _ = std::fs::remove_dir("freqs");
    // end inline
    let _ = std::fs::create_dir("pts");
    let config = Config::load("testfiles/cart.toml");
    let queue = LocalQueue {
        chunk_size: 128,
        dir: "pts".to_string(),
    };

    let now = Instant::now();
    let (spectro, output) = <Cart as CoordType<Stdout, LocalQueue, Mopac>>::run(
        &Cart,
        &mut std::io::stdout(),
        &queue,
        &config,
    );
    println!(
        "Finished Cart run after {} seconds.",
        now.elapsed().as_secs()
    );

    spectro.write_output(&mut std::io::stdout(), output)?;
    println!("normal termination of pbqff");

    Ok(())
}
