use std::{env::temp_dir, io::Stdout, time::Instant};

use pbqff::{
    config::Config,
    coord_type::{Cart, CoordType},
};
use psqs::{program::mopac::Mopac, queue::local::Local};

fn main() -> Result<(), std::io::Error> {
    let config = Config {
        template: "scfcrt=1.D-21 aux(precision=14 comp xp xs xw) PM6 THREADS=1"
            .to_owned(),
        ..Config::load("testfiles/cart.toml")
    };
    let queue = Local {
        chunk_size: 128,
        dir: "pts".to_string(),
        ..Default::default()
    };

    let dir = temp_dir();
    for d in ["opt", "pts", "freqs"] {
        let dir = dir.join(d);
        let _ = std::fs::remove_dir(&dir);
        let _ = std::fs::create_dir(&dir);
    }

    let now = Instant::now();
    let (spectro, output) = <Cart as CoordType<Stdout, Local, Mopac>>::run(
        Cart,
        dir,
        &mut std::io::stdout(),
        &queue,
        &config,
    );
    println!(
        "Finished Cart run after {} seconds.",
        now.elapsed().as_secs()
    );

    spectro.write_output(&mut std::io::stdout(), &output)?;
    println!("normal termination of pbqff");

    Ok(())
}
