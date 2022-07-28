use psqs::queue::local::LocalQueue;
use rust_pbqff::{
    config::Config,
    coord_type::{Cart, CoordType},
    Spectro,
};
use summarize::Summary;

fn main() {
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
    Cart.run(&mut std::io::stdout(), &queue, &config, &spectro);
    let summ = Summary::new("freqs/spectro2.out");
    println!("\nvibrational frequencies:\n{}", summ);
    println!("normal termination of pbqff");
}
