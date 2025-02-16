use std::fs::File;
use std::io::Write;

use intder::Intder;

fn main() {
    let intder = Intder::load_file("testfiles/c3h2.full");
    let new_carts = intder.convert_disps().unwrap();
    let mut file07 = File::create("file07").expect("failed to create file07");
    for cart in new_carts {
        writeln!(file07, "# GEOMUP #################").unwrap();
        Intder::print_cart(&mut file07, &cart);
    }
}
