use psqs::program::{mopac::Mopac, Program};

fn main() {
    let mut res = Vec::new();
    for _ in 0..1000 {
        res.push(Mopac::read_output("testfiles/job"));
    }
}
