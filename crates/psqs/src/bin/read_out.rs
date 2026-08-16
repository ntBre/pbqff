use psqs::program::{Program, mopac::Mopac};

fn main() {
    let program = Mopac;
    let mut res = Vec::new();
    for _ in 0..1000 {
        res.push(program.read_output("testfiles/job"));
    }
}
