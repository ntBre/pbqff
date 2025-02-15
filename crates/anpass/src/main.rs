use anpass::Anpass;

fn main() {
    let args: Vec<_> = std::env::args().collect();
    let infile = args.get(1);
    let anpass = match infile {
        Some(s) => Anpass::load_file(s),
        None => Anpass::load(std::io::stdin()),
    };
    let (f9903, bias, res, kind) = anpass.run().unwrap();
    println!("bias: {bias}");
    println!("Sum of squared residuals: {res:12.6e}");
    println!("stationary point is a {kind}");
    let filename = "fort.9903";
    let mut f = match std::fs::File::create(filename) {
        Ok(f) => f,
        Err(e) => panic!("failed to create {filename} with {e}"),
    };
    anpass.write9903(&mut f, &f9903);
}
