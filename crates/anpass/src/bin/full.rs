use anpass::{Anpass, Bias};

fn main() {
    let anpass = Anpass::load_file("testfiles/c3h2.in");
    // initial fitting
    let (coeffs, _) = anpass.fit();
    // find stationary point
    let (x, _) = anpass.newton(&coeffs).unwrap();
    // determine energy at stationary point
    let e = anpass.eval(&x, &coeffs);
    // bias the displacements and energies to the new stationary point
    let anpass = anpass.bias(&Bias { disp: x, energy: e });
    // perform the refitting
    let (coeffs, _) = anpass.fit();
    for c in &coeffs {
        println!("{c}");
    }
}
