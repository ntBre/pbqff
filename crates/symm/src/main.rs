use std::{
    error::Error,
    io::{stdin, Read},
    str::FromStr,
};

use symm::Molecule;

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args();
    let tol = if args.any(|s| s == "-t") {
        args.next().unwrap().parse().unwrap()
    } else {
        1e-6
    };
    let mut buf = String::new();
    stdin().read_to_string(&mut buf)?;
    let mut mol = Molecule::from_str(&buf)?;
    mol.normalize();
    println!("normalized molecule:\n{mol:.8}");
    let pg = mol.point_group_approx(tol);
    println!("point group: {pg}");
    Ok(())
}
