use intder::Intder;
use symm::Molecule;

fn main() {
    let args: Vec<_> = std::env::args().collect();
    let mut intder = Intder::load_file(&args[1]);
    let labels: Vec<_> = args[2..].iter().map(|s| s.parse().unwrap()).collect();
    assert_eq!(intder.geom.len(), labels.len());
    let mut coords = Vec::new();
    for atom in &intder.geom {
        coords.extend(atom.as_slice());
    }
    let mut mol = Molecule::from_slices(&labels, &coords);
    mol.normalize();
    println!("{mol}");
    let pg = mol.point_group_approx(1e-6);
    println!("Point group = {pg}");
    // load the initial intder
    let nsic = intder.symmetry_internals.len();
    // generate a displacement for each SIC
    let mut disps = Vec::new();
    for i in 0..nsic {
        let mut disp = vec![0.0; nsic];
        disp[i] = 0.005;
        disps.push(disp);
    }
    intder.geom = intder::geom::Geom::from(mol.clone());
    intder.geom.to_bohr();
    intder.disps = disps;
    let disps = intder.convert_disps().unwrap();
    // convert displacements -> symm::Molecules and determine irrep
    let mut irreps = Vec::new();
    for (i, disp) in disps.iter().enumerate() {
        let disp = disp.as_slice();
        let m = Molecule::from_slices(&labels, disp);
        let irrep = match m.irrep_approx(&pg, 1e-6) {
            Ok(rep) => rep,
            Err(e) => {
                eprintln!("failed on coord {} with {}", i, e.msg());
                symm::Irrep::A
            }
        };
        irreps.push((i, irrep));
    }
    println!("\nSymmetry Internal Coordinates:");
    let just_irreps: Vec<_> = irreps.iter().map(|s| s.1).collect();
    intder.print_sics(&mut std::io::stdout(), &just_irreps);
}
