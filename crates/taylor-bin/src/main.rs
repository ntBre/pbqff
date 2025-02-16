use clap::Parser;
use intder::Intder;
use symm::{Atom, Molecule, Pg};
use taylor::Taylor;

// borrowed from summarize-bin
fn irrep(ir: &symm::Irrep) -> &'static str {
    match ir {
        symm::Irrep::A => "a",
        symm::Irrep::B => "b",
        symm::Irrep::Ap => "a'",
        symm::Irrep::App => "a''",
        symm::Irrep::A1 => "a_1",
        symm::Irrep::B2 => "b_2",
        symm::Irrep::B1 => "b_1",
        symm::Irrep::A2 => "a_2",
        symm::Irrep::Ag => "a_g",
        symm::Irrep::B1g => "b_{1g}",
        symm::Irrep::B2g => "b_{2g}",
        symm::Irrep::B3g => "b_{3g}",
        symm::Irrep::Au => "a_u",
        symm::Irrep::B1u => "b_{1u}",
        symm::Irrep::B2u => "b_{2u}",
        symm::Irrep::B3u => "b_{3u}",
        symm::Irrep::A1p => "a_1'",
        symm::Irrep::A2p => "a_2'",
        symm::Irrep::Ep => "e'",
        symm::Irrep::A1pp => "a_1''",
        symm::Irrep::A2pp => "a_2''",
        symm::Irrep::Epp => "e''",
        symm::Irrep::E1 => todo!(),
        symm::Irrep::E2 => todo!(),
        symm::Irrep::Bg => todo!(),
        symm::Irrep::Bu => todo!(),
        symm::Irrep::E1p => todo!(),
        symm::Irrep::E2p => todo!(),
        symm::Irrep::E => todo!(),
    }
}

/// generate intder and anpass input files with a taylor series expansion
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// intder input file header with internal coordinates and geometry
    #[arg(value_parser)]
    infile: String,

    /// write intder and anpass input files using the SICs
    #[arg(short, long, default_value_t = false)]
    write: bool,

    /// tolerance to use for symmetry equivalence
    #[arg(short, long, default_value_t = 1e-6)]
    eps: f64,

    /// step size
    #[arg(short, long, default_value_t = 0.005)]
    step_size: f64,

    /// print the SICs in LaTeX for papers
    #[arg(short, long, default_value_t = false)]
    tex: bool,
}

// this is pieced together from parts of pbqff, but it's not clear how to reuse
// any of the parts
fn main() -> std::io::Result<()> {
    let cfg = Args::parse();
    // expects an Intder without dummy atoms
    let mut intder = Intder::load_file(&cfg.infile);
    let pairs = intder.geom.0.iter().zip(&intder.atoms);
    let mut atoms = Vec::new();
    for (g, a) in pairs {
        atoms.push(Atom::new_from_label(&a.label, g[0], g[1], g[2]));
    }
    let mol = {
        let mut mol = Molecule::new(atoms);
        // NOTE could actually check for linearity here by returning rotor type
        // from normalize
        mol.normalize();
        mol
    };
    let pg = {
        let mut pg = mol.point_group_approx(cfg.eps);
        if cfg.write && pg.is_d2h() {
            eprintln!(
                "full point group is D2h, using C2v subgroup. \
		       rerun without -w to print full symmetry"
            );
            pg = pg.subgroup(Pg::C2v).unwrap();
        };
        pg
    };

    println!("Normalized Geometry:\n{:20.12}", mol);
    println!("Point Group = {}", pg);

    let nsic = intder.symmetry_internals.len();
    // generate a displacement for each SIC
    let mut disps = Vec::new();
    for i in 0..nsic {
        let mut disp = vec![0.0; nsic];
        disp[i] = cfg.step_size;
        disps.push(disp);
    }
    intder.disps = disps;
    intder.geom = mol.clone().into();
    let ndum = if let Some(axis) = pg.axis() {
        intder.add_dummies(axis)
    } else {
        0
    };
    let disps = intder.convert_disps().unwrap();

    let atomic_numbers = mol.atomic_numbers();
    let mut irreps = Vec::new();
    for (i, disp) in disps.iter().enumerate() {
        let disp = disp.as_slice();
        let m = Molecule::from_slices(
            &atomic_numbers,
            &disp[..disp.len() - 3 * ndum],
        );
        let irrep = match m.irrep_approx(&pg, cfg.eps) {
            Ok(rep) => rep,
            Err(e) => panic!("failed on coord {} with {}", i, e.msg()),
        };
        irreps.push((i, irrep));
    }
    // sort by irrep symmetry
    irreps.sort_by_key(|k| k.1);

    let just_irreps: Vec<_> = irreps.iter().map(|s| s.1).collect();

    let mut new_sics = Vec::new();
    for irrep in &irreps {
        new_sics.push(intder.symmetry_internals[irrep.0].clone());
    }
    intder.symmetry_internals = new_sics;

    println!("\nSymmetry Internal Coordinates:");
    if cfg.tex {
        println!(r"\begin{{align}}");
        let nsic = intder.symmetry_internals.len();
        assert_eq!(intder.symmetry_internals.len(), just_irreps.len());
        for (i, sic) in intder.symmetry_internals.iter().enumerate() {
            let len = sic.iter().filter(|&&s| s != 0.0).count();
            let (frac, close) = match len {
                1 => ("", ""),
                2 => (r"\frac{1}{\sqrt{2}}[", "]"),
                _ => panic!("unrecognized number of simple internals, {len}"),
            };
            print!(
                "S_{{{:<2}}}({}) &= & {}",
                i + 1,
                irrep(&just_irreps[i]),
                frac
            );
            // number of siics printed so far
            let mut nprt = 0;
            for (j, s) in sic.iter().enumerate() {
                if *s != 0.0 {
                    if nprt > 0 {
                        let sign = match s.signum() as isize {
                            -1 => "-",
                            1 => "+",
                            _ => panic!("it's NaN"),
                        };
                        print!(" {sign} ");
                    }
                    let l = &intder.simple_internals[j];
                    print!(
                        "{}",
                        match &l {
                            intder::Siic::Stretch(a, b) => format!(
                                "r(\\text{{{}}}_{}-\\text{{{}}}_{})",
                                intder.atoms[*a].label,
                                a + 1,
                                intder.atoms[*b].label,
                                b + 1
                            ),
                            intder::Siic::Bend(a, b, c) =>
				format!(
                                "\\angle(\\text{{{}}}_{}-\\text{{{}}}_{}-\\text{{{}}}_{})",
                                intder.atoms[*a].label,
                                a + 1,
                                intder.atoms[*b].label,
                                b + 1,
                                intder.atoms[*c].label,
                                c + 1
                            ),
                            intder::Siic::Torsion(a, b, c, d) =>
				format!(
                                    "\\tau(\\text{{{}}}_{}-\\text{{{}}}_{}\
				     -\\text{{{}}}_{}-\\text{{{}}}_{})",
                                intder.atoms[*a].label,
                                a + 1,
                                intder.atoms[*b].label,
                                b + 1,
                                intder.atoms[*c].label,
                                c + 1,
                                intder.atoms[*d].label,
                                d + 1
                            ),
			    _ => panic!("tell brent to support {}", l),
                            // intder::Siic::Lin1(_, _, _, _) => todo!(),
                            // intder::Siic::Out(_, _, _, _) => todo!(),
                        }
                    );
                    nprt += 1;
                }
            }
            println!("{close}{}", if i < nsic - 1 { "\\\\" } else { "" });
        }
        println!(r"\end{{align}}");
    } else {
        intder.print_sics(&mut std::io::stdout(), &just_irreps);
    }

    if cfg.write {
        // generate checks
        let checks = Taylor::make_checks(irreps, &pg);
        // run taylor.py to get fcs and disps
        let taylor = Taylor::new(5, nsic, checks.0, checks.1);
        let taylor_disps = taylor.disps();

        intder.disps = taylor_disps.to_intder(cfg.step_size);

        let mut f = std::fs::File::create("intder.in")?;
        use std::io::Write;
        writeln!(f, "{}", intder)?;

        let anpass = taylor.to_anpass(
            &taylor_disps,
            &vec![0.0; taylor_disps.len()],
            cfg.step_size,
        );
        let mut f = std::fs::File::create("anpass.in")?;
        writeln!(f, "{}", anpass)?;
    }

    Ok(())
}
