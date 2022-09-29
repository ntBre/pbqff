use std::{path::Path, rc::Rc};

use intder::Intder;
use na::vector;
use nalgebra as na;
use psqs::{
    geom::Geom,
    program::{mopac::Mopac, Template},
    queue::Queue,
};
use rust_anpass::Anpass;
use spectro::{Output, Spectro};
use symm::{Irrep, Molecule, PointGroup};
use taylor::{Checks, Taylor};

use super::CoordType;
use crate::{config::Config, optimize};

/// whether or not to print the input files used for intder, anpass, and spectro
static DEBUG: bool = false;

/// the precision to call things symmetric
const SYMM_EPS: f64 = 1e-6;

pub struct SIC {
    intder: Intder,
}

impl SIC {
    pub fn new(intder: Intder) -> Self {
        Self { intder }
    }
}

impl<W: std::io::Write, Q: Queue<Mopac>> CoordType<W, Q> for SIC {
    fn run(
        &self,
        w: &mut W,
        queue: &Q,
        config: &Config,
        spectro: &Spectro,
    ) -> Output {
        let template = Template::from(&config.template);
        writeln!(w, "{}", config).unwrap();
        // optimize the geometry
        let geom = if config.optimize {
            let geom = optimize(
                queue,
                config.geometry.clone(),
                template.clone(),
                config.charge,
            )
            .expect("optimization failed");
            let geom = Geom::Xyz(geom.cart_geom);
            writeln!(w, "Optimized Geometry:\n{}", geom).unwrap();
            geom
        } else {
            todo!();
        };

        let mol = {
            let mut mol = Molecule::new(geom.xyz().unwrap().to_vec());
            mol.normalize();
            mol
        };
        let pg = mol.point_group_approx(SYMM_EPS);

        writeln!(w, "Normalized Geometry:\n{:20.12}", mol).unwrap();
        writeln!(w, "Point Group = {}", pg).unwrap();

        let mut intder = self.intder.clone();
        let (geoms, taylor, taylor_disps, atomic_numbers) =
            generate_pts(w, &mol, &pg, &mut intder, config.step_size, &vec![]);

        // TODO switch on Program type eventually

        let mut jobs = Mopac::build_jobs(
            &geoms,
            None,
            "pts",
            0,
            1.0,
            0,
            config.charge,
            template,
        );

        writeln!(w, "\n{} atoms require {} jobs", mol.atoms.len(), jobs.len())
            .unwrap();

        let mut energies = vec![0.0; jobs.len()];
        queue
            .drain(&mut jobs, &mut energies)
            .expect("single-point energies failed");

        let _ = std::fs::create_dir("freqs");
        freqs(
            w,
            "freqs",
            &mut energies,
            &mut intder,
            &taylor,
            &taylor_disps,
            &atomic_numbers,
            spectro,
            config.step_size,
        )
    }
}

/// generate the Taylor series mod and equivalence checks from `irreps` in `pg`
fn make_taylor_checks(
    irreps: Vec<(usize, Irrep)>,
    pg: &PointGroup,
) -> (Option<Checks>, Option<Checks>) {
    use symm::Irrep::*;
    use symm::PointGroup::*;
    match pg {
        C1 => (None, None),
        C2 { axis: _ } => {
            todo!();
        }
        Cs { plane: _ } => {
            // only A'' modes go in checks[0], other two checks are 0-0
            let mut checks = Checks::new();
            for i in irreps {
                match i.1 {
                    Ap => (),
                    App => {
                        if checks[(0, 0)] == 0 {
                            checks[(0, 0)] = i.0 + 1;
                            checks[(0, 1)] = i.0 + 1;
                        } else if i.0 + 1 > checks[(0, 1)] {
                            checks[(0, 1)] = i.0 + 1;
                        }
                    }
                    _ => panic!("non-Cs irrep found in Cs point group"),
                }
            }
            (Some(checks.clone()), Some(checks))
        }
        C2v { axis: _, planes: _ } => {
            let mut checks = Checks::new();
            // first one you hit goes in checks.0, second goes in checks.1
            for i in irreps {
                match i.1 {
                    A1 => (),
                    B1 => {
                        if checks[(0, 0)] == 0 {
                            checks[(0, 0)] = i.0 + 1;
                            checks[(0, 1)] = i.0 + 1;
                        } else if i.0 + 1 > checks[(0, 1)] {
                            checks[(0, 1)] = i.0 + 1;
                        }
                    }
                    B2 => {
                        if checks[(1, 0)] == 0 {
                            checks[(1, 0)] = i.0 + 1;
                            checks[(1, 1)] = i.0 + 1;
                        } else if i.0 + 1 > checks[(1, 1)] {
                            checks[(1, 1)] = i.0 + 1;
                        }
                    }
                    A2 => {
                        if checks[(2, 0)] == 0 {
                            checks[(2, 0)] = i.0 + 1;
                            checks[(2, 1)] = i.0 + 1;
                        } else if i.0 + 1 > checks[(2, 1)] {
                            checks[(2, 1)] = i.0 + 1;
                        }
                    }
                    _ => panic!("non-C2v irrep found in C2v point group"),
                }
            }
            (Some(checks.clone()), Some(checks))
        }
        D2h { axes: _, planes: _ } => todo!(),
        C3v { axis: _, plane: _ } => todo!(),
    }
}

fn taylor_to_anpass(
    taylor: &Taylor,
    taylor_disps: &Vec<Vec<isize>>,
    energies: &[f64],
    step_size: f64,
) -> Anpass {
    let mut disps = Vec::new();
    for disp in taylor_disps {
        for coord in disp {
            disps.push(*coord as f64 * step_size);
        }
    }
    let tdl = taylor_disps.len();
    let fl = taylor.forces.len();
    let mut fs = Vec::new();
    for row in &taylor.forces {
        for c in row {
            fs.push(*c as i32);
        }
    }
    Anpass {
        disps: na::DMatrix::from_row_slice(tdl, disps.len() / tdl, &disps),
        energies: na::DVector::from_row_slice(energies),
        exponents: na::DMatrix::from_column_slice(
            taylor.forces[0].len(),
            fl,
            &fs,
        ),
        bias: None,
    }
}

fn disp_to_intder(disps: &Vec<Vec<isize>>, step_size: f64) -> Vec<Vec<f64>> {
    let mut ret = Vec::new();
    for disp in disps {
        let disp: Vec<_> = disp.iter().map(|i| *i as f64 * step_size).collect();
        ret.push(disp);
    }
    ret
}

type TaylorDisps = Vec<Vec<isize>>;
type AtomicNumbers = Vec<usize>;

// TODO clean up this dummy atom interface. this is the difficulty of not
// reading a template geometry. Somehow I have to add dummy atoms when they are
// needed *after* the optimization, which obviously doesn't include them. Right
// now I'm specifying the dummy atoms in a pretty terrible format in the
// rust-semp config file and they aren't really used within this package
pub fn generate_pts<W: std::io::Write>(
    w: &mut W,
    mol: &Molecule,
    pg: &PointGroup,
    intder: &mut Intder,
    step_size: f64,
    dummies: &Vec<(usize, usize)>,
) -> (Vec<Rc<Geom>>, Taylor, TaylorDisps, AtomicNumbers) {
    let atomic_numbers = mol.atomic_numbers();

    // load the initial intder
    let nsic = intder.symmetry_internals.len();
    // generate a displacement for each SIC
    let mut disps = Vec::new();
    for i in 0..nsic {
        let mut disp = vec![0.0; nsic];
        disp[i] = step_size;
        disps.push(disp);
    }
    intder.geom = intder::geom::Geom::from(mol.clone());
    intder.geom.to_bohr();
    intder.disps = disps;
    // add the dummy atoms
    let mut ndum = 0;
    const ZERO: f64 = 1e-8;
    for dummy in dummies {
        // atom the dummy attaches to
        let real_coord = intder.geom[dummy.1];
        // two of them should be zero and one is non-zero
        let mut zeros = vec![];
        let mut nonzero = 0;
        for (i, c) in real_coord.iter().enumerate() {
            if c.abs() < ZERO {
                zeros.push(i);
            } else {
                nonzero = i;
            }
        }
        if zeros.len() != 2 {
            dbg!(real_coord);
            dbg!(zeros);
            dbg!(nonzero);
            panic!("dummy atom for non-linear molecule");
        }

        let mut coord = [0.0; 3];
        // match the nonzero field in the real geometry
        coord[nonzero] = intder.geom[dummy.1][nonzero];
        coord[zeros[0]] = 1.1111111111;
        coord[zeros[1]] = 0.0;
        intder.geom.push(na::Vector3::from(coord));

        let mut coord = [0.0; 3];
        // match the nonzero field in the real geometry
        coord[nonzero] = intder.geom[dummy.1][nonzero];
        coord[zeros[1]] = 1.1111111111;
        coord[zeros[0]] = 0.0;
        intder.geom.push(na::Vector3::from(coord));

        // push dummy atoms perpendicular in both directions

        ndum += 2;
    }
    // convert them to Cartesian coordinates
    let disps = intder.convert_disps();
    // convert displacements -> symm::Molecules and determine irrep
    let mut irreps = Vec::new();
    for (i, disp) in disps.iter().enumerate() {
        let disp = disp.as_slice();
        let m = Molecule::from_slices(
            atomic_numbers.clone(),
            &disp[..disp.len() - 3 * ndum],
        );
        let irrep = match m.irrep_approx(pg, SYMM_EPS) {
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

    writeln!(w, "\nSymmetry Internal Coordinates:").unwrap();
    intder.print_sics(w, &just_irreps);

    // generate checks
    let checks = make_taylor_checks(irreps, pg);
    // run taylor.py to get fcs and disps
    let taylor = Taylor::new(5, nsic, checks.0, checks.1);
    let taylor_disps = taylor.disps();
    intder.disps = disp_to_intder(&taylor_disps, step_size);

    if DEBUG {
        writeln!(w, "\nIntder Input:\n{}", intder).unwrap();
    }

    // these are the displacements that go in file07, but I'll use them from
    // memory to build the jobs
    let file07 = intder.convert_disps();

    // build and run the points using psqs
    // TODO handle error
    let _ = std::fs::create_dir("pts");
    let mut geoms = Vec::with_capacity(file07.len());
    for geom in file07 {
        // this is a bit unsightly, but I also don't want to duplicate the
        // `from_slices` code in psqs
        geoms.push(Rc::new(Geom::from(Molecule::from_slices(
            atomic_numbers.clone(),
            &geom.as_slice()[..geom.len() - 3 * ndum],
        ))));
    }
    (geoms, taylor, taylor_disps, atomic_numbers)
}

/// run the frequency portion of a QFF in `dir`. The caller is responsible for
/// ensuring this directory exists.
#[allow(clippy::too_many_arguments)]
pub fn freqs<W: std::io::Write>(
    w: &mut W,
    dir: &str,
    energies: &mut [f64],
    intder: &mut Intder,
    taylor: &Taylor,
    taylor_disps: &TaylorDisps,
    atomic_numbers: &AtomicNumbers,
    spectro: &Spectro,
    step_size: f64,
) -> Output {
    let min = energies.iter().cloned().reduce(f64::min).unwrap();
    for energy in energies.iter_mut() {
        *energy -= min;
    }

    // run anpass
    let anpass = taylor_to_anpass(taylor, taylor_disps, energies, step_size);
    let (fcs, long_line) = if DEBUG {
        writeln!(w, "Anpass Input:\n{}", anpass).unwrap();
        let (fcs, long_line) = anpass.run_debug(w);
        writeln!(w, "\nStationary Point:\n{}", long_line).unwrap();
        (fcs, long_line)
    } else {
        anpass.run()
    };

    // intder_geom
    intder.disps = vec![long_line.disp.as_slice().to_vec()];
    let refit_geom = intder.convert_disps();
    let refit_geom = refit_geom[0].as_slice();
    let l = refit_geom.len() - 3 * intder.ndum();
    let dummies = &refit_geom[l..];
    let mol = Molecule::from_slices(atomic_numbers.clone(), &refit_geom[..l]);

    writeln!(w, "\nRefit Geometry\n{:20.12}", mol).unwrap();

    intder.geom = intder::geom::Geom::from(mol.clone());
    for dummy in dummies.chunks(3) {
        intder.geom.push(vector![dummy[0], dummy[1], dummy[2]]);
    }

    intder.disps = vec![];

    // intder freqs
    for fc in fcs {
        // skip zeroth and first derivatives
        if (fc.1, fc.2, fc.3) != (0, 0, 0) {
            intder.add_fc(vec![fc.0, fc.1, fc.2, fc.3], fc.4);
        }
    }

    if DEBUG {
        writeln!(w, "Intder Freqs Input\n{}", intder).unwrap();
    }

    let (f2, f3, f4) = intder.convert_fcs();
    Intder::dump_fcs(dir, &f2, &f3, &f4);

    // spectro
    let mut spectro = spectro.clone();
    spectro.geom = mol;
    let input = format!("{}/spectro.in", dir);
    if DEBUG {
        writeln!(w, "Spectro Input:\n{}", spectro).unwrap();
    }
    spectro.write(&input).unwrap();

    let dir = Path::new(dir);
    let (output, _) = spectro.run(
        dir.join("fort.15"),
        dir.join("fort.30"),
        dir.join("fort.40"),
    );

    output
}
