use std::{fmt::Display, marker::Sync, path::Path};

pub use intder::IntderError;
use intder::{Intder, Siic};
use na::vector;
use nalgebra as na;
use psqs::{
    geom::Geom,
    program::{Program, Template},
    queue::Queue,
};
use spectro::{Output, Spectro};
use symm::{Irrep, Molecule, Pg, PointGroup};
use taylor::Taylor;

use super::{CoordType, SPECTRO_HEADER};
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

impl<
        W: std::io::Write,
        Q: Queue<P> + Sync,
        P: Program + Clone + Send + Sync,
    > CoordType<W, Q, P> for SIC
{
    fn run(self, w: &mut W, queue: &Q, config: &Config) -> (Spectro, Output) {
        let template = Template::from(&config.template);
        writeln!(w, "{}", config).unwrap();
        // optimize the geometry
        let geom = if config.optimize {
            let res = optimize(
                queue,
                config.geometry.clone(),
                template.clone(),
                config.charge,
            );
            let geom = if let Err(e) = res {
                panic!("optimize failed with {:?}", e);
            } else {
                res.unwrap()
            };
            let geom = Geom::Xyz(geom.cart_geom.unwrap());
            writeln!(w, "Optimized Geometry:\n{}", geom).unwrap();
            geom
        } else {
            // expecting cartesian geometry in angstroms
            assert!(config.geometry.is_xyz());
            config.geometry.clone()
        };

        let mol = {
            let mut mol = Molecule::new(geom.xyz().unwrap().to_vec());
            mol.normalize();
            mol
        };
        let mut pg = mol.point_group_approx(SYMM_EPS);

        // use c2v subgroup for d2h
        if pg.is_d2h() {
            writeln!(w, "full point group is D2h, using C2v subgroup").unwrap();
            pg = pg.subgroup(Pg::C2v).unwrap();
        };

        writeln!(w, "Normalized Geometry:\n{:20.12}", mol).unwrap();
        writeln!(w, "Point Group = {}", pg).unwrap();

        let mut intder = self.intder;
        let (geoms, taylor, taylor_disps, atomic_numbers) =
            generate_pts(w, &mol, &pg, &mut intder, config.step_size).unwrap();

        let dir = "pts/inp";
        let jobs =
            P::build_jobs(&geoms, dir, 0, 1.0, 0, config.charge, template);

        writeln!(w, "\n{} atoms require {} jobs", mol.atoms.len(), jobs.len())
            .unwrap();

        let mut energies = vec![0.0; jobs.len()];
        let time = queue
            .drain(dir, jobs, &mut energies)
            .expect("single-point energies failed");
        eprintln!("total job time: {time} sec");

        let _ = std::fs::create_dir("freqs");
        freqs(
            w,
            "freqs",
            &mut energies,
            &mut intder,
            &taylor,
            &taylor_disps,
            &atomic_numbers,
            config.step_size,
        )
        .unwrap()
    }
}

type AtomicNumbers = Vec<usize>;

pub fn generate_pts<W: std::io::Write>(
    w: &mut W,
    mol: &Molecule,
    pg: &PointGroup,
    intder: &mut Intder,
    step_size: f64,
) -> Result<(Vec<Geom>, Taylor, taylor::Disps, AtomicNumbers), IntderError> {
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

    // LIN1 is the only coordinate with dummies
    let ndum = if intder
        .simple_internals
        .iter()
        .any(|s| matches!(s, Siic::Lin1(..)))
    {
        // default to Z, hopefully this only applies for semp disaster
        intder.add_dummies(pg.axis().unwrap_or_else(|| {
            eprintln!("LIN1 but no axis in point group, using z");
            symm::Axis::Z
        }))
    } else {
        0
    };

    // convert them to Cartesian coordinates
    let disps = intder.convert_disps()?;
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
            Err(e) => {
                eprintln!(
                    "irrep determination failed on coord {i}/{} with {}",
                    disps.len(),
                    e.msg()
                );
                Irrep::A
            }
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
    let checks = Taylor::make_checks(irreps, pg);
    // run taylor.py to get fcs and disps
    let taylor = Taylor::new(5, nsic, checks.0, checks.1);
    let taylor_disps = taylor.disps();
    intder.disps = taylor_disps.to_intder(step_size);

    if DEBUG {
        writeln!(w, "\nIntder Input:\n{}", intder).unwrap();
    }

    // build and run the points using psqs
    // TODO handle error
    let _ = std::fs::create_dir_all("pts/inp");

    write_file("pts/intder.in", &intder).unwrap();

    // these are the displacements that go in file07, but I'll use them from
    // memory to build the jobs
    let file07 = intder.convert_disps()?;

    let mut geoms = Vec::with_capacity(file07.len());
    for geom in file07 {
        // this is a bit unsightly, but I also don't want to duplicate the
        // `from_slices` code in psqs
        let mut mol = Molecule::from_slices(
            atomic_numbers.clone(),
            &geom.as_slice()[..geom.len() - 3 * ndum],
        );
        mol.to_angstrom();
        geoms.push(Geom::from(mol));
    }
    Ok((geoms, taylor, taylor_disps, atomic_numbers))
}

fn write_file(f: impl AsRef<Path>, d: impl Display) -> std::io::Result<()> {
    use std::io::Write;
    let mut f = std::fs::File::create(f)?;
    writeln!(f, "{}", d)
}

/// an Error type containing information about a failure to run `freqs`
#[derive(Debug)]
pub struct FreqError(pub String);

/// run the frequency portion of a QFF in `dir`. The caller is responsible for
/// ensuring this directory exists.
#[allow(clippy::too_many_arguments)]
pub fn freqs<W: std::io::Write>(
    w: &mut W,
    dir: &str,
    energies: &mut [f64],
    intder: &mut Intder,
    taylor: &Taylor,
    taylor_disps: &taylor::Disps,
    atomic_numbers: &AtomicNumbers,
    step_size: f64,
) -> Result<(Spectro, Output), FreqError> {
    let min = energies.iter().cloned().reduce(f64::min).unwrap();
    for energy in energies.iter_mut() {
        *energy -= min;
    }

    // run anpass
    let anpass = Taylor::to_anpass(taylor, taylor_disps, energies, step_size);
    write_file(format!("{dir}/anpass.in"), &anpass).unwrap();
    let (fcs, long_line, res) = if DEBUG {
        writeln!(w, "Anpass Input:\n{}", anpass).unwrap();
        let (fcs, long_line, res) = match anpass.run_debug(w) {
            Ok(v) => v,
            Err(e) => return Err(FreqError(e.0)),
        };
        writeln!(w, "\nStationary Point:\n{}", long_line).unwrap();
        (fcs, long_line, res)
    } else {
        match anpass.run() {
            Ok(v) => v,
            Err(e) => return Err(FreqError(e.0)),
        }
    };

    writeln!(w, "anpass sum of squared residuals: {:17.8e}", res).unwrap();

    // intder_geom
    intder.disps = vec![long_line.disp.as_slice().to_vec()];
    write_file(format!("{dir}/intder_geom.in"), &intder).unwrap();
    let refit_geom = intder.convert_disps().unwrap();
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

    intder.atoms = mol.atoms.iter().map(intder::Atom::from).collect();
    intder.input_options[3] = 4;
    intder.input_options[6] = 2;
    intder.input_options[10] = 3;
    intder.input_options[13] = 0;
    intder.input_options[14] = 0;
    write_file(format!("{dir}/intder.in"), &intder).unwrap();

    let (f2, f3, f4) = intder.convert_fcs();
    Intder::dump_fcs(dir, &f2, &f3, &f4);

    // spectro
    let mut spectro = Spectro::from(mol);
    spectro.header = SPECTRO_HEADER.to_vec();

    let fc3 = spectro::new_fc3(spectro.n3n, &f3);
    let fc4 = spectro::new_fc4(spectro.n3n, &f4);

    let input = format!("{}/spectro.in", dir);
    if DEBUG {
        writeln!(w, "Spectro Input:\n{}", spectro).unwrap();
    }
    spectro.write(&input).unwrap();

    let (output, _) = spectro.run(spectro::Derivative::Quartic(f2, fc3, fc4));

    Ok((spectro, output))
}
