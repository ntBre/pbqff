//! Symmetry-internal coordinate (SIC) QFFs.

use std::{fmt::Display, io::Write, marker::Sync, option::Option, path::Path};

pub use intder::IntderError;
use intder::{Intder, Siic};
use na::vector;
use nalgebra as na;
use psqs::{
    geom::Geom,
    program::{Program, Template},
    queue::Queue,
};
use serde::{Deserialize, Serialize};
use spectro::{Output, Spectro};
use symm::{Irrep, Molecule, Pg, PointGroup};
use taylor::Taylor;

use super::{
    fitted::{AtomicNumbers, Fitted},
    CoordType, Load, SPECTRO_HEADER,
};
use crate::{config::Config, coord_type::CHK_NAME, make_check, optimize};

/// whether or not to print the input files used for intder, anpass, and spectro
pub(crate) static DEBUG: bool = false;

/// the precision to call things symmetric
const SYMM_EPS: f64 = 1e-6;

/// The main type for running SIC QFFs.
pub struct Sic {
    /// An [Intder] used for generating the displacements
    pub intder: Intder,
}

impl Sic {
    /// Construct a [Sic] from an [Intder]
    pub fn new(intder: Intder) -> Self {
        Self { intder }
    }
}

impl<W, Q, P> CoordType<W, Q, P> for Sic
where
    W: Write,
    Q: Queue<P> + Sync,
    P: Program + Clone + Send + Sync + Serialize + for<'a> Deserialize<'a>,
{
    fn run(
        mut self,
        dir: impl AsRef<Path>,
        w: &mut W,
        queue: &Q,
        config: &Config,
    ) -> (Spectro, Output) {
        let template = Template::from(&config.template);
        writeln!(w, "{config}").unwrap();
        // optimize the geometry
        let geom = if config.optimize {
            let res = optimize(
                &dir,
                queue,
                config.geometry.clone(),
                template.clone(),
                config.charge,
            );
            let geom = if let Err(e) = res {
                panic!("optimize failed with {e:?}");
            } else {
                res.unwrap()
            };
            let geom = Geom::Xyz(geom.cart_geom.unwrap());
            writeln!(w, "Optimized Geometry:\n{geom}").unwrap();
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
            writeln!(w, "warning: full point group is D2h, using C2v subgroup")
                .unwrap();
            pg = pg.subgroup(Pg::C2v).unwrap();
        };

        writeln!(w, "Normalized Geometry:\n{mol:20.12}").unwrap();
        writeln!(w, "Point Group = {pg}").unwrap();

        let (geoms, taylor, atomic_numbers) = self
            .generate_pts(&dir, w, &mol, &pg, config.step_size)
            .unwrap();

        let freqs_dir = dir.as_ref().join("freqs");

        let pts_dir = dir.as_ref().join("pts").join("inp");
        let pts_dir = pts_dir.to_string_lossy().to_string();
        let jobs =
            P::build_jobs(geoms, &pts_dir, 0, 1.0, 0, config.charge, template);

        let resume = Resume::new(
            self.intder.clone(),
            taylor,
            atomic_numbers,
            config.step_size,
        );
        resume.dump(dir.as_ref().join(CHK_NAME));

        let (mut energies, time) = queue
            .drain(&pts_dir, jobs, make_check(config.check_int, &dir))
            .expect("single-point energies failed");
        eprintln!("total job time: {time:.1} sec");

        let _ = std::fs::create_dir(&freqs_dir);
        self.freqs(
            w,
            Some(freqs_dir),
            &mut energies,
            &resume.taylor,
            &resume.atomic_numbers,
            resume.step_size,
        )
        .unwrap()
    }

    type Resume = Resume;

    fn resume(
        mut self,
        dir: impl AsRef<Path>,
        w: &mut W,
        queue: &Q,
        config: &Config,
        Resume {
            intder,
            taylor,
            atomic_numbers,
            step_size,
        }: Resume,
    ) -> (Spectro, Output) {
        let freqs_dir = dir.as_ref().join("freqs");
        let dir = dir.as_ref().join("pts").join("inp");
        let _ = std::fs::create_dir_all(&dir);
        let dir = dir.to_str().unwrap().to_owned();
        let chk = format!("{dir}/chk.json");
        let (mut energies, time) = queue
            .resume(&dir, &chk, make_check(config.check_int, &dir))
            .expect("single-point energies failed");
        eprintln!("total job time: {time:.1} sec");

        let _ = std::fs::create_dir(&freqs_dir);
        self.intder = intder;
        self.freqs(
            w,
            Some(freqs_dir),
            &mut energies,
            &taylor,
            &atomic_numbers,
            step_size,
        )
        .unwrap()
    }
}

#[derive(Serialize, Deserialize)]
pub struct Resume {
    intder: Intder,
    taylor: Taylor,
    atomic_numbers: Vec<usize>,
    step_size: f64,
}

impl Resume {
    pub fn new(
        intder: Intder,
        taylor: Taylor,
        atomic_numbers: Vec<usize>,
        step_size: f64,
    ) -> Self {
        Self {
            intder,
            taylor,
            atomic_numbers,
            step_size,
        }
    }
}

impl Load for Resume {}

/// returned by [Sic::prepare_points]. see its documentation for details
struct Prep {
    atomic_numbers: Vec<usize>,

    /// number of symmetry internal coordinates
    nsic: usize,

    /// number of dummy atoms in the geometry
    ndum: usize,

    /// indices of SICs and their corresponding irreps
    irreps: Vec<(usize, Irrep)>,
}

impl Fitted for Sic {
    type Error = IntderError;

    fn generate_pts<W: Write>(
        &mut self,
        dir: impl AsRef<Path>,
        w: &mut W,
        mol: &Molecule,
        pg: &PointGroup,
        step_size: f64,
    ) -> Result<(Vec<Geom>, Taylor, AtomicNumbers), IntderError> {
        let Prep {
            atomic_numbers,
            nsic,
            ndum,
            irreps,
        } = self.prepare_points(mol, step_size, pg, w)?;

        // TODO think we can start here with normals and reuse basically all of
        // the rest of the code

        // generate checks
        let checks = Taylor::make_checks(irreps, pg);
        // run taylor.py to get fcs and disps
        let taylor = Taylor::new(5, nsic, checks.0, checks.1);
        let taylor_disps = taylor.disps();
        self.intder.disps = taylor_disps.to_intder(step_size);

        if DEBUG {
            writeln!(w, "\nIntder Input:\n{}", self.intder).unwrap();
        }

        // build and run the points using psqs
        // TODO handle error
        let _ = std::fs::create_dir_all(dir.as_ref().join("pts/inp"));

        write_file(dir.as_ref().join("pts/intder.in"), &self.intder).unwrap();

        // these are the displacements that go in file07, but I'll use them from
        // memory to build the jobs
        let file07 = self.intder.convert_disps()?;

        let mut geoms = Vec::with_capacity(file07.len());
        for geom in file07 {
            // this is a bit unsightly, but I also don't want to duplicate the
            // `from_slices` code in psqs
            let mut mol = Molecule::from_slices(
                &atomic_numbers,
                &geom.as_slice()[..geom.len() - 3 * ndum],
            );
            mol.to_angstrom();
            geoms.push(Geom::from(mol));
        }
        Ok((geoms, taylor, atomic_numbers))
    }

    fn anpass<W: Write>(
        &self,
        dir: Option<impl AsRef<Path>>,
        energies: &mut [f64],
        taylor: &Taylor,
        step_size: f64,
        w: &mut W,
    ) -> Result<
        (Vec<rust_anpass::fc::Fc>, rust_anpass::Bias),
        Box<Result<(Spectro, Output), FreqError>>,
    > {
        make_rel(dir.as_ref(), energies);
        let anpass =
            Taylor::to_anpass(taylor, &taylor.disps(), energies, step_size);
        if let Some(dir) = dir.as_ref() {
            write_file(dir.as_ref().join("anpass.in"), &anpass).unwrap();
        }
        let (fcs, long_line, res) = if DEBUG {
            writeln!(w, "Anpass Input:\n{anpass}").unwrap();
            let (fcs, long_line, res) = match anpass.run_debug(w) {
                Ok(v) => v,
                Err(e) => return Err(Box::new(Err(FreqError(e.0)))),
            };
            writeln!(w, "\nStationary Point:\n{long_line}").unwrap();
            (fcs, long_line, res)
        } else {
            match anpass.run() {
                Ok(v) => (v.0, v.1, v.2),
                Err(e) => return Err(Box::new(Err(FreqError(e.0)))),
            }
        };
        writeln!(w, "anpass sum of squared residuals: {res:17.8e}").unwrap();
        Ok((fcs, long_line))
    }
}

/// find the minimum energy in `energies` and make the others relative to it
/// (subtract it from the others). also create `energy.dat` and `rel.dat` in the
/// current directory
pub(crate) fn make_rel(dir: Option<impl AsRef<Path>>, energies: &mut [f64]) {
    let min = energies
        .iter()
        .min_by(|a, b| a.total_cmp(b))
        .cloned()
        .unwrap();
    if let Some(dir) = dir {
        let dir = dir.as_ref();
        let mut efile = std::fs::File::create(dir.join("energy.dat")).unwrap();
        let mut rel = std::fs::File::create(dir.join("rel.dat")).unwrap();
        for energy in energies.iter() {
            writeln!(efile, "{energy:20.12}").unwrap();
            writeln!(rel, "{:20.12}", *energy - min).unwrap();
        }
    }
    for energy in energies.iter_mut() {
        *energy -= min;
    }
}

pub(crate) fn write_file(
    f: impl AsRef<Path>,
    d: impl Display,
) -> std::io::Result<()> {
    let mut f = std::fs::File::create(f)?;
    writeln!(f, "{d}")
}

/// an Error type containing information about a failure to run `freqs`
#[derive(Debug)]
pub struct FreqError(pub String);

impl Sic {
    /// prepare to [Self::generate_pts] by determining the symmetries of each
    /// SIC and sorting them into the correct order for [Taylor::make_checks]
    fn prepare_points<W: Write>(
        &mut self,
        mol: &Molecule,
        step_size: f64,
        pg: &PointGroup,
        w: &mut W,
    ) -> Result<Prep, IntderError> {
        let intder = &mut self.intder;
        let atomic_numbers = mol.atomic_numbers();
        let nsic = intder.symmetry_internals.len();
        let mut disps = Vec::new();
        for i in 0..nsic {
            let mut disp = vec![0.0; nsic];
            disp[i] = step_size;
            disps.push(disp);
        }
        intder.geom = intder::geom::Geom::from(mol.clone());
        intder.geom.to_bohr();
        intder.disps = disps;
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
        let disps = intder.convert_disps()?;
        let mut irreps = Vec::new();
        for (i, disp) in disps.iter().enumerate() {
            let disp = disp.as_slice();
            let m = Molecule::from_slices(
                &atomic_numbers,
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
        irreps.sort_by_key(|k| k.1);
        let just_irreps: Vec<_> = irreps.iter().map(|s| s.1).collect();
        let mut new_sics = Vec::new();
        for irrep in &irreps {
            new_sics.push(intder.symmetry_internals[irrep.0].clone());
        }
        intder.symmetry_internals = new_sics;
        writeln!(w, "\nSymmetry Internal Coordinates:").unwrap();
        intder.print_sics(w, &just_irreps);
        Ok(Prep {
            atomic_numbers,
            nsic,
            ndum,
            irreps,
        })
    }
    /// run the frequency portion of a QFF in `dir`. The caller is responsible
    /// for ensuring this directory exists.
    #[allow(clippy::too_many_arguments)]
    pub fn freqs<W: Write>(
        &mut self,
        w: &mut W,
        dir: Option<impl AsRef<Path>>,
        energies: &mut [f64],
        taylor: &Taylor,
        atomic_numbers: &AtomicNumbers,
        step_size: f64,
    ) -> Result<(Spectro, Output), FreqError> {
        let (fcs, long_line) =
            match self.anpass(dir.as_ref(), energies, taylor, step_size, w) {
                Ok(value) => value,
                Err(value) => return *value,
            };

        // intder_geom
        self.intder.disps = vec![long_line.disp.as_slice().to_vec()];
        if let Some(dir) = &dir {
            write_file(dir.as_ref().join("intder_geom.in"), &self.intder)
                .unwrap();
        }
        let refit_geom = self.intder.convert_disps().unwrap();
        let refit_geom = refit_geom[0].as_slice();
        let l = refit_geom.len() - 3 * self.intder.ndum();
        let dummies = &refit_geom[l..];
        let mol = Molecule::from_slices(atomic_numbers, &refit_geom[..l]);

        writeln!(w, "\nRefit Geometry\n{mol:20.12}").unwrap();

        self.intder.geom = intder::geom::Geom::from(mol.clone());
        for dummy in dummies.chunks(3) {
            self.intder.geom.push(vector![dummy[0], dummy[1], dummy[2]]);
        }

        self.intder.disps = vec![];

        // self.intder freqs
        for fc in fcs {
            // skip zeroth and first derivatives
            if (fc.1, fc.2, fc.3) != (0, 0, 0) {
                self.intder.add_fc(vec![fc.0, fc.1, fc.2, fc.3], fc.4);
            }
        }

        if DEBUG {
            writeln!(w, "Self.Intder Freqs Input\n{}", self.intder).unwrap();
        }

        self.intder.atoms = mol.atoms.iter().map(intder::Atom::from).collect();
        self.intder.input_options[3] = 4;
        self.intder.input_options[6] = 2;
        self.intder.input_options[10] = 3;
        self.intder.input_options[13] = 0;
        self.intder.input_options[14] = 0;

        if let Some(dir) = &dir {
            write_file(dir.as_ref().join("self.intder.in"), &self.intder)
                .unwrap();
        }

        let (f2, f3, f4) = self.intder.convert_fcs();
        if let Some(dir) = &dir {
            Intder::dump_fcs(dir.as_ref().to_str().unwrap(), &f2, &f3, &f4);
        }

        // spectro
        let mut spectro = Spectro::from(mol);
        spectro.header = SPECTRO_HEADER.to_vec();

        let fc3 = spectro::new_fc3(spectro.n3n, &f3);
        let fc4 = spectro::new_fc4(spectro.n3n, &f4);

        if DEBUG {
            writeln!(w, "Spectro Input:\n{spectro}").unwrap();
        }

        if let Some(dir) = dir {
            let input = dir.as_ref().join("spectro.in");
            spectro.write(input).unwrap();
        }

        let (output, _) =
            spectro.run(spectro::Derivative::Quartic(f2, fc3, fc4));

        Ok((spectro, output))
    }
}
