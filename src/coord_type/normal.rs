//! automatic normal coordinates. the normal coordinates are obtained by
//! running a harmonic, cartesian force field and feeding that to spectro to
//! obtain the LXM matrix. then the cartesian QFF machinery is used to
//! generate the normal coordinate displacements, which can be fed back in
//! to spectro at the end as f3qcm and f4qcm

use std::{io::Write, marker::Sync};

use intder::{fc3_index, fc4_index};
use nalgebra::DVector;
use psqs::{
    geom::Geom,
    program::{Job, Program, Template},
    queue::Queue,
};
use rust_anpass::fc::Fc;
use spectro::{Output, Spectro};
use symm::{Irrep, Molecule, PointGroup};
use taylor::{Disps, Taylor};

use crate::{cleanup, config::Config, coord_type::write_file};

use super::{
    findiff::{atom_parts, bighash::BigHash, zip_atoms, FiniteDifference},
    fitting::{AtomicNumbers, Fitted},
    sic::DEBUG,
    Cart, CoordType, Derivative, FreqError, Nderiv, SPECTRO_HEADER,
};

#[derive(Default)]
pub struct Normal {
    /// the normal coordinates, called the LXM matrix in spectro
    lxm: Option<nalgebra::DMatrix<f64>>,

    /// 1/√m where m is the atomic mass
    m12: Vec<f64>,

    /// the number of normal coordinates. convenient to have here so I don't
    /// have to keep thinking about 3n-6/5 stuff
    ncoords: usize,

    /// whether to use finite differences or least-squares fitting for computing
    /// the derivatives. defaults to false => use the fitting
    findiff: bool,

    irreps: Option<Vec<Irrep>>,
}

impl Normal {
    pub fn findiff(findiff: bool) -> Self {
        Self {
            findiff,
            ..Default::default()
        }
    }
}

impl<W, Q, P> CoordType<W, Q, P> for Normal
where
    W: Write,
    Q: Queue<P> + Sync,
    P: Program + Clone + Send + Sync,
{
    fn run(
        mut self,
        w: &mut W,
        queue: &Q,
        config: &Config,
    ) -> (Spectro, Output) {
        let (mut s, o, ref_energy, pg) = self.cart_part(config, queue, w);
        cleanup();
        let _ = std::fs::create_dir("pts");
        // pretty sure I assert lxm is square somewhere in spectro though. lxm
        // should be in column-major order so I think this is all right
        let cols = o.lxm.len();
        let rows = o.lxm[0].len();

        // actual initialization of self, must happen before build_points
        self.lxm = Some(nalgebra::DMatrix::from_iterator(
            rows,
            cols,
            o.lxm.iter().flatten().cloned(),
        ));
        self.m12 = o.geom.weights().iter().map(|w| 1.0 / w.sqrt()).collect();
        // 3n - 6 + 1 if linear = 3n - 5
        self.ncoords =
            3 * o.geom.atoms.len() - 6 + s.rotor.is_linear() as usize;
        let template = Template::from(&config.template);

        let (f3qcm, f4qcm) = if self.findiff {
            let n = self.ncoords;
            let nfc2 = n * n;
            let nfc3 = n * (n + 1) * (n + 2) / 6;
            let nfc4 = n * (n + 1) * (n + 2) * (n + 3) / 24;
            let mut fcs = vec![0.0; nfc2 + nfc3 + nfc4];
            let mut map = BigHash::new(o.geom.clone(), pg);
            let geoms = self.build_points(
                Geom::Xyz(o.geom.atoms.clone()),
                config.step_size,
                ref_energy,
                Derivative::Quartic(nfc2, nfc3, nfc4),
                &mut fcs,
                &mut map,
                n,
            );
            let dir = "pts";
            let jobs: Vec<_> = geoms
                .into_iter()
                .enumerate()
                .map(|(job_num, mol)| {
                    let filename = format!("{dir}/job.{:08}", job_num);
                    Job::new(
                        P::new(
                            filename,
                            template.clone(),
                            config.charge,
                            mol.geom,
                        ),
                        mol.index,
                    )
                })
                .collect();
            writeln!(
                w,
                "{n} normal coordinates requires {} points",
                jobs.len()
            )
            .unwrap();
            time!(w, "draining points",
                  // drain into energies
                  let mut energies = vec![0.0; jobs.len()];
                  let time = queue
                  .drain(dir, jobs, &mut energies)
                  .expect("single-point calculations failed");
            );
            eprintln!("total job time: {time:.1} sec");

            self.map_energies(&map, &energies, &mut fcs);

            let cubs = &fcs[nfc2..nfc2 + nfc3];
            let quarts = &fcs[nfc2 + nfc3..];

            to_qcm(&o.harms, n, cubs, quarts, intder::HART)
        } else {
            self.irreps = Some(o.irreps.clone());
            let (geoms, taylor, taylor_disps, _atomic_numbers) = self
                .generate_pts(w, &o.geom, &pg, config.step_size)
                .unwrap();
            let dir = "pts";
            let jobs =
                P::build_jobs(&geoms, dir, 0, 1.0, 0, config.charge, template);

            writeln!(
                w,
                "\n{} atoms require {} jobs",
                o.geom.atoms.len(),
                jobs.len()
            )
            .unwrap();

            let mut energies = vec![0.0; jobs.len()];
            let time = queue
                .drain(dir, jobs, &mut energies)
                .expect("single-point energies failed");
            eprintln!("total job time: {time:.1} sec");

            let _ = std::fs::create_dir("freqs");
            let (fcs, long_line) = self
                .anpass(
                    "freqs",
                    &mut energies,
                    &taylor,
                    &taylor_disps,
                    config.step_size,
                    w,
                )
                .unwrap();

            println!("before refit: s.geom={:.8}", s.geom);
            let ll = long_line.disp;
            let lxm = self.lxm.as_ref().unwrap();
            let n3n = o.geom.atoms.len() * 3;
            let (names, mut coords) = atom_parts(&o.geom.atoms);
            for k in 0..n3n {
                for n in 0..self.ncoords {
                    coords[k] += self.m12[k / 3] * lxm[(k, n)] * ll[n];
                }
            }
            s.geom = Molecule::new(zip_atoms(&names, coords.into()));
            println!("after refit: s.geom={:.8}", s.geom);

            let mut f3qcm = Vec::new();
            let mut f4qcm = Vec::new();
            for Fc(i, j, k, l, val) in fcs {
                // adapted from Intder::add_fc
                match (k, l) {
                    (0, 0) => {
                        // harmonic, skip it
                    }
                    (_, 0) => {
                        let idx = fc3_index(i, j, k);
                        if f3qcm.len() <= idx {
                            f3qcm.resize(idx + 1, 0.0);
                        }
                        f3qcm[idx] = val;
                    }
                    (_, _) => {
                        let idx = fc4_index(i, j, k, l);
                        if f4qcm.len() <= idx {
                            f4qcm.resize(idx + 1, 0.0);
                        }
                        f4qcm[idx] = val;
                    }
                }
            }

            to_qcm(&o.harms, self.ncoords, &f3qcm, &f4qcm, 1.0)
        };
        let (o, _) = s.finish(
            DVector::from(o.harms.clone()),
            spectro::F3qcm::new(f3qcm),
            spectro::F4qcm::new(f4qcm),
            o.irreps,
            self.lxm.unwrap(),
        );

        (s, o)
    }
}

/// convert the force constants in `cubs` and `quarts` to wavenumbers, as
/// expected by spectro. `fac` should be [intder::HART] for finite difference
/// fcs and 1.0 for fitted ones
fn to_qcm(
    harms: &[f64],
    n: usize,
    cubs: &[f64],
    quarts: &[f64],
    fac: f64,
) -> (Vec<f64>, Vec<f64>) {
    let mut ijk = 0;
    let mut ijkl = 0;
    let mut f3qcm = Vec::new();
    let mut f4qcm = Vec::new();
    for i in 0..n {
        let wi = harms[i];
        for j in 0..=i {
            let wj = harms[j];
            for k in 0..=j {
                let wk = harms[k];
                let wijk = wi * wj * wk;
                let fact = fac * spectro::consts::FACT3 / wijk.sqrt();
                let val = cubs[ijk];
                f3qcm.push(val * fact);
                ijk += 1;
                (0..=k).for_each(|l| {
                    let wl = harms[l];
                    let wijkl = wijk * wl;
                    let fact = fac * spectro::consts::FACT4 / wijkl.sqrt();
                    let val = quarts[ijkl];
                    f4qcm.push(val * fact);
                    ijkl += 1;
                });
            }
        }
    }
    (f3qcm, f4qcm)
}

impl Fitted for Normal {
    type Prep = ();

    type Error = ();

    fn prepare_points<W: Write>(
        &mut self,
        _mol: &Molecule,
        _step_size: f64,
        _pg: &PointGroup,
        _w: &mut W,
    ) -> Result<Self::Prep, Self::Error> {
        Ok(())
    }

    fn generate_pts<W: Write>(
        &mut self,
        _w: &mut W,
        mol: &Molecule,
        pg: &PointGroup,
        step_size: f64,
    ) -> Result<(Vec<Geom>, Taylor, Disps, AtomicNumbers), Self::Error> {
        let mut irreps: Vec<_> = self
            .irreps
            .as_ref()
            .unwrap()
            .iter()
            .cloned()
            .enumerate()
            .collect();
        irreps.sort_by_key(|(_, irrep)| *irrep);
        let checks = Taylor::make_checks(irreps, pg);
        let taylor = Taylor::new(5, self.ncoords, checks.0, checks.1);
        let taylor_disps = taylor.disps();
        let disps = taylor_disps.to_intder(step_size);
        let (names, coords) = atom_parts(&mol.atoms);
        let lxm = self.lxm.as_ref().unwrap();
        let mut geoms = Vec::new();
        for dq in disps {
            let mut coords = coords.clone();
            let nc = coords.len();
            for k in 0..nc {
                for n in 0..self.ncoords {
                    coords[k] += self.m12[k / 3] * lxm[(k, n)] * dq[n];
                }
            }
            geoms.push(Geom::Xyz(zip_atoms(&names, coords.into())))
        }
        Ok((geoms, taylor, taylor_disps, mol.atomic_numbers()))
    }

    /// perform the initial anpass fitting, but we can't refit on the geometry
    /// because it would change the geometry, which would change the normal
    /// coordinates themselves and invalidate all of the force constants
    fn anpass<W: Write>(
        &self,
        dir: &str,
        energies: &mut [f64],
        taylor: &Taylor,
        taylor_disps: &taylor::Disps,
        step_size: f64,
        w: &mut W,
    ) -> Result<
        (Vec<rust_anpass::fc::Fc>, rust_anpass::Bias),
        Box<Result<(Spectro, Output), super::FreqError>>,
    > {
        let min = energies.iter().cloned().reduce(f64::min).unwrap();
        for energy in energies.iter_mut() {
            *energy -= min;
        }
        let anpass =
            Taylor::to_anpass(taylor, taylor_disps, energies, step_size);
        write_file(format!("{dir}/anpass.in"), &anpass).unwrap();
        let (fcs, long_line, res) = if DEBUG {
            writeln!(w, "Anpass Input:\n{}", anpass).unwrap();
            let (fcs, long_line, res) = match anpass.run_debug(w) {
                Ok(v) => v,
                Err(e) => return Err(Box::new(Err(FreqError(e.0)))),
            };
            writeln!(w, "\nStationary Point:\n{}", long_line).unwrap();
            (fcs, long_line, res)
        } else {
            match anpass.run() {
                Ok(v) => v,
                Err(e) => return Err(Box::new(Err(FreqError(e.0)))),
            }
        };
        writeln!(w, "anpass sum of squared residuals: {:17.8e}", res).unwrap();
        Ok((fcs, long_line))
    }
}

impl FiniteDifference for Normal {
    fn scale(&self, nderiv: usize, step_size: f64) -> f64 {
        match nderiv {
            2 => 1.0 / (4.0 * step_size * step_size),
            3 => 1.0 / (8.0 * step_size * step_size * step_size),
            4 => 1.0 / (16.0 * step_size * step_size * step_size * step_size),
            _ => panic!("unrecognized derivative level"),
        }
    }

    fn new_geom(
        &self,
        names: &[&str],
        coords: nalgebra::DVector<f64>,
        step_size: f64,
        steps: Vec<isize>,
    ) -> psqs::geom::Geom {
        // python version:
        // N = 3
        // v = [0.0] * (3 * N)
        // for k in range(3 * N):
        //     for n in range(3 * N - 6):
        //         v[k] += m12[k] * lxm[k, n] * dq[n]
        // the normal coordinate displacement
        let mut dq = vec![0.0; self.ncoords];
        for step in steps {
            if step < 1 {
                dq[(-step - 1) as usize] -= step_size;
            } else {
                dq[(step - 1) as usize] += step_size;
            }
        }

        let lxm = self.lxm.as_ref().unwrap();

        let nc = coords.len();
        let mut v = DVector::zeros(nc);
        // TODO do this as a mat mul, but get it to work first
        for k in 0..nc {
            for n in 0..self.ncoords {
                v[k] += self.m12[k / 3] * lxm[(k, n)] * dq[n];
            }
        }
        let coords = coords + v;
        Geom::Xyz(zip_atoms(names, coords))
    }
}

impl Normal {
    /// run the Cartesian harmonic force field and return the spectro output,
    /// from which we can extract the geometry and normal coordinates (lxm)
    fn cart_part<P, Q, W>(
        &self,
        config: &Config,
        queue: &Q,
        w: &mut W,
    ) -> (Spectro, Output, f64, PointGroup)
    where
        P: Program + Clone + Send + Sync,
        Q: Queue<P> + Sync,
        W: Write,
    {
        let (
            n,
            nfc2,
            nfc3,
            mut fcs,
            mol,
            energies,
            mut target_map,
            ref_energy,
            pg,
        ) = Cart.first_part(w, config, queue, Nderiv::Two);
        let (fc2, _, _) = self.make_fcs(
            &mut target_map,
            &energies,
            &mut fcs,
            n,
            Derivative::Quartic(nfc2, nfc3, 0),
            "freqs",
        );
        let (a, b) = self.harm_freqs("freqs", &mol, fc2);
        (a, b, ref_energy, pg)
    }

    /// run the harmonic frequencies through spectro
    fn harm_freqs(
        &self,
        dir: &str,
        mol: &Molecule,
        fc2: nalgebra::DMatrix<f64>,
    ) -> (Spectro, Output) {
        let mut mol = mol.clone();
        mol.to_bohr();
        let mut spectro = Spectro::from(mol);
        spectro.header = SPECTRO_HEADER.to_vec();

        // write input
        let input = format!("{}/spectro.in", dir);
        spectro.write(&input).unwrap();

        let (output, _) = spectro.run(spectro::Derivative::Harmonic(fc2));
        (spectro, output)
    }
}
