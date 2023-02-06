//! automatic normal coordinates. the normal coordinates are obtained by
//! running a harmonic, cartesian force field and feeding that to spectro to
//! obtain the LXM matrix. then the cartesian QFF machinery is used to
//! generate the normal coordinate displacements, which can be fed back in
//! to spectro at the end as f3qcm and f4qcm

use std::{io::Write, marker::Sync};

use intder::Intder;
pub use intder::{fc3_index, fc4_index};
use nalgebra::DVector;
use psqs::{
    geom::Geom,
    program::{Job, Program, Template},
    queue::Queue,
};
use rust_anpass::Dmat;
pub use rust_anpass::{fc::Fc, Bias};
use serde::{Deserialize, Serialize};
pub use spectro::{F3qcm, F4qcm, Output, Spectro};
use symm::{Irrep, Molecule, PointGroup};
use taylor::{Disps, Taylor};

use crate::{
    cleanup,
    config::Config,
    coord_type::{write_file, CHK_NAME},
};

use super::{
    findiff::{atom_parts, bighash::BigHash, zip_atoms, FiniteDifference},
    fitting::{AtomicNumbers, Fitted},
    make_rel, Cart, CoordType, Derivative, FirstPart, Load, Nderiv,
    SPECTRO_HEADER,
};

#[derive(Clone, Default, Debug, Serialize, Deserialize, PartialEq)]
pub struct Normal {
    /// the normal coordinates, called the LXM matrix in spectro
    pub lxm: Option<nalgebra::DMatrix<f64>>,

    /// lx matrix, also passed to spectro::finish
    pub lx: Option<nalgebra::DMatrix<f64>>,

    /// 1/√m where m is the atomic mass
    pub m12: Vec<f64>,

    /// the number of normal coordinates. convenient to have here so I don't
    /// have to keep thinking about 3n-6/5 stuff
    pub ncoords: usize,

    /// whether to use finite differences or least-squares fitting for computing
    /// the derivatives. defaults to false => use the fitting
    pub findiff: bool,

    pub irreps: Option<Vec<Irrep>>,
}

impl Normal {
    pub fn findiff(findiff: bool) -> Self {
        Self {
            findiff,
            ..Default::default()
        }
    }

    /// actual initialization of self with the results of the initial harmonic
    /// FF, must happen before build_points. TODO use phantomdata/generic trick
    /// to enforce that
    pub fn prep_qff<W>(&mut self, w: &mut W, o: &Output)
    where
        W: Write,
    {
        // pretty sure I assert lxm is square somewhere in spectro though. lxm
        // should be in column-major order so I think this is all right
        let cols = o.lxm.len();
        let rows = o.lxm[0].len();
        let lxm = nalgebra::DMatrix::from_iterator(
            rows,
            cols,
            o.lxm.iter().flatten().cloned(),
        );
        let lx = nalgebra::DMatrix::from_iterator(
            rows,
            cols,
            o.lx.iter().flatten().cloned(),
        );
        writeln!(w, "Normal Coordinates:{lxm:.8}").unwrap();
        writeln!(w, "Harmonic Frequencies:").unwrap();
        for (i, (r, h)) in o.irreps.iter().zip(&o.harms).enumerate() {
            writeln!(w, "{i:5}{r:>5}{h:8.1}").unwrap();
        }
        self.lxm = Some(lxm);
        self.lx = Some(lx);
        self.m12 = o.geom.weights().iter().map(|w| 1.0 / w.sqrt()).collect();
        self.ncoords = o.harms.len();
        self.irreps = Some(o.irreps.clone());
    }

    /// run the QFF using a least-squares fitting
    #[allow(clippy::too_many_arguments)]
    fn run_fitted<P, W, Q>(
        &mut self,
        o: &Output,
        s: &Spectro,
        w: &mut W,
        pg: PointGroup,
        config: &Config,
        template: Template,
        queue: &Q,
    ) -> (Vec<f64>, Vec<f64>)
    where
        W: Write,
        Q: Queue<P> + Sync,
        P: Program + Clone + Send + Sync + Serialize + for<'a> Deserialize<'a>,
    {
        let (geoms, taylor, taylor_disps, _atomic_numbers) = self
            .generate_pts(w, &o.geom, &pg, config.step_size)
            .unwrap();
        let dir = "pts";
        let jobs =
            P::build_jobs(&geoms, dir, 0, 1.0, 0, config.charge, template);
        writeln!(
            w,
            "{} normal coordinates require {} points",
            self.ncoords,
            jobs.len()
        )
        .unwrap();

        let resume = Resume {
            normal: self.clone(),
            njobs: jobs.len(),
            deriv: DerivType::Fitted {
                taylor,
                taylor_disps,
                step_size: config.step_size,
            },
            output: o.clone(),
            spectro: s.clone(),
        };
        resume.dump(CHK_NAME);

        let DerivType::Fitted { taylor, taylor_disps, step_size, ..} =
	    resume.deriv else {
		unreachable!();
	    };

        let mut energies = vec![0.0; jobs.len()];
        let time = queue
            .drain(dir, jobs, &mut energies, config.check_int)
            .expect("single-point energies failed");
        eprintln!("total job time: {time:.1} sec");
        let (fcs, _) = self
            .anpass(
                Some("freqs"),
                &mut energies,
                &taylor,
                &taylor_disps,
                step_size,
                w,
            )
            .unwrap();
        // needed in case taylor eliminated some of the higher derivatives
        // by symmetry. this should give the maximum, full sizes without
        // resizing
        let n = self.ncoords;
        let mut f3qcm = vec![0.0; fc3_index(n, n, n) + 1];
        let mut f4qcm = vec![0.0; fc4_index(n, n, n, n) + 1];
        for Fc(i, j, k, l, val) in fcs {
            // adapted from Intder::add_fc
            match (k, l) {
                (0, 0) => {
                    // harmonic, skip it
                }
                (_, 0) => {
                    let idx = fc3_index(i, j, k);
                    f3qcm[idx] = val;
                }
                (_, _) => {
                    let idx = fc4_index(i, j, k, l);
                    f4qcm[idx] = val;
                }
            }
        }
        to_qcm(&o.harms, self.ncoords, &f3qcm, &f4qcm, 1.0)
    }

    /// run the QFF using finite differences
    #[allow(clippy::too_many_arguments)]
    fn run_findiff<P, W, Q>(
        &self,
        o: &Output,
        s: &Spectro,
        pg: PointGroup,
        config: &Config,
        ref_energy: f64,
        template: &Template,
        w: &mut W,
        queue: &Q,
    ) -> (Vec<f64>, Vec<f64>)
    where
        W: Write,
        Q: Queue<P> + Sync,
        P: Program + Clone + Send + Sync + Serialize + for<'a> Deserialize<'a>,
    {
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
                let filename = format!("{dir}/job.{job_num:08}");
                Job::new(
                    P::new(filename, template.clone(), config.charge, mol.geom),
                    mol.index,
                )
            })
            .collect();
        writeln!(w, "{n} normal coordinates require {} points", jobs.len())
            .unwrap();

        let resume = Resume {
            normal: self.clone(),
            njobs: jobs.len(),
            deriv: DerivType::Findiff {
                map,
                fcs,
                n,
                nfc2,
                nfc3,
            },
            output: o.clone(),
            spectro: s.clone(),
        };
        resume.dump(CHK_NAME);

        time!(w, "draining points",
              // drain into energies
              let mut energies = vec![0.0; jobs.len()];
              let time = queue
              .drain(dir, jobs, &mut energies, config.check_int)
              .expect("single-point calculations failed");
        );
        eprintln!("total job time: {time:.1} sec");
        let DerivType::Findiff{ map, mut fcs, .. } = resume.deriv  else  {
	    unreachable!()
	};
        self.map_energies(&map, &energies, &mut fcs);
        let cubs = &fcs[nfc2..nfc2 + nfc3];
        let quarts = &fcs[nfc2 + nfc3..];
        to_qcm(&o.harms, n, cubs, quarts, intder::HART)
    }
}

impl<W, Q, P> CoordType<W, Q, P> for Normal
where
    W: Write,
    Q: Queue<P> + Sync,
    P: Program + Clone + Send + Sync + Serialize + for<'a> Deserialize<'a>,
{
    fn run(
        mut self,
        w: &mut W,
        queue: &Q,
        config: &Config,
    ) -> (Spectro, Output) {
        let (s, o, ref_energy, pg, fc2) =
            self.cart_part(&FirstPart::from(config.clone()), queue, w, "pts");
        cleanup();
        let _ = std::fs::create_dir("pts");
        self.prep_qff(w, &o);

        let tmpl = config.template.clone().into();

        // TODO have to split these run_* methods into a prep_* and run_* so I
        // can save the Resume between them
        let _ = std::fs::create_dir("freqs");
        let (f3qcm, f4qcm) = if self.findiff {
            self.run_findiff(&o, &s, pg, config, ref_energy, &tmpl, w, queue)
        } else {
            self.run_fitted(&o, &s, w, pg, config, tmpl, queue)
        };
        Intder::dump_fcs("freqs", &fc2, &f3qcm, &f4qcm);
        s.write("freqs/spectro.in").unwrap();
        let fin = spectro::SpectroFinish::new(
            DVector::from(o.harms.clone()),
            F3qcm::new(f3qcm),
            F4qcm::new(f4qcm),
            o.irreps,
            self.lxm.unwrap(),
            self.lx.unwrap(),
        );
        fin.dump("freqs/finish.spectro").unwrap_or_else(|e| {
            eprintln!("failed to dump finish.spectro with `{e}`")
        });
        let (o, _) = s.finish(fin);

        (s, o)
    }

    type Resume = Resume;

    fn resume(
        #[allow(unused_assignments)] mut self,
        w: &mut W,
        queue: &Q,
        Resume {
            normal,
            njobs,
            deriv,
            output,
            spectro,
        }: Resume,
    ) -> (Spectro, Output) {
        self = normal;
        let dir = "pts";
        time!(w, "draining points",
              // drain into energies
              let mut energies = vec![0.0; njobs];
              let time = queue
              .resume(dir, "chk.json", &mut energies, 0)
              .expect("single-point calculations failed");
        );
        eprintln!("total job time: {time:.1} sec");
        let (f3qcm, f4qcm) = match deriv {
            DerivType::Findiff {
                map,
                mut fcs,
                n,
                nfc2,
                nfc3,
            } => {
                self.map_energies(&map, &energies, &mut fcs);
                let cubs = &fcs[nfc2..nfc2 + nfc3];
                let quarts = &fcs[nfc2 + nfc3..];
                to_qcm(&output.harms, n, cubs, quarts, intder::HART)
            }
            DerivType::Fitted {
                taylor,
                taylor_disps,
                step_size,
            } => {
                let _ = std::fs::create_dir("freqs");
                let (fcs, _) = self
                    .anpass(
                        Some("freqs"),
                        &mut energies,
                        &taylor,
                        &taylor_disps,
                        step_size,
                        w,
                    )
                    .unwrap();
                // needed in case taylor eliminated some of the higher
                // derivatives by symmetry. this should give the maximum, full
                // sizes without resizing
                let n = self.ncoords;
                let mut f3qcm = vec![0.0; fc3_index(n, n, n) + 1];
                let mut f4qcm = vec![0.0; fc4_index(n, n, n, n) + 1];
                for Fc(i, j, k, l, val) in fcs {
                    // adapted from Intder::add_fc
                    match (k, l) {
                        (0, 0) => {
                            // harmonic, skip it
                        }
                        (_, 0) => {
                            let idx = fc3_index(i, j, k);
                            f3qcm[idx] = val;
                        }
                        (_, _) => {
                            let idx = fc4_index(i, j, k, l);
                            f4qcm[idx] = val;
                        }
                    }
                }
                to_qcm(&output.harms, self.ncoords, &f3qcm, &f4qcm, 1.0)
            }
        };

        let fin = spectro::SpectroFinish::new(
            DVector::from(output.harms.clone()),
            F3qcm::new(f3qcm),
            F4qcm::new(f4qcm),
            output.irreps,
            self.lxm.unwrap(),
            self.lx.unwrap(),
        );

        fin.dump("freqs/finish.spectro").unwrap_or_else(|e| {
            eprintln!("failed to dump finish.spectro with `{e}`")
        });
        let (o, _) = spectro.finish(fin);

        (spectro, o)
    }
}

#[derive(Serialize, Deserialize, PartialEq, Debug)]
pub enum DerivType {
    Findiff {
        map: BigHash,
        fcs: Vec<f64>,
        n: usize,
        nfc2: usize,
        nfc3: usize,
    },
    Fitted {
        taylor: Taylor,
        taylor_disps: Disps,
        step_size: f64,
    },
}

#[derive(Serialize, Deserialize, PartialEq, Debug)]
pub struct Resume {
    pub(crate) normal: Normal,
    pub(crate) njobs: usize,
    pub(crate) output: Output,
    pub(crate) spectro: Spectro,
    pub(crate) deriv: DerivType,
}

impl Resume {
    pub fn new(
        normal: Normal,
        njobs: usize,
        output: Output,
        spectro: Spectro,
        deriv: DerivType,
    ) -> Self {
        Self {
            normal,
            njobs,
            output,
            spectro,
            deriv,
        }
    }
}

impl Load for Resume {}

/// convert the force constants in `cubs` and `quarts` to wavenumbers, as
/// expected by spectro. `fac` should be [intder::HART] for finite difference
/// fcs and 1.0 for fitted ones
pub fn to_qcm(
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
        dir: Option<&str>,
        energies: &mut [f64],
        taylor: &Taylor,
        taylor_disps: &taylor::Disps,
        step_size: f64,
        w: &mut W,
    ) -> Result<
        (Vec<rust_anpass::fc::Fc>, rust_anpass::Bias),
        Box<Result<(Spectro, Output), super::FreqError>>,
    > {
        make_rel(energies);
        let anpass =
            Taylor::to_anpass(taylor, taylor_disps, energies, step_size);
        if let Some(dir) = dir {
            write_file(format!("{dir}/anpass.in"), &anpass).unwrap();
        }
        let (fcs, f) = anpass.fit();
        writeln!(
            w,
            "anpass sum of squared residuals: {:17.8e}",
            anpass.residuals(&fcs, &f)
        )
        .unwrap();
        Ok((anpass.make9903(&fcs), Bias::default()))
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
    pub fn cart_part<P, Q, W>(
        &self,
        config: &FirstPart,
        queue: &Q,
        w: &mut W,
        dir: &str,
    ) -> (Spectro, Output, f64, PointGroup, Dmat)
    where
        P: Program + Clone + Send + Sync + Serialize + for<'a> Deserialize<'a>,
        Q: Queue<P> + Sync,
        W: Write,
    {
        let (
            n,
            nfc2,
            _,
            mut fcs,
            mol,
            energies,
            mut target_map,
            ref_energy,
            pg,
        ) = Cart.first_part(w, config, queue, Nderiv::Two, dir);
        let (fc2, _, _) = self.make_fcs(
            &mut target_map,
            &energies,
            &mut fcs,
            n,
            Derivative::Harmonic(nfc2),
            "freqs",
        );
        let (spectro, output) = self.harm_freqs("freqs", &mol, fc2.clone());
        (spectro, output, ref_energy, pg, fc2)
    }

    /// run the harmonic frequencies through spectro
    pub fn harm_freqs(
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
        let input = format!("{dir}/spectro.in");
        spectro.write(&input).unwrap();

        let (output, _) = spectro.run(spectro::Derivative::Harmonic(fc2));
        (spectro, output)
    }
}
