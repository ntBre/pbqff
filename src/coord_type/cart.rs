//! Cartesian coordinate QFFs.

use std::{error::Error, io, marker::Sync, path::Path};

use psqs::{
    geom::Geom,
    program::{Job, Program, Template},
    queue::Queue,
};
use serde::{Deserialize, Serialize};
use spectro::{Output, Spectro};
use symm::{Molecule, PointGroup};

use crate::{config::Config, optimize, ref_energy};

use super::{
    findiff::bighash::BigHash,
    findiff::{bighash::Target, zip_atoms, FiniteDifference},
    CoordType, Load, SPECTRO_HEADER,
};

/// debugging options. currently supported options: disp, fcs, none
pub(crate) static DEBUG: &str = "none";

pub struct Cart;

#[derive(Debug)]
pub struct CartGeom {
    pub geom: Geom,
    pub coeff: f64,
    pub index: usize,
}

pub enum Derivative {
    Harmonic(usize),
    Cubic(usize, usize),
    Quartic(usize, usize, usize),
}

pub enum Nderiv {
    Two,
    Four,
}

pub fn freqs(
    dir: Option<impl AsRef<Path>>,
    mol: &Molecule,
    fc2: nalgebra::DMatrix<f64>,
    f3: &[f64],
    f4: &[f64],
) -> (Spectro, Output) {
    let mut mol = mol.clone();
    mol.to_bohr();
    let mut spectro = Spectro::from(mol);
    spectro.header = SPECTRO_HEADER.to_vec();

    // write input
    if let Some(dir) = dir {
        let input = dir.as_ref().join("spectro.in");
        spectro.write(input).unwrap();
    }

    let fc3 = spectro::new_fc3(spectro.n3n, f3);
    let fc4 = spectro::new_fc4(spectro.n3n, f4);

    let output = spectro.run(spectro::Derivative::Quartic(fc2, fc3, fc4));
    (spectro, output)
}

impl<W, Q, P> CoordType<W, Q, P> for Cart
where
    W: io::Write,
    Q: Queue<P> + Sync,
    P: Program + Clone + Send + Sync + Serialize + for<'a> Deserialize<'a>,
{
    fn run(
        self,
        dir: impl AsRef<std::path::Path>,
        w: &mut W,
        queue: &Q,
        config: &Config,
    ) -> (Spectro, Output) {
        let FirstOutput {
            n,
            nfc2,
            nfc3,
            mut fcs,
            mol,
            energies,
            targets,
            ..
        } = Cart
            .first_part(
                w,
                &FirstPart::from(config.clone()),
                queue,
                Nderiv::Four,
                dir.as_ref().join("pts"),
            )
            .unwrap();

        let freq_dir = &dir.as_ref().join("freqs");
        time!(w, "freqs",
          let (fc2, f3, f4) = self.make_fcs(
          targets,
          &energies,
          &mut fcs,
          n,
          Derivative::Quartic(
          nfc2,
          nfc3,
          0),
          Some(freq_dir),
          );

          let r = freqs(Some(freq_dir), &mol, fc2, f3, f4);
        );
        r
    }

    type Resume = Resume;

    fn resume(
        self,
        _dir: impl AsRef<std::path::Path>,
        _w: &mut W,
        _queue: &Q,
        _config: &Config,
        _resume: Resume,
    ) -> (Spectro, Output) {
        todo!()
    }
}

#[derive(Serialize, Deserialize)]
pub struct Resume;

impl Load for Resume {}

impl FiniteDifference for Cart {
    fn new_geom(
        &self,
        names: &[&str],
        mut coords: Vec<f64>,
        step_size: f64,
        steps: Vec<isize>,
    ) -> Geom {
        for step in steps {
            if step < 1 {
                coords[(-step - 1) as usize] -= step_size;
            } else {
                coords[(step - 1) as usize] += step_size;
            }
        }
        Geom::Xyz(zip_atoms(names, coords))
    }
}

/// contains just the fields of Config needed for running [Cart::first_part]
pub struct FirstPart {
    pub template: String,
    pub optimize: bool,
    pub geometry: Geom,
    pub charge: isize,
    pub step_size: f64,
    pub weights: Option<Vec<f64>>,
}

impl From<Config> for FirstPart {
    fn from(config: Config) -> Self {
        Self {
            template: config.template,
            optimize: config.optimize,
            geometry: config.geometry,
            charge: config.charge,
            step_size: config.step_size,
            weights: config.weights,
        }
    }
}

impl Cart {
    /// run the "first part" of the Cartesian QFF, including the optimization if
    /// requested and the generation and running of the single-point energies
    pub fn first_part<W, Q, P>(
        &self,
        w: &mut W,
        config: &FirstPart,
        queue: &Q,
        nderiv: Nderiv,
        dir: impl AsRef<Path>,
    ) -> Result<FirstOutput, Box<dyn Error>>
    where
        W: io::Write,
        Q: Queue<P> + Sync,
        P: Program + Clone + Send + Sync + Serialize + for<'a> Deserialize<'a>,
    {
        time!(w, "opt",
            let template = Template::from(&config.template);
            let (geom, ref_energy) = if config.optimize {
                let res = optimize(
            &dir,
                    queue,
                    config.geometry.clone(),
                    template.clone(),
                    config.charge,
                )
                .expect("optimization failed");
                let Some(cart) = res.cart_geom else {
                    panic!("failed to extract cart geom from {res:?}");
                };
                (Geom::Xyz(cart), res.energy)
            } else {
                let ref_energy = ref_energy(
                    queue,
                    config.geometry.clone(),
                    template.clone(),
                    config.charge,
                );
                (config.geometry.clone(), ref_energy)
            };
        );
        let geom = geom.xyz().expect("expected an XYZ geometry, not Zmat");
        // 3 * #atoms
        let n = 3 * geom.len();
        let nfc2 = n * n;
        let nfc3 = n * (n + 1) * (n + 2) / 6;
        let nfc4 = n * (n + 1) * (n + 2) * (n + 3) / 24;
        let deriv = match nderiv {
            Nderiv::Two => Derivative::Harmonic(nfc2),
            Nderiv::Four => Derivative::Quartic(nfc2, nfc3, nfc4),
        };
        let mut fcs = vec![0.0; nfc2 + nfc3 + nfc4];
        let mut mol = Molecule::new(geom.to_vec());
        if let Some(ws) = &config.weights {
            for (i, w) in ws.iter().enumerate() {
                mol.atoms[i].weight = Some(*w);
            }
        }
        mol.normalize();
        let pg = mol.point_group();
        writeln!(w, "normalized geometry:\n{mol}").unwrap();
        writeln!(w, "point group:{pg}").unwrap();
        let mut target_map = BigHash::new(mol.clone(), pg);
        time! (w, "building points",
               let geoms = self.build_points(
               Geom::Xyz(mol.atoms.clone()),
               config.step_size,
           ref_energy,
               deriv,
               &mut fcs,
               &mut target_map,
           n,
               );
        );
        let targets = target_map.values();
        let jobs: Vec<_> = geoms
            .into_iter()
            .enumerate()
            .map(|(job_num, mol)| {
                let filename = format!("job.{job_num:08}");
                let filename =
                    dir.as_ref().join(filename).to_string_lossy().to_string();
                Job::new(
                    P::new(filename, template.clone(), config.charge, mol.geom),
                    mol.index,
                )
            })
            .collect();
        writeln!(
            w,
            "{n} Cartesian coordinates requires {} points",
            jobs.len()
        )
        .unwrap();
        time!(w, "draining points",
              // drain into energies
              let mut energies = vec![0.0; jobs.len()];
              let time = queue
              .drain(dir.as_ref().to_str().unwrap(), jobs, &mut energies, psqs::queue::Check::None)?;
        );

        eprintln!("total job time: {time:.1} sec");

        Ok(FirstOutput {
            n,
            nfc2,
            nfc3,
            fcs,
            mol,
            energies,
            targets,
            ref_energy,
            pg,
        })
    }
}

pub struct FirstOutput {
    pub n: usize,
    pub nfc2: usize,
    pub nfc3: usize,
    pub fcs: Vec<f64>,
    pub mol: Molecule,
    pub energies: Vec<f64>,
    pub targets: Vec<Target>,
    pub ref_energy: f64,
    pub pg: PointGroup,
}
