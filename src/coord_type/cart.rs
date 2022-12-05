use std::{io, marker::Sync};

use psqs::{
    geom::Geom,
    program::{Job, Program, Template},
    queue::Queue,
};
use spectro::{Output, Spectro};
use symm::Molecule;

use crate::{config::Config, optimize, ref_energy};

use super::{
    findiff::bighash::BigHash, findiff::FiniteDifference, CoordType,
    SPECTRO_HEADER,
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

impl Cart {}

pub fn freqs(
    dir: &str,
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
    let input = format!("{}/spectro.in", dir);
    spectro.write(&input).unwrap();

    let fc3 = spectro::new_fc3(spectro.n3n, f3);
    let fc4 = spectro::new_fc4(spectro.n3n, f4);

    let (output, _) = spectro.run(spectro::Derivative::Quartic(fc2, fc3, fc4));
    (spectro, output)
}

impl<W, Q, P> CoordType<W, Q, P> for Cart
where
    W: io::Write,
    Q: Queue<P> + Sync,
    P: Program + Clone + Send + Sync,
{
    fn run(&self, w: &mut W, queue: &Q, config: &Config) -> (Spectro, Output) {
        let (n, nfc2, nfc3, mut fcs, mol, energies, mut target_map) =
            Cart.first_part(w, config, queue);

        time!(w, "freqs",
          let (fc2, f3, f4) = <Self as FiniteDifference<W, Q, P>>::make_fcs(
          self,
          &mut target_map,
          &energies,
          &mut fcs,
          n,
          Derivative::Quartic(
          nfc2,
          nfc3,
          0),
          "freqs",
          );

          let r = freqs("freqs", &mol, fc2, f3, f4);
        );
        r
    }
}

impl<W, Q, P> FiniteDifference<W, Q, P> for Cart
where
    W: io::Write,
    Q: Queue<P> + Sync,
    P: Program + Clone + Send + Sync,
{
}

impl Cart {
    /// run the "first part" of the Cartesian QFF, including the optimization if
    /// requested and the generation and running of the single-point energies
    pub(crate) fn first_part<W, Q, P>(
        &self,
        w: &mut W,
        config: &Config,
        queue: &Q,
    ) -> (usize, usize, usize, Vec<f64>, Molecule, Vec<f64>, BigHash)
    where
        W: io::Write,
        Q: Queue<P> + Sync,
        P: Program + Clone + Send + Sync,
    {
        time!(w, "opt",
            let template = Template::from(&config.template);
            let (geom, ref_energy) = if config.optimize {
                let res = optimize(
                    queue,
                    config.geometry.clone(),
                    template.clone(),
                    config.charge,
                )
                .expect("optimization failed");
                (Geom::Xyz(res.cart_geom.unwrap()), res.energy)
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
        let mut fcs = vec![0.0; nfc2 + nfc3 + nfc4];
        let mut mol = Molecule::new(geom.to_vec());
        mol.normalize();
        let pg = mol.point_group();
        writeln!(w, "normalized geometry:\n{}", mol).unwrap();
        writeln!(w, "point group:{}", pg).unwrap();
        let mut target_map = BigHash::new(mol.clone(), pg);
        time! (w, "building points",
               let geoms = <Self as FiniteDifference<W,Q, P>>::build_points(
           self,
               Geom::Xyz(mol.atoms.clone()),
               config.step_size,
               ref_energy,
               Derivative::Quartic(nfc2, nfc3, nfc4),
               &mut fcs,
               &mut target_map,
               );
        );
        let dir = "pts";
        let mut jobs = {
            let mut jobs = Vec::new();
            for (job_num, mol) in geoms.into_iter().enumerate() {
                let filename = format!("{dir}/job.{:08}", job_num);
                let job = Job::new(
                    P::new(filename, template.clone(), config.charge, mol.geom),
                    mol.index,
                );
                jobs.push(job);
            }
            jobs
        };
        writeln!(
            w,
            "{n} Cartesian coordinates requires {} points",
            jobs.len()
        )
        .unwrap();
        time!(w, "draining points",
              // drain into energies
              let mut energies = vec![0.0; jobs.len()];
              queue
              .drain(dir, &mut jobs, &mut energies)
              .expect("single-point calculations failed");
        );
        (n, nfc2, nfc3, fcs, mol, energies, target_map)
    }
}

extern crate test;
#[bench]
fn bench_to_keys(b: &mut test::Bencher) {
    use std::str::FromStr;
    let mol = Molecule::from_str(
        "
    C        0.000000   -0.888844    0.000000
    C       -0.662697    0.368254    0.000000
    C        0.662697    0.368254    0.000000
    H       -1.595193    0.906925    0.000000
    H        1.595193    0.906925    0.000000
",
    )
    .unwrap();
    b.iter(|| BigHash::to_keys(&mol))
}
