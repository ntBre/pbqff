//! automatic normal coordinates. the normal coordinates are obtained by
//! running a harmonic, cartesian force field and feeding that to spectro to
//! obtain the LXM matrix. then the cartesian QFF machinery is used to
//! generate the normal coordinate displacements, which can be fed back in
//! to spectro at the end as f3qcm and f4qcm

use std::{io::Write, marker::Sync};

use psqs::{
    geom::Geom,
    program::{Job, Program, Template},
    queue::Queue,
};
use spectro::{Output, Spectro};
use symm::{Molecule, PointGroup};

use crate::config::Config;

use super::{
    findiff::{bighash::BigHash, zip_atoms, FiniteDifference},
    Cart, CoordType, Derivative, Nderiv, SPECTRO_HEADER,
};

#[derive(Default)]
pub struct Normal {
    lxm: Option<nalgebra::DMatrix<f64>>,
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
        let (s, o, ref_energy, pg) = self.cart_part(config, queue, w);
        // pretty sure I assert lxm is square somewhere in spectro though. lxm
        // should be in column-major order so I think this is all right
        let cols = o.lxm.len();
        let rows = o.lxm[0].len();
        self.lxm = Some(nalgebra::DMatrix::from_iterator(
            rows,
            cols,
            o.lxm.iter().flatten().cloned(),
        ));

        // 3n - 6 + 1 if linear = 3n - 5
        let n = 3 * o.geom.atoms.len() - 6 + s.rotor.is_linear() as usize;
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
        let template = Template::from(&config.template);
        let jobs: Vec<_> = geoms
            .into_iter()
            .enumerate()
            .map(|(job_num, mol)| {
                let filename = format!("{dir}/job.{:08}", job_num);
                Job::new(
                    P::new(filename, template.clone(), config.charge, mol.geom),
                    mol.index,
                )
            })
            .collect();
        writeln!(w, "{n} normal coordinates requires {} points", jobs.len())
            .unwrap();
        time!(w, "draining points",
              // drain into energies
              let mut energies = vec![0.0; jobs.len()];
              queue
              .drain(dir, jobs, &mut energies)
              .expect("single-point calculations failed");
        );
        dbg!(energies);
        (s, o)
    }
}

impl FiniteDifference for Normal {
    fn new_geom(
        &self,
        names: &[&str],
        coords: nalgebra::DVector<f64>,
        step_size: f64,
        steps: Vec<isize>,
    ) -> psqs::geom::Geom {
        let mut v = vec![0.0; coords.len()];
        for step in steps {
            if step < 1 {
                v[(-step - 1) as usize] -= step_size;
            } else {
                v[(step - 1) as usize] += step_size;
            }
        }
        let coords = coords + nalgebra::DVector::from(v);
        Geom::Xyz(zip_atoms(names, coords))
    }
}

impl Normal {
    /// run the Cartesian harmonic force field and return the spectro output, from
    /// which we can extract the geometry and normal coordinates (lxm)
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
