//! automatic normal coordinates. the normal coordinates are obtained by
//! running a harmonic, cartesian force field and feeding that to spectro to
//! obtain the LXM matrix. then the cartesian QFF machinery is used to
//! generate the normal coordinate displacements, which can be fed back in
//! to spectro at the end as f3qcm and f4qcm

use std::{io::Write, marker::Sync};

use nalgebra::DVector;
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
    /// the normal coordinates, called the LXM matrix in spectro
    lxm: Option<nalgebra::DMatrix<f64>>,

    /// 1/√m where m is the atomic mass
    m12: Vec<f64>,

    /// the number of normal coordinates. convenient to have here so I don't
    /// have to keep thinking about 3n-6/5 stuff
    ncoords: usize,
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

        self.map_energies(&map, &energies, &mut fcs);

        let cubs = &fcs[nfc2..nfc2 + nfc3];
        let quarts = &fcs[nfc2 + nfc3..];

        let freq = &o.harms;
        let mut ijk = 0;
        let mut ijkl = 0;
        let mut f3qcm = Vec::new();
        let mut f4qcm = Vec::new();
        for i in 0..n {
            let wi = freq[(i)];
            for j in 0..=i {
                let wj = freq[(j)];
                for k in 0..=j {
                    let wk = freq[(k)];
                    let wijk = wi * wj * wk;
                    let fact =
                        intder::HART * spectro::consts::FACT3 / wijk.sqrt();
                    let val = cubs[ijk];
                    f3qcm.push(val * fact);
                    ijk += 1;
                    for l in 0..=k {
                        let wl = freq[l];
                        let wijkl = wijk * wl;
                        let fact = intder::HART * spectro::consts::FACT4
                            / wijkl.sqrt();
                        let val = quarts[ijkl];
                        f4qcm.push(val * fact);
                        ijkl += 1;
                    }
                }
            }
        }

        let (o, _) = s.finish(
            DVector::from(o.harms.clone()),
            spectro::F3qcm::new(f3qcm),
            spectro::F4qcm::new(f4qcm),
            o.irreps,
            self.lxm.unwrap(),
        );

        (s, o)

        // this is f3qcm from the cart run:
        //  519.79254777
        //    0.00000000
        //   42.91584633
        //    0.00000000
        //  323.20772069
        //   -0.00000000
        //  102.27900299
        // -544.13505623
        //    0.00000000
        //  710.18634048

        // and these are f4qcm from the cart run
        //  279.69429451
        //   -0.00000000
        //  -97.10265238
        //   -0.00000000
        // -147.94627196
        //   86.48307573
        //    0.00000000
        //   59.85703101
        //    0.00000000
        // -166.94324978
        //    0.00000000
        // -194.08270483
        //   87.87581824
        //    0.00000000
        // -237.54261635
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
