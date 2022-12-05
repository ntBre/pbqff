//! automatic normal coordinates. the normal coordinates are obtained by
//! running a harmonic, cartesian force field and feeding that to spectro to
//! obtain the LXM matrix. then the cartesian QFF machinery is used to
//! generate the normal coordinate displacements, which can be fed back in
//! to spectro at the end as f3qcm and f4qcm

use std::{io::Write, time::Instant};

use psqs::{
    geom::Geom,
    program::{Job, Program, Template},
    queue::Queue,
};
use spectro::{Output, Spectro};
use symm::Molecule;

use crate::{config::Config, optimize, ref_energy};

use super::{make_fcs, BigHash, Cart, CoordType, SPECTRO_HEADER};

pub struct Normal;

pub fn harm_freqs(
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

impl<
        W: Write,
        Q: Queue<P> + Sync,
        P: Program + Clone + Send + std::marker::Sync,
    > CoordType<W, Q, P> for Normal
{
    fn run(&self, w: &mut W, queue: &Q, config: &Config) -> (Spectro, Output) {
        let mut now = Instant::now();

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

        writeln!(
            w,
            "finished opt after {:.1} sec",
            now.elapsed().as_millis() as f64 / 1000.0
        )
        .unwrap();
        now = Instant::now();

        let geom = geom.xyz().expect("expected an XYZ geometry, not Zmat");
        // 3 * #atoms
        let n = 3 * geom.len();
        let nfc2 = n * n;
        let mut fcs = vec![0.0; nfc2];

        let mut mol = Molecule::new(geom.to_vec());
        mol.normalize();
        let pg = mol.point_group();
        writeln!(w, "normalized geometry:\n{}", mol).unwrap();
        writeln!(w, "point group:{}", pg).unwrap();
        let mut target_map = BigHash::new(mol.clone(), pg);

        let geoms = Cart.build_points(
            Geom::Xyz(mol.atoms.clone()),
            config.step_size,
            ref_energy,
            super::Derivative::Harmonic(nfc2),
            &mut fcs,
            &mut target_map,
        );

        writeln!(
            w,
            "finished building points after {:.1} sec",
            now.elapsed().as_millis() as f64 / 1000.0
        )
        .unwrap();
        now = Instant::now();

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

        // drain into energies
        let mut energies = vec![0.0; jobs.len()];
        queue
            .drain(dir, &mut jobs, &mut energies)
            .expect("single-point calculations failed");

        writeln!(
            w,
            "finished draining points after {:.1} sec",
            now.elapsed().as_millis() as f64 / 1000.0
        )
        .unwrap();
        now = Instant::now();

        let (fc2, _, _) =
            make_fcs(&mut target_map, &energies, &mut fcs, n, nfc2, 0, "freqs");

        let r = harm_freqs("freqs", &mol, fc2);
        println!("r.1.lxm={:#?}", r.1.lxm);
        writeln!(
            w,
            "finished freqs after {:.1} sec",
            now.elapsed().as_millis() as f64 / 1000.0
        )
        .unwrap();
        r
    }
}
