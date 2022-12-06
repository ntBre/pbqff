//! automatic normal coordinates. the normal coordinates are obtained by
//! running a harmonic, cartesian force field and feeding that to spectro to
//! obtain the LXM matrix. then the cartesian QFF machinery is used to
//! generate the normal coordinate displacements, which can be fed back in
//! to spectro at the end as f3qcm and f4qcm

use std::{io::Write, marker::Sync};

use psqs::{program::Program, queue::Queue};
use spectro::{Output, Spectro};
use symm::Molecule;

use crate::config::Config;

use super::{
    findiff::FiniteDifference, Cart, CoordType, Derivative, Nderiv,
    SPECTRO_HEADER,
};

pub struct Normal;

impl<W, Q, P> CoordType<W, Q, P> for Normal
where
    W: Write,
    Q: Queue<P> + Sync,
    P: Program + Clone + Send + Sync,
{
    fn run(&self, w: &mut W, queue: &Q, config: &Config) -> (Spectro, Output) {
        let (r, o) = self.cart_part(config, queue, w);
        // pretty sure I assert lxm is square somewhere in spectro though. lxm
        // should be in column-major order so I think this is all right
        let cols = o.lxm.len();
        let rows = o.lxm[0].len();
        let lxm = nalgebra::DMatrix::from_iterator(
            rows,
            cols,
            o.lxm.iter().flatten().cloned(),
        );

        println!("lxm={:.8}", lxm);

        // TODO at this point we have what we need for the normal coordinate
        // displacements. I need to write Normal.build_points, trying to reuse
        // as much from the Cart implementation as possible. I especially want
        // to avoid duplicating all of the finite difference code. it looks like
        // the main difference is the new_geom implementation, which just needs
        // to accept LXM. maybe lxm should be a field on Normal and then in the
        // macros I can call self.new_geom instead of just new_geom. that has
        // the slight downside of having to pass self to the macros, but it
        // might be the best bet

        // Cart and Normal share a whole lot in common. I think I should define
        // a trait like FinDiff to represent that commonality

        (r, o)
    }
}

impl FiniteDifference for Normal {
    fn new_geom(
        &self,
        _names: &[&str],
        _coords: nalgebra::DVector<f64>,
        _step_size: f64,
        _steps: Vec<isize>,
    ) -> psqs::geom::Geom {
        todo!()
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
    ) -> (Spectro, Output)
    where
        P: Program + Clone + Send + Sync,
        Q: Queue<P> + Sync,
        W: Write,
    {
        let (n, nfc2, nfc3, mut fcs, mol, energies, mut target_map) =
            Cart.first_part(w, config, queue, Nderiv::Two);
        let (fc2, _, _) = self.make_fcs(
            &mut target_map,
            &energies,
            &mut fcs,
            n,
            Derivative::Quartic(nfc2, nfc3, 0),
            "freqs",
        );
        self.harm_freqs("freqs", &mol, fc2)
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
