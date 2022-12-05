//! automatic normal coordinates. the normal coordinates are obtained by
//! running a harmonic, cartesian force field and feeding that to spectro to
//! obtain the LXM matrix. then the cartesian QFF machinery is used to
//! generate the normal coordinate displacements, which can be fed back in
//! to spectro at the end as f3qcm and f4qcm

use std::io::Write;

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
        let r = cart_part(config, queue, w);
        println!("r.1.lxm={:#?}", r.1.lxm);

        // TODO at this point we have what we need for the normal coordinate
        // displacements. I need to write Normal.build_points, trying to reuse
        // as much from the Cart implementation as possible. I especially want
        // to avoid duplicating all of the finite difference code. it looks like
        // the main difference is the new_geom implementation, which just needs
        // to accept LXM. maybe lxm should be a field on Normal and then in the
        // macros I can call self.new_geom instead of just new_geom. that has
        // the slight downside of having to pass self to the macros, but it
        // might be the best bet

        r
    }
}

/// run the Cartesian harmonic force field and return the spectro output, from
/// which we can extract the geometry and normal coordinates (lxm)
fn cart_part<
    P: Program + Clone + Send + std::marker::Sync,
    Q: psqs::queue::Queue<P> + std::marker::Sync,
    W: Write,
>(
    config: &Config,
    queue: &Q,
    w: &mut W,
) -> (Spectro, Output) {
    // TODO call first part of run from carts
    // harm_freqs("freqs", &mol, fc2):
    todo!()
}
