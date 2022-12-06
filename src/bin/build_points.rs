use psqs::{
    geom::Geom,
    program::{mopac::Mopac, Template},
    queue::local::LocalQueue,
};
use rust_pbqff::{
    config::Config,
    coord_type::{
        findiff::{bighash::BigHash, FiniteDifference},
        Cart, Derivative,
    },
    optimize,
};
use symm::Molecule;

fn main() {
    let config = Config::load("testfiles/cart.toml");
    let template = Template::from(&config.template);

    let queue = LocalQueue {
        chunk_size: 128,
        dir: "pts".to_string(),
        ..Default::default()
    };
    let (geom, ref_energy) = if config.optimize {
        let res = optimize::<LocalQueue, Mopac>(
            &queue,
            config.geometry.clone(),
            template,
            config.charge,
        )
        .expect("optimization failed");
        (Geom::Xyz(res.cart_geom.unwrap()), res.energy)
    } else {
        todo!()
    };

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
    let mut target_map = BigHash::new(mol.clone(), pg);

    Cart.build_points(
        Geom::Xyz(mol.atoms.clone()),
        config.step_size,
        ref_energy,
        Derivative::Quartic(nfc2, nfc3, nfc4),
        &mut fcs,
        &mut target_map,
    );
}
