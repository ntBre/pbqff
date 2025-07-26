use pbqff::{
    config::Config,
    coord_type::{
        Cart, Derivative,
        findiff::{FiniteDifference, bighash::BigHash},
    },
    optimize,
};
use psqs::{
    geom::Geom,
    program::{Template, mopac::Mopac},
    queue::local::Local,
};
use symm::Molecule;

fn main() {
    let config = Config::load("testfiles/cart.toml");
    let template = Template::from(&config.template);

    let queue = Local {
        chunk_size: 128,
        dir: "pts".to_string(),
        ..Default::default()
    };
    let (geom, ref_energy) = if config.optimize {
        let res = optimize::<Local, Mopac>(
            ".",
            &queue,
            config.geometry.clone(),
            template,
            config.charge,
        )
        .expect("optimization failed");
        (Geom::Xyz(res.cart_geom.unwrap()), res.energy)
    } else {
        unimplemented!()
    };

    let geom = geom.xyz().expect("expected an XYZ geometry, not Zmat");
    // 3 * #atoms
    let n = 3 * geom.len();
    let deriv = Derivative::quartic(n);
    let mut fcs = vec![0.0; deriv.nfcs()];

    let mut mol = Molecule::new(geom.to_vec());
    mol.normalize();
    let pg = mol.point_group();
    let mut target_map = BigHash::new(mol.clone(), pg);

    Cart.build_points(
        Geom::Xyz(mol.atoms.clone()),
        config.step_size,
        ref_energy,
        deriv,
        &mut fcs,
        &mut target_map,
        n,
    );
}
