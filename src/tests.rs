use std::io::Stdout;

use approx::abs_diff_ne;
use approx::assert_abs_diff_eq;
use intder::Intder;
use nalgebra::dvector;
use psqs::geom::Geom;
use psqs::program::mopac::Mopac;
use psqs::program::Template;
use psqs::queue::local::LocalQueue;
use rust_anpass::Dvec;
use symm::Molecule;

use crate::cleanup;
use crate::config::Config;
use crate::coord_type::normal::Normal;
use crate::coord_type::BigHash;
use crate::coord_type::Cart;
use crate::coord_type::CoordType;
use crate::coord_type::SIC;
use crate::optimize;

#[test]
#[ignore]
fn normal() {
    cleanup();
    let _ = std::fs::create_dir("opt");
    let _ = std::fs::create_dir("pts");
    let _ = std::fs::create_dir("freqs");
    let config = Config::load("testfiles/water.toml");
    let queue = LocalQueue {
        dir: "pts".to_string(),
        chunk_size: 512,
        ..Default::default()
    };
    let (_, summ) = <Normal as CoordType<Stdout, LocalQueue, Mopac>>::run(
        &Normal,
        &mut std::io::stdout(),
        &queue,
        &config,
    );

    // these match the Go version from
    // ~/chem/c3h2/reparam_cart/16/qffs/000/freqs/spectro2.out on eland

    // harmonics
    assert_abs_diff_eq!(
        Dvec::from(summ.harms),
        dvector![2602.732027773963, 2523.471527899651, 1315.4979243007172],
        epsilon = 2e-3
    );
    let got = Dvec::from(summ.corrs);
    let want = dvector![2584.12640107, 2467.79949790, 1118.91300446];
    // corr
    if abs_diff_ne!(got, want, epsilon = 2.6e-1) {
        println!("got={:.8}", got);
        println!("want={:.8}", want);
        println!("diff={:.8}", got - want);
        panic!("mismatch");
    }
}

#[test]
#[ignore]
fn h2o_sic() {
    cleanup();
    let config = Config::load("testfiles/water.toml");
    let coord = SIC::new(Intder::load_file("testfiles/h2o.intder"));
    let queue = LocalQueue {
        dir: "pts".to_string(),
        chunk_size: 512,
        ..Default::default()
    };
    let (_, summ) = <SIC as CoordType<Stdout, LocalQueue, Mopac>>::run(
        &coord,
        &mut std::io::stdout(),
        &queue,
        &config,
    );

    // these match the Go version from
    // ~/chem/c3h2/reparam_cart/16/qffs/000/freqs/spectro2.out on eland

    // harmonics
    assert_abs_diff_eq!(
        Dvec::from(summ.harms),
        dvector![2603.972868873955, 2523.265599529732, 1315.3583799359571],
        epsilon = 2e-3
    );
    let got = Dvec::from(summ.corrs);
    let want = dvector![2636.08089954, 2437.94844727, 1283.99909982];
    // corr
    if abs_diff_ne!(got, want, epsilon = 2.6e-1) {
        println!("got={:.8}", got);
        println!("want={:.8}", want);
        println!("diff={:.8}", got - want);
        panic!("mismatch");
    }
}

#[test]
#[ignore]
fn sic() {
    cleanup();
    let config = Config::load("testfiles/test.toml");
    let coord = SIC::new(Intder::load_file("testfiles/intder.in"));
    let queue = LocalQueue {
        dir: "pts".to_string(),
        chunk_size: 512,
        ..Default::default()
    };
    let (_, summ) = <SIC as CoordType<Stdout, LocalQueue, Mopac>>::run(
        &coord,
        &mut std::io::stdout(),
        &queue,
        &config,
    );

    // these match the Go version from
    // ~/chem/c3h2/reparam_cart/16/qffs/000/freqs/spectro2.out on eland

    // harmonics
    assert_abs_diff_eq!(
        Dvec::from(summ.harms),
        dvector![
            2820.227, 2799.282, 1819.161, 1198.887, 1060.531, 963.513, 931.318,
            929.900, 912.358
        ],
        epsilon = 2e-3
    );
    let got = Dvec::from(summ.corrs);
    let want = Dvec::from(vec![
        2783.9552, 2764.3024, 1775.6603, 1177.1131, 1040.6267, 960.1012,
        919.9009, 926.9755, 905.3032,
    ]);
    println!("{:.8}", got.clone() - want.clone());
    // corr
    assert_abs_diff_eq!(got, want, epsilon = 2.6e-1);
}

#[test]
#[ignore]
fn cart() {
    cleanup();
    let _ = std::fs::create_dir("opt");
    let _ = std::fs::create_dir("pts");
    let _ = std::fs::create_dir("freqs");

    let config = Config::load("testfiles/cart.toml");
    let queue = LocalQueue {
        dir: "pts".to_string(),
        chunk_size: 512,
        ..Default::default()
    };
    let (_, summ) = <Cart as CoordType<Stdout, LocalQueue, Mopac>>::run(
        &Cart,
        &mut std::io::stdout(),
        &queue,
        &config,
    );
    assert_eq!(summ.harms.len(), 9);
    // harmonics
    assert_abs_diff_eq!(
        Dvec::from(summ.harms),
        dvector![
            2819.297, 2798.273, 1819.846, 1199.526, 1061.197, 964.357, 932.103,
            930.917, 913.221
        ],
        epsilon = 1.0
    );
    // corr
    assert_abs_diff_eq!(
        Dvec::from(summ.corrs),
        Dvec::from(vec![
            2783.1, 2763.3, 1776.4, 1177.8, 1041.3, 960.0, 920.7, 927.3, 906.1
        ]),
        epsilon = 1.0
    );
}

/// test approximately to check that I don't break the symmetry stuff without
/// running the whole QFF.
#[test]
fn build_pts() {
    let config = Config::load("testfiles/cart.toml");
    let queue = LocalQueue {
        dir: "pts".to_string(),
        chunk_size: 512,
        ..Default::default()
    };
    let (geom, ref_energy) = {
        let res = optimize::<LocalQueue, Mopac>(
            &queue,
            config.geometry.clone(),
            Template::from(&config.template),
            config.charge,
        )
        .expect("optimization failed in build_pts");
        (Geom::Xyz(res.cart_geom.unwrap()), res.energy)
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
    println!("normalized geometry:\n{}", mol);
    let mut target_map = BigHash::new(mol.clone(), pg);

    let geoms = Cart.build_points(
        Geom::Xyz(mol.atoms.clone()),
        config.step_size,
        ref_energy,
        nfc2,
        nfc3,
        &mut fcs,
        &mut target_map,
    );
    assert_eq!(geoms.len(), 11952);
}
