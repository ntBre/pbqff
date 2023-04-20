use std::fs;
use std::io::Stdout;
use std::str::FromStr;

use approx::abs_diff_ne;
use approx::assert_abs_diff_eq;
use intder::Intder;
use nalgebra::dvector;
use psqs::geom::Geom;
use psqs::program::mopac::Mopac;
use psqs::queue::local::Local;
use rust_anpass::Dvec;
use spectro::Output;
use spectro::Spectro;
use symm::Molecule;

use crate::cleanup;
use crate::config::Config;
use crate::coord_type::findiff::bighash::BigHash;
use crate::coord_type::findiff::FiniteDifference;
use crate::coord_type::normal;
use crate::coord_type::normal::DerivType;
use crate::coord_type::normal::Normal;
use crate::coord_type::normal::Resume;
use crate::coord_type::Cart;
use crate::coord_type::CoordType;
use crate::coord_type::Load;
use crate::coord_type::Sic;

#[test]
#[ignore]
fn h2o_normal() {
    cleanup();
    init();
    let config = Config::load("testfiles/water.toml");
    let queue = Local {
        dir: "pts".to_string(),
        chunk_size: 512,
        ..Default::default()
    };
    let (_, summ) = <Normal as CoordType<Stdout, Local, Mopac>>::run(
        Normal::default(),
        &mut std::io::stdout(),
        &queue,
        &config,
    );

    // harmonics
    assert_abs_diff_eq!(
        Dvec::from(summ.harms),
        dvector![2602.732027773963, 2523.471527899651, 1315.4979243007172],
        epsilon = 2e-3
    );
    let got = Dvec::from(summ.corrs);
    let want = dvector![2584.12640107, 2467.79949790, 1118.91300446];
    assert_eq!(got.len(), want.len());
    // corr
    if abs_diff_ne!(got, want, epsilon = 2.6e-1) {
        println!("got={got:.8}");
        println!("want={want:.8}");
        println!("diff={:.8}", got - want);
        panic!("mismatch");
    }
}

macro_rules! check {
    ($got:expr, $want:expr, $eps:expr) => {
        if abs_diff_ne!($got, $want, epsilon = $eps) {
            println!("got={:.8}", $got);
            println!("want={:.8}", $want);
            println!("diff={:.8}", &$got - &$want);
            println!("max diff={:.2e}", ($got - $want).abs().max());
            panic!("mismatch at {}", line!());
        }
    };
}

#[test]
#[ignore]
fn c3h2_normal() {
    cleanup();
    init();
    let config = Config {
        check_int: 1,
        ..Config::load("testfiles/cart.toml")
    };
    let queue = Local {
        dir: "pts".to_string(),
        chunk_size: 512,
        ..Default::default()
    };
    let (_, summ) = <Normal as CoordType<Stdout, Local, Mopac>>::run(
        Normal::findiff(false),
        &mut std::io::stdout(),
        &queue,
        &config,
    );
    // harmonics
    let got = Dvec::from(summ.harms);
    let want = dvector![
        2819.297, 2798.273, 1819.846, 1199.526, 1061.197, 964.357, 932.103,
        930.917, 913.221
    ];
    // higher eps because comparing to the pure cart wants
    check!(got, want, 1e-1);
    let got = Dvec::from(summ.corrs);
    let want = dvector![
        2783.1, 2763.3, 1776.4, 1177.8, 1041.3, 960.0, 920.7, 927.3, 906.1
    ];
    assert_eq!(got.len(), want.len());
    // corr
    check!(got, want, 2.6e-1);

    let (_, summ) = <Normal as CoordType<Stdout, Local, Mopac>>::resume(
        Normal::findiff(false),
        &mut std::io::stdout(),
        &queue,
        &config,
        normal::Resume::load("res.chk"),
    );
    // harmonics
    let got = Dvec::from(summ.harms);
    let want = dvector![
        2819.297, 2798.273, 1819.846, 1199.526, 1061.197, 964.357, 932.103,
        930.917, 913.221
    ];
    // higher eps because comparing to the pure cart wants
    check!(got, want, 1e-1);
    let got = Dvec::from(summ.corrs);
    let want = dvector![
        2783.1, 2763.3, 1776.4, 1177.8, 1041.3, 960.0, 920.7, 927.3, 906.1
    ];
    assert_eq!(got.len(), want.len());
    // corr
    check!(got, want, 2.6e-1);
}

#[test]
#[ignore]
fn c3h2_normal_findiff() {
    cleanup();
    init();
    let config = Config {
        check_int: 1,
        ..Config::load("testfiles/cart.toml")
    };
    let queue = Local {
        dir: "pts".to_string(),
        chunk_size: 512,
        ..Default::default()
    };
    let (_, summ) = <Normal as CoordType<Stdout, Local, Mopac>>::run(
        Normal::findiff(true),
        &mut std::io::stdout(),
        &queue,
        &config,
    );
    // harmonics
    let got = Dvec::from(summ.harms);
    let want = dvector![
        2819.297, 2798.273, 1819.846, 1199.526, 1061.197, 964.357, 932.103,
        930.917, 913.221
    ];
    // higher eps because comparing to the pure cart wants
    check!(got, want, 1e-1);
    let got = Dvec::from(summ.corrs);
    let want = dvector![
        2783.1, 2763.3, 1776.4, 1177.8, 1041.3, 960.0, 920.7, 927.3, 906.1
    ];
    assert_eq!(got.len(), want.len());
    // corr
    check!(got, want, 2.6e-1);
}

#[test]
#[ignore]
fn h2o_cart() {
    cleanup();
    init();
    let config = Config::load("testfiles/water.toml");
    let queue = Local {
        dir: "pts".to_string(),
        chunk_size: 512,
        ..Default::default()
    };
    let (_, summ) = <Cart as CoordType<Stdout, Local, Mopac>>::run(
        Cart,
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
        println!("got={got:.8}");
        println!("want={want:.8}");
        println!("diff={:.8}", got - want);
        panic!("mismatch");
    }
}

#[test]
#[ignore]
fn h2o_sic() {
    cleanup();
    let config = Config::load("testfiles/water.toml");
    let coord = Sic::new(Intder::load_file("testfiles/h2o.intder"));
    let queue = Local {
        dir: "pts".to_string(),
        chunk_size: 512,
        ..Default::default()
    };
    let (_, summ) = <Sic as CoordType<Stdout, Local, Mopac>>::run(
        coord,
        &mut std::io::stdout(),
        &queue,
        &config,
    );

    // these match the Go version from
    // ~/chem/c3h2/reparam_cart/16/qffs/000/freqs/spectro2.out on eland

    // harmonics
    assert_abs_diff_eq!(
        Dvec::from(summ.harms),
        dvector![2613.383565566701, 2526.510364189273, 1334.5036927050803],
        epsilon = 2e-3
    );
    let got = Dvec::from(summ.corrs);
    let want = dvector![2601.05343567, 2473.32961982, 1157.98123703];
    // corr
    if abs_diff_ne!(got, want, epsilon = 2.6e-1) {
        println!("got={got:.8}");
        println!("want={want:.8}");
        println!("diff={:.8}", got - want);
        panic!("mismatch");
    }
}

#[test]
#[ignore]
fn sic() {
    cleanup();
    let config = Config::load("testfiles/test.toml");
    let coord = Sic::new(Intder::load_file("testfiles/intder.in"));
    let queue = Local {
        dir: "pts".to_string(),
        chunk_size: 512,
        ..Default::default()
    };
    let (_, summ) = <Sic as CoordType<Stdout, Local, Mopac>>::run(
        coord,
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

fn init() {
    psqs::max_threads(8);
    let _ = std::fs::create_dir("opt");
    let _ = std::fs::create_dir("pts");
    let _ = std::fs::create_dir("freqs");
}

#[test]
#[ignore]
fn cart() {
    cleanup();
    init();
    let config = Config::load("testfiles/cart.toml");
    let queue = Local {
        dir: "pts".to_string(),
        chunk_size: 512,
        ..Default::default()
    };
    let (_, summ) = <Cart as CoordType<Stdout, Local, Mopac>>::run(
        Cart,
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
/// running the whole QFF. While we're at it, make sure that Resume is
/// Serializable and Deserializable
#[test]
fn build_pts() {
    let config = Config::load("testfiles/cart.toml");
    let geom = Molecule::from_str(
        "
C       0.0000000000   0.0000000000   0.0000000000
C       1.4361996710   0.0000000000   0.0000000000
C       0.7993316434   1.1932051034   0.0000000000
H       2.3607105070  -0.5060383826   0.0000000000
H       0.8934572640   2.2429362586  -0.0000000000",
    )
    .unwrap()
    .atoms;
    let ref_energy = 0.0;

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

    let geoms = Cart.build_points(
        Geom::Xyz(mol.atoms.clone()),
        config.step_size,
        ref_energy,
        crate::coord_type::Derivative::Quartic(nfc2, nfc3, nfc4),
        &mut fcs,
        &mut target_map,
        n,
    );
    assert_eq!(geoms.len(), 11952);

    let resume = Resume::new(
        Normal::findiff(true),
        geoms.len(),
        Output::default(),
        Spectro::default(),
        DerivType::Findiff {
            map: target_map,
            fcs,
            n,
            nfc2,
            nfc3,
        },
    );
    let chk = "/tmp/build_pts.json";
    resume.dump(chk);
    let got = Resume::load(chk);
    if got != resume {
        assert_eq!(got.njobs, resume.njobs);
        assert_eq!(got.output, resume.output);
        assert_eq!(got.spectro, resume.spectro);
        got.dump("/tmp/got.json");
        let DerivType::Findiff {
	map: gmap, fcs: gfcs, n: gn, nfc2: gnfc2, nfc3: gnfc3
    } = got.deriv else { unreachable!() };
        let DerivType::Findiff {
	map: wmap, fcs: wfcs, n: wn, nfc2: wnfc2, nfc3: wnfc3
    } = resume.deriv else { unreachable!() };
        assert_eq!(gn, wn);
        assert_eq!(gnfc2, wnfc2);
        assert_eq!(gnfc3, wnfc3);
        assert_eq!(gfcs, wfcs);
        assert_eq!(gmap.pg, wmap.pg);
        assert_eq!(gmap.buddy, wmap.buddy);
        // asserting these as eq is very ugly when false, so check that the len
        // is equal, that all of the keys in gmap are contained in wmap and vice
        // versa, and then check that all of the values are the same to 6e-11
        // tolerance since there is some deviation, surprisingly. 5e-11 fails,
        // so this is pretty tight
        assert_eq!(gmap.map.len(), wmap.map.len());
        assert!(gmap.map.keys().all(|k| wmap.map.contains_key(k)));
        assert!(wmap.map.keys().all(|k| gmap.map.contains_key(k)));
        for (k, v) in gmap.map {
            assert_abs_diff_eq!(wmap.map.get(&k).unwrap(), &v, epsilon = 6e-11);
        }
    }
    fs::remove_file(chk).unwrap();
}
