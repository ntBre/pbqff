use intder::Intder;
use rust_anpass::Dvec;
use spectro::Spectro;

use crate::config::Config;
use crate::coord_type::Cart;
use crate::coord_type::CoordType;
use crate::coord_type::SIC;

#[test]
#[ignore]
fn sic() {
    let config = Config::load("testfiles/test.toml");
    let coord = SIC::new(Intder::load_file("testfiles/intder.in"));
    let spectro = Spectro::load("testfiles/spectro.in");
    let summ = coord.run(&config, &spectro);

    // harmonics
    approx::assert_abs_diff_eq!(
        Dvec::from(summ.harm),
        Dvec::from(vec![
            2820.2, 2799.3, 1819.2, 1198.9, 1060.5, 963.5, 931.3, 929.9, 912.4,
        ]),
        epsilon = 1.0
    );
    // corr
    approx::assert_abs_diff_eq!(
        Dvec::from(summ.corr),
        Dvec::from(vec![
            2784.0, 2764.3, 1775.7, 1177.1, 1040.6, 960.1, 920.0, 927.0, 905.3,
        ]),
        epsilon = 1.0
    );
}

#[test]
#[ignore]
fn cart() {
    let _ = std::fs::remove_dir_all("opt");
    let _ = std::fs::remove_dir_all("pts");
    let _ = std::fs::remove_dir_all("freqs");

    let _ = std::fs::create_dir("opt");
    let _ = std::fs::create_dir("pts");
    let _ = std::fs::create_dir("freqs");

    let config = Config::load("testfiles/test.toml");
    let spectro = Spectro::load("testfiles/spectro.in");
    let summ = Cart.run(&config, &spectro);
    assert_eq!(summ.harm.len(), 9);
    // harmonics
    approx::assert_abs_diff_eq!(
        Dvec::from(summ.harm),
        Dvec::from(vec![
            2819.297, 2798.273, 1819.846, 1199.526, 1061.197, 964.357, 932.103,
            930.917, 913.221,
        ]),
        epsilon = 1e-3
    );
    // corr
    approx::assert_abs_diff_eq!(
        Dvec::from(summ.corr),
        // this is the result from Go on eland
        // 2783.1, 2763.3, 1776.4, 1177.8, 1041.3, 960.0, 920.7, 927.3, 906.1

        // this is the result on cactus
        Dvec::from(vec![
            2791.6, 2773.8, 1769.7, 1172.5, 1049.2, 971.8, 931.1, 934.4, 908.3,
        ]),
        epsilon = 1.0
    );
}
