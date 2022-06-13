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
    let summ = coord.run(&mut std::io::stdout(), &config, &spectro);

    // these match the Go version from
    // ~/chem/c3h2/reparam_cart/16/qffs/000/freqs/spectro2.out on eland

    // harmonics
    approx::assert_abs_diff_eq!(
        Dvec::from(summ.harm),
        Dvec::from(vec![
            2820.2, 2799.3, 1819.2, 1198.9, 1060.5, 963.5, 931.3, 929.9, 912.4,
        ]),
        epsilon = 0.1
    );
    // corr
    approx::assert_abs_diff_eq!(
        Dvec::from(summ.corr),
        Dvec::from(vec![
            2784.0, 2764.3, 1775.7, 1177.1, 1040.6, 960.1, 920.0, 927.0, 905.3,
        ]),
        epsilon = 0.1
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
    let summ = Cart.run(&mut std::io::stdout(), &config, &spectro);
    assert_eq!(summ.harm.len(), 9);
    // harmonics
    approx::assert_abs_diff_eq!(
        Dvec::from(summ.harm),
        Dvec::from(vec![
            2819.297, 2798.273, 1819.846, 1199.526, 1061.197, 964.357, 932.103,
            930.917, 913.221,
        ]),
        epsilon = 1.0
    );
    // corr
    approx::assert_abs_diff_eq!(
        Dvec::from(summ.corr),
        Dvec::from(vec![
            2783.1, 2763.3, 1776.4, 1177.8, 1041.3, 960.0, 920.7, 927.3, 906.1
        ]),
        epsilon = 1.0
    );

    // TODO test on this eventually: looks good for now though

    // A_e   34262.7
    // B_e   31847.3
    // C_e   16505.4
    // A_0   34199.6
    // B_0   31812.0
    // C_0   16430.8
    // A_1   34195.4
    // B_1   31789.3
    // C_1   16424.0
    // A_2   34196.4
    // B_2   31796.8
    // C_2   16426.2
    // A_3   34167.7
    // B_3   31685.5
    // C_3   16390.6
    // A_4   34319.1
    // B_4   31767.2
    // C_4   16427.6
    // A_5   34054.7
    // B_5   31940.6
    // C_5   16390.0
    // A_6   34205.9
    // B_6   31725.9
    // C_6   16450.8
    // A_7   34170.5
    // B_7   31904.4
    // C_7   16413.8
}
