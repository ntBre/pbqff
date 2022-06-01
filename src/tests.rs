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
    approx::assert_abs_diff_eq!(
        Dvec::from(summ.corr),
        Dvec::from(vec![
            2784.0, 2764.3, 1775.7, 1177.1, 1040.6, 960.1, 920.0, 927.0, 905.3,
        ]),
        epsilon = 1.0
    )
}

#[test]
#[ignore]
fn cart() {
    let config = Config::load("testfiles/test.toml");
    let spectro = Spectro::load("testfiles/spectro.in");
    let summ = Cart.run(&config, &spectro);
    approx::assert_abs_diff_eq!(
        Dvec::from(summ.corr),
	// this is the result from Go on eland
	// 2783.1, 2763.3, 1776.4, 1177.8, 1041.3, 960.0, 920.7, 927.3, 906.1

	// this is the result on cactus
        Dvec::from(vec![
	    2790.7, 2771.1, 1777.4, 1179.2, 1043.6, 967.7, 923.6, 934.9, 909.7
        ]),
        epsilon = 1.0
    )
}
