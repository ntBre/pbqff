use rust_anpass::Dvec;

use crate::run;

#[test]
fn full() {
    let summ = run(
        "testfiles/test.toml",
        "testfiles/intder.in",
        "testfiles/spectro.in",
    );
    approx::assert_abs_diff_eq!(
        Dvec::from(summ.corr),
        Dvec::from(vec![
            2784.0, 2764.3, 1775.7, 1177.1, 1040.6, 960.1, 920.0, 927.0, 905.3,
        ]),
        epsilon = 1.0
    )
}
