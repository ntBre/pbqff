use symm::Molecule;

use crate::coord_type::findiff::bighash::BigHash;

extern crate test;
#[bench]
fn bench_to_keys(b: &mut test::Bencher) {
    use std::str::FromStr;
    let mol = Molecule::from_str(
        "
    C        0.000000   -0.888844    0.000000
    C       -0.662697    0.368254    0.000000
    C        0.662697    0.368254    0.000000
    H       -1.595193    0.906925    0.000000
    H        1.595193    0.906925    0.000000
",
    )
    .unwrap();
    b.iter(|| BigHash::to_keys(&mol))
}
