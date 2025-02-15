use criterion::{criterion_group, criterion_main, Criterion};
use psqs::program::{mopac::Mopac, Program, Template};
use symm::molecule;

pub fn write_input(c: &mut Criterion) {
    let mol = molecule![
        C      0.000000000000      0.003768239200     -1.686245109400
        C      0.000000000000      1.243805099800      0.688097726900
        C      0.000000000000     -1.242134139500      0.700109341300
        H      0.000000000000      2.998368900600      1.719180578800
        H      0.000000000000     -3.003808100200      1.718996398400
    ];
    let mut mop = Mopac::new(
        "/tmp/job".to_owned(),
        Template::from("scfcrt=1.D-21 aux(precision=14) PM6 SINGLET THREADS=1"),
        0,
        psqs::geom::Geom::Xyz(mol.atoms),
    );

    c.bench_function("write_input", |b| {
        b.iter(|| mop.write_input(psqs::program::Procedure::SinglePt));
    });
}

criterion_group!(benches, write_input);
criterion_main!(benches);
