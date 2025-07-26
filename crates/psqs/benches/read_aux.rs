use criterion::{Criterion, criterion_group, criterion_main};
use psqs::program::mopac::Mopac;

pub fn read_aux(c: &mut Criterion) {
    c.bench_function("read aux", |b| {
        b.iter(|| Mopac::read_aux("testfiles/job"))
    });
}

criterion_group!(benches, read_aux);
criterion_main!(benches);
