use approx::abs_diff_ne;

use crate::{consts::FACT2, *};

#[derive(Clone)]
struct Test {
    infile: String,
    fort15: String,
    want: Dmat,
}

impl Test {
    fn new(dir: &'static str, n3n: usize) -> Self {
        let start = Path::new("testfiles");
        // this lets me load fxm straight from the fortran version
        let v = load_vec(start.join(dir).join("fxm").to_str().unwrap());
        let mut want = Dmat::zeros(n3n, n3n);
        let mut i = 0;
        let mut j = 0;
        for a in 0..v.len() {
            if j > i {
                j = 0;
                i += 1;
            }
            if i >= n3n {
                // should be done before getting here
                panic!("{}", a);
            }
            want[(i, j)] = v[a];
            want[(j, i)] = v[a];
            j += 1;
        }
        Self {
            infile: String::from(
                start.join(dir).join("spectro.in").to_str().unwrap(),
            ),
            fort15: String::from(
                start.join(dir).join("fort.15").to_str().unwrap(),
            ),
            want,
        }
    }
}

#[test]
pub(crate) fn test_sec() {
    let tests = [
        Test::new("h2o", 9),
        Test::new("h2co", 12),
        Test::new("c3h2", 15),
        Test::new("c3hf", 15),
        Test::new("c3hcn", 18),
        Test::new("c3hcn010", 18),
    ];
    for test in tests {
        let s = Spectro::load(test.infile);
        let fc2 = load_fc2(test.fort15, s.n3n);
        let fc2 = s.rot2nd(&fc2);
        let fc2 = FACT2 * fc2;
        let w = s.geom.weights();
        let sqm: Vec<_> = w.iter().map(|w| 1.0 / w.sqrt()).collect();
        let got = s.form_sec(fc2, &sqm);
        if abs_diff_ne!(got.abs(), test.want.abs(), epsilon = 1.73e-8) {
            println!("got\n{got:.8}");
            println!("want\n{:.8}", test.want);
            println!(
                "max diff = {:.2e}",
                (got.clone().abs() - test.want.clone().abs()).abs().max()
            );
        }
    }
}
