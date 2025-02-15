use super::load_dmat;
use crate::*;

#[derive(Clone)]
struct Test {
    infile: String,
    fort15: String,
    want: Dmat,
}

impl Test {
    fn new(dir: &'static str, n3n: usize) -> Self {
        let start = Path::new("testfiles");
        Self {
            infile: String::from(
                start.join(dir).join("spectro.in").to_str().unwrap(),
            ),
            fort15: String::from(
                start.join(dir).join("fort.15").to_str().unwrap(),
            ),
            want: load_dmat(
                start.join(dir).join("fx").to_str().unwrap(),
                n3n,
                n3n,
            ),
        }
    }
}

#[test]
fn rot2nd() {
    let tests = [
        Test::new("h2o", 9),
        Test::new("h2co", 12),
        Test::new("c3h2", 15),
        Test::new("c3hf", 15),
        Test::new("c3hcn", 18),
        Test::new("c3hcn010", 18),
        Test::new("c2h-", 9),
    ];
    for test in tests {
        let s = Spectro::load(&test.infile);
        let fc2 = load_fc2(test.fort15, s.n3n);
        let got = s.rot2nd(&fc2);
        check_mat!(&got.abs(), &test.want.abs(), 1.08e-10, test.infile);
    }
}
