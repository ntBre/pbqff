use std::fs::File;
use std::io::Write;

use intder::Intder;

fn parse_args<I: Iterator<Item = String>>(args: &mut I) -> Option<String> {
    let mut hold = Vec::new();
    for arg in args {
        if arg == "-v" {
            unsafe {
                std::env::set_var("INTDER_DEBUG", "1");
            }
        } else {
            hold.push(arg);
        }
    }
    hold.get(1).map(String::from)
}

#[test]
fn test_parse_args() {
    assert_eq!(
        parse_args(
            &mut vec!["intder", "-v", "intder.in"]
                .into_iter()
                .map(|s| s.to_string())
        ),
        Some("intder.in".to_string())
    );

    assert_eq!(
        parse_args(
            &mut vec!["intder", "intder.in"]
                .into_iter()
                .map(|s| s.to_string())
        ),
        Some("intder.in".to_string())
    );

    assert_eq!(
        parse_args(&mut vec!["intder"].into_iter().map(|s| s.to_string())),
        None,
    );
}

fn main() {
    let mut args = std::env::args();
    let infile = parse_args(&mut args);
    let intder = match infile {
        Some(s) => Intder::load_file(&s),
        None => Intder::load(std::io::stdin()),
    };
    if intder.input_options[14] != 0 {
        let new_carts = intder.convert_disps().unwrap();
        let mut file07 =
            File::create("file07").expect("failed to create file07");
        for cart in new_carts {
            writeln!(file07, "# GEOMUP #################").unwrap();
            Intder::print_cart(&mut file07, &cart);
        }
    } else {
        let (f2, f3, f4) = intder.convert_fcs();
        Intder::dump_fcs(".", &f2, &f3, &f4);
    }
}
